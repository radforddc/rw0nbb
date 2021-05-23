/*
   eventprocess.c

   event processing code for MJD
   - up to 16 digitizer modules for the Ge detectors, plus
   - QDCs etc for the veto

   David Radford   Nov 2016
*/
/* Process the global events that have been built in eventbuild()
   - extract energy, PSA cuts, etc
   - create histograms
   returns: -1 on error,
             1 on pulser event,
             2 on bad event (bad = cannot be processed; not dirty)
             0 otherwise (good event, can be clean or dirty)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define DEBUG        0
#define VERBOSE      0
#define EVENTLIST    1800  // min energy in keV for output to event list, or 0 for no list
#define EVENTLISTMAX 2350  // max energy in keV for output to event list
#define GRAN_CUT     1     // cut events due to granularity? (set to 0 for calibrations etc)
#define WRITE_SIGS   1     // write selected signals to an otput file
#define DBIT(j)    (1 - (1&(dirty_sig>>j)))


/*  ------------------------------------------------------------ */

int eventprocess(MJDetInfo *Dets, MJRunInfo *runInfo, int nChData, BdEvent *ChData[]) {

  int     i, j, k, bad[NCHS], pulser = 0, tmax, wsig_len;
  int     ievt, idet, ch, crate, slot, chan, evlen, siglen, board_type;
  int     dirty_sig = 0, granularity;  // data cleaning result
  long long int  time;
  // unsigned short *head2;
  short          *signal;

  static int *dataId, dataIdRun;
  static int module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  static int totevts=0, builtevts=0, badevts=0, goodevts = 0, pulserevts=0;    // counters
  static int badevts_bytype[32];        // counters
  static FILE *pulserlist = NULL;
  static long long int t1 = 0;

  static PZinfo PZI;
  static CTCinfo CTC;
  static PSAinfo PSA;
  double e_ctc, e_adc, e_raw;
  float  emax, drift;
  float  fsignal[8192], aovere, dcr, lamda;
  int    t0, t95, a_e_good, dcr_good, lamda_good;

  /*
    detector channels:
         0 -  57:  Ge high gain
       100 - 157:  Ge low gain
       158 - 168:  pulser tag
   */

  static int   first = 1;
  static FILE *f_out=0, *f_out2, *f_evl;
  static int   run_t_stamp_0  = -1;   // first timestamp seen in this run
  static int   doPT = 1;              // pulser_tag_init() return value
  static PTag  ptInfo;                // pulser energy and delta-time info for tagging
  static DataClean dcInfo;            // data-cleaning info/limits
  static int   doDC;                  // data_clean_init() return value
  static char  *bad_evt_type[10] =
    {"Unknown module or channel out of range",
     "Unknown detector or channel",
     "Data from disabled channel",
     "Timestamp out of range",
     "Saturated signal",
     "Granularity cut",
     "?","?","?","?"};

  static short  *mat[4096], ee[100], nee, sum, roihit;


  /* --------------- initialization --------------- */
  if (first) {
    first = 0;
    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;

    /* malloc and clear matrix space */
    if ((mat[0] = calloc(4096*4096, sizeof(short))) == NULL) {
      printf("ERROR; cannot malloc mat!\n");
      exit(-1);
    }
    for (i=1; i<4096; i++) mat[i] = mat[i-1] + 4096;

    /* identify data types that we want to decode */
    dataId = runInfo->dataId;
    for (i=0; i<runInfo->idNum; i++) {
      if (strstr(runInfo->decoder[i], "ORRunDecoderForRun")) dataIdRun = dataId[i];
    }

    doPT = pulser_tag_init(Dets, runInfo, &ptInfo);
    doDC = data_clean_init(runInfo, &dcInfo);

    for (i=0; i<32; i++) badevts_bytype[i] = 0;
    if (EVENTLIST) {
      f_evl = fopen("coinc_evl.txt", "w");
      fprintf(f_evl, "#chan e_ctc    timestamp    LDA.bits.bits  detector"
                     " Qx      A/E     DCR  lamda    Run  t95-t0\n");
    }

    /* open ouput file for signals */
    if (WRITE_SIGS) {
      f_out = fopen("s.rms", "w");
    }

    run_t_stamp_0 = ChData[0]->time / 100000000;       // first timestamp, in seconds

    /* read PZ correction info */
    if (!PZ_info_read(runInfo, &PZI)) {
      printf("\n ERROR: No initial pole-zero data read. Does PZ.input exist?\n");
      exit(-1);
    }
    /* read energy correction factors from ctc.input */
    if (!CTC_info_read(runInfo, &CTC)) {
      printf("\n Warning: No initial charge-trapping correction data read. Does ctc.input exist?\n");
      exit(-1);
    }
    /* read A/E, DCR, and lamda values from psa.input */
    if (!PSA_info_read(runInfo, &PSA)) {
      printf("\n ERROR: No PSA data read. Does psa.input exist?\n");
      exit(-1);
    }
  }                   // end of initialization

  /* some versions of the event builder allow for a special call
     to trigger end-of-run cleanup */
  if (runInfo->analysisPass < 0) {
    //return 0;  // do not return; next section (nChData == -99) needs to run.
  }

  /* ------------ finalization / end-of-run handling ------------- */
  if (nChData == -99) {  // special value as flag to trigger finalization

    printf("\n Processed %d built evts, %d ch-evts, %d bad ch-evts (%d%%)  ",
           builtevts, totevts, badevts, (badevts*100)/totevts);
    for (i=1; i<6; i++) printf(" %d", badevts_bytype[i]);
    printf(" %6d pulser evts\n", pulserevts);
    for (i=1; i<6; i++) {
      if (badevts_bytype[i] > 0)
        printf("   >>>> %d bad events of type %d: %s\n",
               badevts_bytype[i], i, bad_evt_type[i]);
    }

    // write out coinc matrix
    //f_out2 = fopen("coinc.mat", "w");
    //f_out2 = fopen("sum3.mat", "w");
    f_out2 = fopen("coinc3.mat", "w");
    fwrite(mat[0], 4096*4096, sizeof(short), f_out2);
    fclose(f_out2);
    return 0;
  }                     // end of post-run processing

  
  /* -=-=-=-=-=-=-=-=-=- process input events  -=-=-=-=-=-=-=-=-=- */

  builtevts++;
  totevts += nChData;

  // fillEvent: fill out entris in BdEvent ChData; chan, det, e, sig...
  fillEvent(runInfo, nChData, ChData, module_lu, det_lu, chan_lu);
  if (doPT) pulser = checkForPulserEvent(Dets, runInfo, nChData, ChData, &ptInfo);
  if (pulser) pulserevts++;
  if (0 && pulser) {  // log pulser events for comparison with GAT analysis
    if (!pulserlist) {
      pulserlist = fopen("pulserlist.out", "w");
      fprintf(pulserlist, "### run   timestamp   number_of_chs\n");
    }
    fprintf(pulserlist, "%4d %12lld %3d\n",
            runInfo->runNumber, ChData[0]->time, nChData);
    if (ChData[0]->time < t1)
      fprintf(pulserlist,
              " >>>>>>> Pulser event times out of order; %lld < %lld; %lld\n",
              ChData[0]->time, t1, t1-ChData[0]->time);
    t1 = ChData[0]->time;
  }

  granularity = checkGranularity(Dets, runInfo, nChData, ChData);
  if (!granularity) return 0;  // process only coincidence events
  nee = 0;
  roihit = 0;

  /* LOOP over channel-events in the built-event again to handle data */
  for (ievt = 0; ievt < nChData; ievt++) {
    evlen  = ChData[ievt]->evlen;
    crate  = ChData[ievt]->crate;
    slot   = ChData[ievt]->slot;
    ch     = ChData[ievt]->ch;
    time   = ChData[ievt]->time;
    idet   = ChData[ievt]->det;
    chan   = ChData[ievt]->chan;
    siglen = ChData[ievt]->siglen;
    board_type = ChData[ievt]->orca_type;

    if (chan > 99) continue;  // ignore LG channels for now
    if (board_type == dataIdRun) {
      continue; // skip all processing of run start/stop/beat records
    }

    /* identify det and chan IDs, and flag some types of bad events */
    bad[ievt] = 0;
    if (ChData[ievt]->mod < 0 || ch > 9) {
      bad[ievt] = 1;                // identify unknown module or channel out of range
    } else {
      if (idet < 0 || chan < 0) {
        bad[ievt] = 2;              // identify unknown detector or channel
      } else if (idet < 99 &&
                 ((chan < 99 && !Dets[idet].HGChEnabled) ||
                  (chan > 99 && !Dets[idet].LGChEnabled))) {
        bad[ievt] = 3;              // identify disabled channel
      }
    }
    if (time < 0 && !bad[ievt]) {
      bad[ievt] = 4;                // identify timestamp out of range
      printf("Bad time: %lld\n", time);
    }
   
    /* ---------------- Ge detector data processing ----------------- */
    if ((board_type == runInfo->dataIdGM ||
         board_type == runInfo->dataIdGA) && !bad[ievt]) {
    /* at this point, we have Gretina4M or 4A data */

      if (chan < 0 || chan > 199) {
        j = module_lu[crate][slot];
        printf("ERROR in eventprocess(): chan = %d!\n"
               "   event number %d; %d of %d; crate, slot, ch: %d %d %d\n"
               "      neighbors: %d %d %d\n",
               chan, totevts, ievt+1, nChData, crate, slot, ch,
               chan_lu[j][ch-1], chan_lu[j][ch], chan_lu[j][ch+1]);
        continue;
      }

      /* extract energies */
      signal = ChData[ievt]->sig;

      /*  ------- do signal-based data cleaning ------- */
      if (doDC && !pulser) {
        dirty_sig = data_clean(signal, chan, &dcInfo);
        /*
          returns 0 for clean signals, or sum of:
                  1 for bad mean baseline value (offset)
                  2 for noisy RMS of baseline
                  4 for bad initial baseline slope
                  8 for bad signal tail slope
                 16 for saturated signal
                 add:
                 32 for granularity cut
                 64 for signal rise too late
        */
        if (granularity) dirty_sig += 32;
      } else {
        dirty_sig = 0;
      }

      // ----------------------------------------------------------------
      /* find t100 and t95*/
      int t100 = 700;                 // FIXME? arbitrary 700?
      for (i = t100+1; i < 1500; i++)
        if (signal[t100] < signal[i]) t100 = i;
      /* get mean baseline value */
      int bl = 0;
      for (i=300; i<400; i++) bl += signal[i];
      bl /= 100;
      for (t95 = t100-1; t95 > 500; t95--)
        if ((signal[t95] - bl) <= (signal[t100] - bl)*19/20) break;

      if (t100   > 1300 ||   // important for cleaning, gets rid of pileup
          t95+50 > 1300) {   // dcr trap extends past end of signal
        if (VERBOSE || DEBUG)
          printf(">>> Error: chan %d signal at timestamp %lld is too late!\n", chan, time);
        dirty_sig += 64;
        continue;
      }

      /* do (optional) INL  correction */
      if (DO_INL) {
        if (inl_correct(Dets, runInfo, signal+10, fsignal+10, siglen-10, chan)) {
          printf(" >>> inl_correct return error for chan %d\n", chan);
          return -1;
        }
      } else {
        for (i=0; i<siglen; i++) fsignal[i] = signal[i];
      }

      /* do fitting of pole-zero parameters to get lamda (~ DCR) */

      float chisq, lamda1, frac2;
      int tlo = t100+50, thi = 1990;            // FIXME; variable length
      // if (thi > tlo + 700) thi = tlo + 700;  // FIXME; check performance
      chisq = pz_fitter(fsignal, tlo, thi, chan, &PZI, &lamda1, &frac2, &lamda);
      if (chisq < 0.01 || chisq > 50.0) continue;      // FIXME     // fit failed, or bad fit chisq
      lamda = (0.01 / PZI.tau[chan] - lamda) * 1.5e6 + 0.3;  // 0.3 fudge to get mean ~ 0
      lamda *= 8.0;                                          // scaled to roughly match DCR at 2 MeV

      /* do (required) PZ correction */
      if (PZ_fcorrect(fsignal+10, siglen-10, chan, &PZI))
        printf(" >>> PZ_fcorrect return error for chan %d\n", chan);

      /* ------- do DT correction ------- */
      emax = float_trap_max(fsignal, &tmax, 401, 250)/401.0;
      e_ctc = get_CTC_energy(fsignal, siglen, chan, Dets, &PSA,
                             &t0, &e_adc, &e_raw, &drift, CTC.e_dt_slope[chan]);
      if (t0 < 0) {
        if (VERBOSE) printf("t0 determination failed; chan %3d, timestamp %lld\n", chan, time);
        if (VERBOSE) printf("chan, t0, tmax: %d %d %d; e: %.2f\n",
                            chan, t0, tmax+TRAP_RISE+TRAP_FLAT, e_raw);
        continue;
      } else if (e_ctc < 0.001) {
        printf("E_ctc = %.1f < 1 eV!\n", e_ctc);
        if (VERBOSE) printf("chan, t0, tmax: %d %d %d; e, drift: %.2f %.2f\n",
                            chan, t0, tmax+TRAP_RISE+TRAP_FLAT, e_raw, drift);
        continue;
      }

      if (EVENTLIST && e_ctc > EVENTLIST && e_ctc < EVENTLISTMAX && chan < 100)
          roihit = 1;

      if (0 && dirty_sig) {  // DIRTY signals
        if (EVENTLIST && e_ctc > EVENTLIST && e_ctc < EVENTLISTMAX && chan < 100) {
          fprintf(f_evl, "%3d %7.1f %15lld 000.%d%d%d%d.%d%d%d%d", chan, e_ctc, time,
                 DBIT(7), DBIT(6), DBIT(5), DBIT(4), DBIT(3), DBIT(2), DBIT(1), DBIT(0));
          if (chan > 99)
            fprintf(f_evl, " %s LG Q0 %24d\n", Dets[chan-100].StrName, runInfo->runNumber);
          else
            fprintf(f_evl, " %s HG Q0 %24d\n", Dets[chan].StrName, runInfo->runNumber);
        }

      } else if (!pulser) {          // CLEAN signals

        /* find A/E */
        if (e_ctc > 50 && t0 > 200 && t0 < 1500) {
          float s2 = float_trap_max_range(fsignal, &tmax, PSA.a_e_rise[chan], 0, t0-20, t0+300); // FIXME: hardwired range??
          aovere = s2 * PSA.a_e_factor[chan] / e_raw;
          /* do quadratic fit/interpolation over +- one sample, to improve max A/E determination
             does this help? seems to work well in some cases, in other cases not so much?
          */
          if (AoE_quad_int) {
            float s1, s3, a, b, aoe;
            int rise = PSA.a_e_rise[chan];
            s1 = s2 - (fsignal[tmax  ] - 2.0*fsignal[tmax   + rise] + fsignal[tmax   + 2*rise]);
            s3 = s2 + (fsignal[tmax+1] - 2.0*fsignal[tmax+1 + rise] + fsignal[tmax+1 + 2*rise]);
            b = s2 - (s1+s3)/2.0;
            a = s2 - s1 + b;
            aoe = s1 + a*a/(4.0*b);
            if (aoe > s2) aovere = aoe * PSA.a_e_factor[chan] / e_raw;
          }
          if (runInfo->flashcam) aovere /= 2.0;
        } else {
          aovere = 0;
          if (e_raw > 1000)
            printf("Error getting A/E: chan %d   e_raw, e_ctc = %.0f, %.0f  t0 = %d\n", chan, e_raw, e_ctc, t0);
        }

        /* ---- This next section calculates the GERDA-style A/E ---- */
        if (PSA.gerda_aoe[chan]) {
          float  ssig[6][2000];
          for (k=100; k<2000; k++) ssig[0][k] = fsignal[k] - fsignal[300];
          for (j=1; j<6; j++) {  // number of cycles
            for (i=200; i < 1500; i++) {
              ssig[j][i] = 0.0;
              for (k=0; k<10; k++) {
                ssig[j][i] += ssig[j-1][i+k-5];
              }
              ssig[j][i] /= 10.0;
            }
          }
          aovere = 0;
          for (k=250; k<1450; k++) {
            if (aovere < ssig[5][k] - ssig[5][k-1]) aovere = ssig[5][k] - ssig[5][k-1];
          }
          aovere *= 12.71*1593/e_raw;
        }
        /* ---- end of GERDA-style A/E ---- */

        /* get DCR from slope of PZ-corrected signal tail */
        dcr = float_trap_fixed(fsignal, t95+50, 100, 500) / 25.0;

        /* do drift-time and energy corrections to A/E, DCR, lamda */
        float dtc = drift*3.0*2614.0/e_adc;   // in (x)us, for DT- correction to A/E
        if (chan%100 == 50) dtc *= 3.0;       // special hack for detector 50; has especially strong variation
        float ec  = e_ctc/1000.0 - 2.0;       // in MeV, for E-correction to A/E
        aovere += PSA.ae_dt_slope[chan]    * dtc;
        aovere += PSA.ae_e_slope[chan]     * ec;
        dcr    -= PSA.dcr_dt_slope[chan]   * drift;
        lamda  -= PSA.lamda_dt_slope[chan] * drift;

        /* ------- do A/E cut ------- */
        /* find A/E */
        a_e_good = 0;
        if (1) {
          // adjust for energy dependence of cut
          // assume FWHM of A distribution from series noise is ~ 2 keV
          aovere -= 2.0 * (2.75 * (1.0 - 1593.0/e_ctc));  // 1593 keV = DE
        }
        if (aovere >= PSA.ae_cut[chan]) a_e_good = 1;

        /* ------- do DCR cut ------- */
        dcr_good = 0;
        if (dcr > -100 && e_ctc > 1000 && dcr < PSA.dcr_lim[chan]) dcr_good = 1;

        /* ------- do lamda cut ------- */
        lamda_good = 0;
        if (lamda > -100 && e_ctc > 1000 && lamda < PSA.lamda_lim[chan]) lamda_good = 1;

        // ----------------------------------------------------------------

        if (EVENTLIST && e_ctc > EVENTLIST && e_ctc < EVENTLISTMAX && chan < 100) {
          fprintf(f_evl, "%3d %7.1f %15lld %d%d%d.%d%d%d%d.%d%d%d%d",
                  chan, e_ctc, time, lamda_good, dcr_good, a_e_good,
                  DBIT(7), DBIT(6), DBIT(5), DBIT(4), DBIT(3), DBIT(2), DBIT(1), DBIT(0));
          if (chan > 99)
            fprintf(f_evl, " %s LG Q%d  %8.1f %6.1f %6.1f   %d %d\n",
                    Dets[chan-100].StrName, a_e_good + 2*dcr_good + 4*lamda_good,
                    aovere, dcr, lamda, runInfo->runNumber, t95-t0);
          else
            fprintf(f_evl, " %s HG Q%d  %8.1f %6.1f %6.1f   %d %d\n",
                    Dets[chan].StrName, a_e_good + 2*dcr_good + 4*lamda_good,
                    aovere, dcr, lamda, runInfo->runNumber, t95-t0);
        }
      }
      // ----------------------------------------------------------------

      if (e_ctc > 5 && e_ctc < 4090)
        ee[nee++] = e_ctc + 0.5;

      /* write out good (or selected) signals */

      /*
        if (WRITE_SIGS) {
        wsig_len = siglen;
        if (0) {
        // optionally subtract initial baseline
        j = (signal[20]+signal[21]+signal[22]+signal[23]+signal[24])/5;
        for (i=0; i<wsig_len; i++) signal[i] -= j;
        }
        if (1) {
        // optionally add some extra information at the start
        signal  -= 8; // 8 extra words (shorts)
        wsig_len += 8;
        for (i=0; i<8; i++) signal[i] = 0;
        signal[2] = (short) (e_ctc + 0.5);
        head2  = ((unsigned short *) ChData[ievt]->evbuf) + 4;
        signal[4] = (short) ((head2[8]>>11)&0xf)*100;
        sigqnal[5] = (short) dirty_sig;
        signal[6] = (short) chan;
        }

        if (e_ctc > 5)
        write_sig(signal, wsig_len, chan, f_out);
        }
      */
    }
    /* ------------- end of Ge detector data processing -------------- */

    /* deal with bad events identified so far */
    if (bad[ievt]) {  // bad event, not just dirty/cleaned
      badevts++;
      badevts_bytype[bad[ievt]]++;
      return 2;
    }

  } // END of loop over channel-events

  // increment coincidence matrix
  // if (nee > 1 && nee < 4) {
  if (nee > 2 && nee < 5) {
    mat[0][nee]++;
    mat[0][0] += nee*(nee-1)/2;
    mat[0][100+nee] += nee*(nee-1)/2;
    for (i=0; i<nee; i++) mat[1][ee[i]]++; 
    sum = ee[0];
    for (i=1; i<nee; i++) sum += ee[i];
    if (sum < 4095) mat[2][sum]++;
    if (nee > 5) printf("nee = %d\n", nee);

    // increment coincidence matrix
    if (0) {
      if (ee[0] < 0 || ee[0] > 4095) printf("aaa %d %d!!!\n", 0, ee[0]);
      for (i=1; i<nee; i++) {
        if (ee[i] < 0 || ee[i] > 4095) printf("aaa %d %d!!!\n", i, ee[i]);
        for (j=0; j<i; j++) {
          mat[ee[i]][ee[j]]++;
          mat[ee[j]][ee[i]]++;
        }
      }
    }
    // increment sum matrix
    if (1 && nee < 5 && sum < 4095) {
      for (i=0; i<nee; i++) mat[sum][ee[i]]++; 
    }
  }

  // if (nee > 6) {
  if (nee < 4 && roihit) {
    for (ievt = 0; ievt < nChData; ievt++) {
      crate  = ChData[ievt]->crate;
      time   = ChData[ievt]->time;
      chan   = ChData[ievt]->chan;
      siglen = ChData[ievt]->siglen;
      signal = ChData[ievt]->sig;
      if (chan > 99) continue;
      if (DO_INL) {
        if (inl_correct(Dets, runInfo, signal, fsignal, siglen, chan))
          printf(" >>> inl_correct return error for chan %d\n", chan);
        if (PZ_fcorrect(fsignal, siglen, chan, &PZI))
          printf(" >>> PZ_correct return error for chan %d\n", chan);
      } else {
        if (PZ_correct(signal, fsignal, siglen, chan, &PZI))
          printf(" >>> PZ_correct return error for chan %d\n", chan);
      }
      emax = float_trap_max(fsignal, &tmax, 401, 250)/401.0;
      e_ctc = get_CTC_energy(fsignal, siglen, chan, Dets, &PSA,
                             &t0, &e_adc, &e_raw, &drift, CTC.e_dt_slope[chan]);
      dirty_sig = data_clean(signal, chan, &dcInfo);
      fprintf(f_evl, "%3d %7.1f %15lld %d%d%d.%d%d%d%d.%d%d%d%d",
              chan, e_ctc, time, lamda_good, dcr_good, a_e_good,
              DBIT(7), DBIT(6), DBIT(5), DBIT(4), DBIT(3), DBIT(2), DBIT(1), DBIT(0));
      fprintf(f_evl, " %s HG   %2d  %2d\n",
              Dets[chan].StrName, ievt, nee);

      wsig_len = siglen;
      if (WRITE_SIGS) {
        if (1) {
          // optionally subtract initial baseline
          j = (signal[20]+signal[21]+signal[22]+signal[23]+signal[24])/5;
          for (i=0; i<wsig_len; i++) signal[i] -= j;
        }
        if (1) {
          // optionally add some extra information at the start
          signal  -= 8; // 8 extra words (shorts)
          wsig_len += 8;
          for (i=0; i<8; i++) signal[i] = 0;
          signal[2] = (short) (e_ctc + 0.5);
          // head2  = ((unsigned short *) ChData[ievt]->evbuf) + 4;
          // signal[4] = (short) ((head2[8]>>11)&0xf)*100;
          signal[5] = (short) dirty_sig;
          signal[6] = (short) chan;
        }
        write_sig(signal, wsig_len, chan, f_out);
      }
    }
  }

  goodevts++;
  if (totevts % 20000 < (totevts-nChData) % 20000) {
    printf(" Processed %d built evts, %d ch-evts, %d bad ch-evts (%d%%)  ",
           builtevts, totevts, badevts, (badevts*100)/totevts);
    for (i=1; i<6; i++) printf(" %d", badevts_bytype[i]);
    printf(" %6d pulser evts\n", pulserevts);
  }

  if (pulser) return 1;
  return  0;
} /* eventprocess() */

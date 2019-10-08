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
/*
 *  NOTE: define compile-time option DO_PSA for version
 *         where we want to apply A/E and DCR cuts
 *         It is a little slower, and requires extra input to define the cuts
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
#define SUM_UP_PULSER_CNTS  0  // set to 1 to sum pulser counts over many data subsets,
                               //     and delete ptagtot.dat to reset sums to zero

#define MAKE_2D      0        // make a file with data for 2d plots
#define WRITE_SIGS   0        // write selected signals to an otput file
#define HIS_COUNT    2200     // number of spectra in BOTH his[][] and his2[][] arrays
#define DBIT(j)      (1 - (1&(dirty_sig>>j)))

#define SUBTR_DCR_MEAN (e_adc > 20 && e_adc < 8000 && \
                        ((chan < 100 && SUBTRACT_MEAN_DCR_HG) || \
                         (chan > 99  && SUBTRACT_MEAN_DCR_LG)))

/*  ------------------------------------------------------------ */

int ep_finalize(MJDetInfo *Dets, MJRunInfo *runInfo, int **his,
                int *on_bd_rise, int totevts, PTag  *pt, DataClean *dcInfo);

/*  ------------------------------------------------------------ */

int eventprocess(MJDetInfo *Dets, MJRunInfo *runInfo, int nChData, BdEvent *ChData[]) {

  int     i, j, k, n, bad[NCHS], pulser = 0, tmax, t_1;
  int     ievt, idet, ch, crate, slot, chan, evlen, siglen, board_type;
  int     e_onbd, e_offline, e_trapmax, de;
  int     dirty_sig = 0, granularity;  // data cleaning result
  double  s1, s2;
  long long int  time;
  unsigned short *head2;
  short          *signal;

  static int *dataId, dataIdRun;
  static int module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  static int totevts=0, builtevts=0, badevts=0, goodevts = 0, pulserevts=0;    // counters
  static int badevts_bytype[32];        // counters
  //static int vetoevts = 0, coincevts = 0;
  static long long int t1 = 0;
  static int *his[HIS_COUNT];           // storage for histograms

  static PTag      ptInfo, pt2;         // pulser energy and delta-time info for tagging
  static DataClean dcInfo;              // data-cleaning info/limits
#ifdef DO_PSA
  static PZinfo PZI;                    // pole-zero info
  static CTCinfo CTC;                   // charge-trapping correction info
  static PSAinfo PSA;                   // pulse-shape analysis info and limits
  static int  *his2[HIS_COUNT];         // storage for final-result histograms
  static int  *dcr_mean[200], mean_dcr_ready = 0;
  static FILE *f_evl, *fp;              // file pointer for event list
#endif
  static FILE *pulserlist = NULL;
  static FILE *f_sig, *f_2d;

  static int   first = 1;
  static FILE *f_out=0, *bf_out=0, *f_in;
  static int   sigout_channum = -1;     // sig channel flag to be defined based on run arguments
  static int   sigout_ptag    = -1;     // pulser-yes/no flag to be defined based on run arguments
  static int   sigout_ctag    = -1;     // clean-yes/no  flag to be defined based on run arguments
  static int   run_t_stamp_0  = -1;     // first timestamp seen in this run
  static int   doPT = 1;                // pulser_tag_init() return value
  static int   doDC;                    // data_clean_init() return value
  static int   on_bd_rise[200], on_bd_flat[200];  // on-board trap rise and flat times (10ns)
  static int   on_bd_change[200];       // flag for consecutive wrong rise-times for on-board trap
  static char  *bad_evt_type[10] =
    {"Unknown module or channel out of range",
     "Unknown detector or channel",
     "Data from disabled channel",
     "Timestamp out of range",
     "Saturated signal",
     "Granularity cut",
     "Bad t0 or E_ctc",
     "Bad chisq for fit to tail",
     "?","?"};

  /*
    detector channels:
         0 -  57:  Ge high gain
       100 - 157:  Ge low gain
       158 - 168:  pulser tag
    his channels:
         0 - 199   pulser-removed trapmax (4-2-4)
       200 - 399   pulser-tagged trapmax (4-2-4)
       400 - 599   dirty energy (e_ctc or trapmax) after data cleaning; ADC units
       600 - 799   clean energy (e_ctc or trapmax) after data cleaning; ADC units
       800 - 999   dirty energy (e_ctc or trapmax) after data cleaning; 0.5 keV/ch
      1000 - 1199  clean energy (e_ctc or trapmax) after data cleaning; 0.5 keV/ch
      1200 - 1399  A/E histograms
      1400 - 1599  DCR and lamda histograms
      1600 - 1799  baseline mean, RMS, and slope (+4000 when trapmax is < 50)
      1800 - 1999  tmax from trapmax (4-2-4) and 2000+t_1 (1% time)
                     and 6000+baseline_trap (E=0, 4-1.5-4) for FWHM, 0.05 ADC units
      2000 - 2199  trap difference (on-board - trapmax) for trapmax > 100
                     in bins of 0.01 (+1000) and 0.2 (+3000) ADC
                     and (t90 - t0) time difference for different A/E and DCR cuts 
   */


  /* --------------- initialization --------------- */
  if (first) {
    first = 0;
    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
    if (MAKE_2D)    f_2d  = fopen("2d.dat", "w");
    if (WRITE_SIGS) f_sig = fopen("sig.mat", "w");

    /* identify data types that we want to decode */
    dataId = runInfo->dataId;
    for (i=0; i<runInfo->idNum; i++) {
      if (strstr(runInfo->decoder[i], "ORRunDecoderForRun")) dataIdRun = dataId[i];
    }

    doPT = pulser_tag_init(Dets, runInfo, &ptInfo);
#ifndef QUIET
    doDC = data_clean_init(runInfo, &dcInfo);
#endif

    for (i=0; i<32; i++) badevts_bytype[i] = 0;
    /* malloc and clear histogram space */
    if ((his[0] = calloc(HIS_COUNT * 8192, sizeof(int))) == NULL) {
      printf("ERROR in eventprocess(); cannot malloc his!\n");
      exit(-1);
    }
    for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1] + 8192;
#ifdef DO_PSA
    if ((his2[0] = calloc(HIS_COUNT * 8192, sizeof(int))) == NULL) {
      printf("ERROR in eventprocess(); cannot malloc his!\n");
      exit(-1);
    }
    if ((dcr_mean[0] = calloc(200*8192, sizeof(int))) == NULL) {
      printf("ERROR in PSAcal.c; cannot malloc dcr_mean!\n");
      exit(-1);
    }
    for (i=1; i<200; i++) dcr_mean[i] = dcr_mean[i-1] + 8192;
    if (SUBTRACT_MEAN_DCR_HG || SUBTRACT_MEAN_DCR_LG) {
      if ((fp = fopen("INL_DCR_input.sec", "r"))) {
        fread(dcr_mean[0], 8192*sizeof(int), 200, fp);        // get or skip smoothed values
        if (SUBTRACT_MEAN_DCR_HG > 1) {
          fread(dcr_mean[0], 8192*sizeof(int), 100, fp);      // get unsmoothed HG values
          if (SUBTRACT_MEAN_DCR_LG > 1)
            fread(dcr_mean[100], 8192*sizeof(int), 100, fp);  // get unsmoothed LG values
        } else if (SUBTRACT_MEAN_DCR_LG > 1) {
          fread(dcr_mean[100], 8192*sizeof(int), 100, fp);    // skip unsmoothed HG values
          fread(dcr_mean[100], 8192*sizeof(int), 100, fp);    // get unsmoothed LG values
        }
        fclose(fp);
        mean_dcr_ready = 1;   // all is ready
        printf("\nRead mean INL_DCR values from INL_DCR_input.sec\n\n");
      } else {
        mean_dcr_ready = -1;  // -1 indicates need to take further action later
        printf("\nReading of mean INL_DCR values from INL_DCR_input.sec was UNSUCCESSFUL\n\n");
      }
    }

    for (i=1; i<HIS_COUNT; i++) his2[i] = his2[i-1] + 8192;
    if (EVENTLIST) {
      f_evl = fopen("evl.txt", "w");
      fprintf(f_evl, "#chan e_ctc    timestamp    LDA.bits.bits  detector"
                     " Qx      A/E     DCR  lamda    Run  t90-t0 t100-t90  A/E,DCR,lamda - limit\n");
    }
#endif

    /* decode/process command-line agruments */
    if (runInfo->argc > 3 &&
        strstr(runInfo->argv[1], "-s")) {
      sigout_channum = atoi(runInfo->argv[2]);  // select channel # for signals output
      if (strstr(runInfo->argv[1], "-s0")) sigout_channum = -2;  // disable signal output
      if (strstr(runInfo->argv[1], "p0")) sigout_ptag = 0;  // select non-pulser
      if (strstr(runInfo->argv[1], "p1")) sigout_ptag = 1;  // select puser
      if (strstr(runInfo->argv[1], "c0")) sigout_ctag = 0;  // select non-clean
      if (strstr(runInfo->argv[1], "c1")) sigout_ctag = 1;  // select clean
    }
    if (sigout_ptag>0 && sigout_ctag>0)
      printf("Selecting clean pulser sigs from ch %d\n", sigout_channum);
    if (sigout_ptag>0 && !sigout_ctag)
      printf("Selecting dirty pulser sigs from ch %d\n", sigout_channum);
    if (!sigout_ptag && sigout_ctag>0)
      printf("Selecting clean non-pulser sigs from ch %d\n", sigout_channum);
    if (!sigout_ptag && !sigout_ctag)
      printf("Selecting dirty non-pulser sigs from ch %d\n", sigout_channum);
    if (sigout_ptag<0 && sigout_ctag>0)
      printf("Selecting clean sigs from ch %d\n", sigout_channum);
    if (sigout_ptag<0 && !sigout_ctag)
      printf("Selecting dirty sigs from ch %d\n", sigout_channum);
    if (sigout_ptag>0 && sigout_ctag<0)
      printf("Selecting pulser sigs from ch %d\n", sigout_channum);
    if (!sigout_ptag && sigout_ctag<0)
      printf("Selecting non-pulser sigs from ch %d\n", sigout_channum);

    /* open ouput file for signals */
#ifndef QUIET
    if (sigout_channum >= -1) {
      f_out = fopen("s.rms", "w");
      if (doDC) bf_out = fopen("bs.rms", "w"); // for bad signals
      if (DEBUG) inl_correct(Dets, runInfo, 0, 0, 0, 0);
    }
#endif

    run_t_stamp_0 = ChData[0]->time / 100000000;       // first timestamp, in seconds

    /* read on-board trapezoid rise and flat-top times from file */
    for (i=0; i<200; i++) {
      on_bd_rise[i] = 401;
      on_bd_flat[i] = 201;
      on_bd_change[i] = 0;
    }
    if ((f_in = fopen(OB_TRAP_FILENAME, "r"))) {
      int   r, f, m=0;
      char  line[256];
      while (fgets(line, sizeof(line), f_in)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d %d %d", &j, &k, &r, &f) != 4) break;
        if (j > k || j < 0 || k > 99) break;
        for (i=j; i<=k; i++) {
          if (i >=0 && i < runInfo->nGe) {
            on_bd_rise[i] = r;
            on_bd_flat[i] = f;
            m++;
          }
        }
      }
      fclose(f_in);
      printf(" %d onboard trapezoids defined from file %s\n", m, OB_TRAP_FILENAME);
    } else {
      printf("No onboard trap times defined, using 401 / 201\n");
    }

#ifdef DO_PSA
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
#endif

  }  /* ---------------------------- end of initialization ----------------------------
      * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

  /* some versions of the event builder allow for a special call
     to trigger end-of-run cleanup */
  if (runInfo->analysisPass < 0) {
    if (SUM_UP_PULSER_CNTS) {     // set to 1 for summation of pulser counts over many data subsets
      FILE *f_out2;
      printf("Reading ptagtot.dat\n"); fflush(stdout);
      if ((f_out2 = fopen("ptagtot.dat", "r"))) {  // delete ptagtot.dat to reset sums to zero
        fread(&pt2, sizeof(pt2), 1, f_out2);
        fclose(f_out2);
        for (chan=0; chan < runInfo->nGe; chan++) {
          for (i=0; i<8; i++) pt2.nevts[chan][i] += ptInfo.nevts[chan][i];
        }
        printf("Writing ptagtot.dat\n"); fflush(stdout);
        f_out2 = fopen("ptagtot.dat", "w");
        fwrite(&pt2, sizeof(pt2), 1, f_out2);
        fclose(f_out2);
        printf("\n"
               "                           Pulser Counts                         LiveTime\n"
               "        Detector        HG    LG  Both  Expected     diff      Either   (Both)\n");
        for (chan=0; chan < runInfo->nGe; chan++) {
  /* ptag.nevts[chan][0]: energy-ungated HG pulser cts
     ptag.nevts[chan][1]: energy-ungated LG pulser cts
     ptag.nevts[chan][2]: energy-ungated HG&&LG pulser cts
     ptag.nevts[chan][3]: energy-gated HG pulser cts
     ptag.nevts[chan][4]: energy-gated LG pulser cts
     ptag.nevts[chan][5]: energy-gated HG&&LG pulser cts
     ptag.nevts[chan][6]: expected energy-ungated pulser cts (from finding max over the CC)
   */
          if (pt2.nevts[chan][6] < 2 ||   // no pulser expected
              (pt2.nevts[chan][0] == 0 && pt2.nevts[chan][1] == 0)) continue;  // no pulser
          n = pt2.nevts[chan][0] + pt2.nevts[chan][1] - pt2.nevts[chan][2]; // either
          printf("%3d %8s %s %6d %5d %5d %7d %6d %4d %10.4f  (%6.4f)\n",
                 chan, Dets[chan].DetName, Dets[chan%100].StrName,
                 pt2.nevts[chan][0], pt2.nevts[chan][1],
                 pt2.nevts[chan][2], pt2.nevts[chan][6],
                 pt2.nevts[chan][6] - pt2.nevts[chan][0],
                 pt2.nevts[chan][6] - pt2.nevts[chan][1],
                 (float) n / (float) pt2.nevts[chan][6],
                 (float) pt2.nevts[chan][2] / (float) pt2.nevts[chan][6]);
        }
      } else {
        printf("Writing ptagtot.dat\n"); fflush(stdout);
        f_out2 = fopen("ptagtot.dat", "w");
        fwrite(&ptInfo, sizeof(pt2), 1, f_out2);
        fclose(f_out2);
      }
    }
    //return 0;  // do not return; next section (nChData == -99) needs to run.
  }

  /* --------------------  finalization / end-of-run handling  ---------------------
   * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
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

    if (f_out) {fclose(f_out); f_out = 0;}
    if (bf_out) {fclose(bf_out); bf_out = 0;}
    
    int ret_value = ep_finalize(Dets, runInfo, his, on_bd_rise, totevts, &ptInfo, &dcInfo);

#ifdef DO_PSA
    // write out final-result histograms (his2) to fin.rms
    FILE *f_out2 = fopen("fin.rms", "w");
    for (i=0; i<HIS_COUNT; i++) {
      char spname[256];
      if (i < 200) {
        sprintf(spname, "%d; E, no cuts [0.5 keV]; ch %d run %d", i, i%200, runInfo->runNumber);
      } else if (i < 400) {
        sprintf(spname, "%d; E, A/E cut [0.5 keV]; ch %d", i, i%200);
      } else if (i < 600) {
        sprintf(spname, "%d; E, DCR cut [0.5 keV]; ch %d", i, i%200);
      } else if (i < 800) {
        sprintf(spname, "%d; E, lamda cut [0.5 keV]; ch %d", i, i%200);
      } else if (i < 1000) {
        sprintf(spname, "%d; E, all cuts [0.5 keV]; ch %d", i, i%200);
      } else if (i < 1200) {
        sprintf(spname, "%d; E, HG/LG only, all cuts [0.5 keV]; ch %d", i, i%200);
      } else if (i < 1400) {
        sprintf(spname, "%d; E, Dirty signals incl. gran.; ch %d", i, i%200);
      } else if (i < 1600) {
        sprintf(spname, "%d; E, Dirty signals excl. gran.; ch %d", i, i%200);
      } else {
        sprintf(spname, "%d;", i);
      }
      write_his(his2[i], 8192, i, spname, f_out2);
    }
    fclose(f_out2);

    if (mean_dcr_ready  < 0)
      printf("\nWarning reminder: No INL_DCR_input.sec was read, but was expected!\n"
             " Maybe mv INL_DCR_output.sec INL_DCR_input.sec ?\n\n");
#endif

    return ret_value;
  }    // end of post-run processing

  
  /* ---------------------------  process input events   ---------------------------
   * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

  builtevts++;
  totevts += nChData;

  // fillEvent: fill out entries in BdEvent ChData; chan, det, e, sig...
  fillEvent(runInfo, nChData, ChData, module_lu, det_lu, chan_lu);
  if (doPT) pulser = checkForPulserEvent(Dets, runInfo, nChData, ChData, &ptInfo);
  if (pulser) pulserevts++;
  if (SUM_UP_PULSER_CNTS && pulser) {  // log pulser events for comparison with GAT analysis
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

  j = (ChData[0]->time / 100000000) - run_t_stamp_0;  // time since start of run, in seconds
  his[95+pulser][j%8000]++;

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
    if (DEBUG) {
      printf("ievt, len, cr, slot, time: %d %d %d %d %lld\n",
             ievt, evlen, crate, slot, time);
      fflush(stdout);
    }

    j = (time / 100000000) - run_t_stamp_0;
    his[97+pulser][j%8000]++;
    his[99][300+crate]++;
 
    /* identify det and chan IDs, and flag some types of bad events */
    bad[ievt] = 0;
    if (board_type == dataIdRun) {
      continue; // skip all processing of run start/stop/beat records
    }

    if (ChData[ievt]->mod < 0 || ch > 9) {
      bad[ievt] = 1;                // identify unknown module or channel out of range
      idet = chan = 300 + bad[ievt];
    } else {
      if (idet < 0 || chan < 0) {
        bad[ievt] = 2;              // identify unknown detector or channel
        idet = chan = 300 + bad[ievt];
      } else if (idet < 99 &&
                 ((chan < 99 && !Dets[idet].HGChEnabled) ||
                  (chan > 99 && !Dets[idet].LGChEnabled))) {
        bad[ievt] = 3;              // identify disabled channel
        idet = chan = 300 + bad[ievt];
      } else {
        his[99][chan]++;
      }
    }
    if (time < 0 && !bad[ievt]) {
      bad[ievt] = 4;                // identify timestamp out of range
      printf("Bad time: %lld\n", time);
    }

    if (VERBOSE)
      printf("len, ch, cr, slot, chan: %d %d %d %d %d\n", evlen, ch, crate, slot, chan);
    
    /* ------------------------ Ge detector data processing --------------------------
     * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

    if ((board_type == runInfo->dataIdGM ||
         board_type == runInfo->dataIdGA) && !bad[ievt]) {
    /* at this point, we have Gretina4M or 4A data */

      if (chan < 0 || chan > 199) {
        printf("ERROR in eventprocess(): chan = %d!\n"
               "   event number %d; %d of %d; crate, slot, ch: %d %d %d\n"
               "      neighbors: %d %d %d\n",
               chan, totevts, ievt+1, nChData, crate, slot, ch,
               chan_lu[j][ch-1], chan_lu[j][ch], chan_lu[j][ch+1]);
        continue;
      }

      /* extract energies */
      signal = ChData[ievt]->sig;
      head2  = ((unsigned short *) ChData[ievt]->evbuf) + 4;
      e_onbd = (((head2[8] & 0x1ff) << 16) | head2[7]);
      /* use stored trapmax (401-200-401) as real energy from here on out */
      e_trapmax = ChData[ievt]->e;

      /* histogram the energy */
      if (e_trapmax >= 0 && e_trapmax < 8192) {
        his[200*pulser+chan][e_trapmax]++;
        /* process pulser signa shapes */
        if (pulser) {
          /*
          // histogram an approximation to the pulser DCR value
          int t90, t100, bl;
          float dcr;
          // find pulser t90
          bl = 0;
          for (i=300; i<400; i++) bl += signal[i];
          bl /= 100;
          t100 = 900;
          for (i=901; i<1400; i++)
            if (signal[t100] < signal[i]) t100 = i;
          for (t90 = t100-1; t90 > 500; t90--)
            if ((signal[t90]-bl) <= (signal[t100] - bl)*19/20) break;
          if (t90+50 <= 1300) {
            dcr = 6000.0 + trap_fixed(signal, t90+250, 100, 300) / 20.0;  // FIXME; see below
            dcr += e_trapmax * 5.0 * 4.0/72.0;  // FIXME; delete?
            // physics signals use:
            // dcr = float_trap_fixed(fsignal, t90+50, 100, 500) / 25.0;
            if (chan > 99) dcr = 3.0*dcr - 12000.0;
            his[200+chan][(int) dcr]++;
          }
          */
        }
      }
      // trig_sign = (head2[8] >> 12) & 1;
      if (DEBUG) printf(" PT%d > chan %3d e_onbd = %4d  e_tmax = %4d at %4d\n",
                        pulser, chan, e_onbd/on_bd_rise[idet], e_trapmax, tmax);

      /* Skip much of the following code for pulser tag chs (idet > runInfo->nGe) */
      if (idet < runInfo->nGe) {
        /* calculate trap_max with same times as on-board trap */
        e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);

        /* ----- check for bad trap rise time; 401 -> 465 or vice versa ----- */
        if (1 && e_trapmax > 200 && e_trapmax < 6000) {
          s1 = (float) e_onbd / (float) e_offline;
          if (on_bd_rise[idet] == 401 && s1 > 1.13 && s1 < 1.17) {
            // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
            if (on_bd_change[idet]++ > 2) {
              on_bd_rise[idet] = 465;
              e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);
              printf(" ****> Det %d on-board rise time: %d -> %d\n", idet, 401, 465);
              on_bd_change[idet] = 0;
            }
          } else if (on_bd_rise[idet] == 465 && 1.0/s1 > 1.13 && 1.0/s1 < 1.17) {
            // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
            if (on_bd_change[idet]++ > 2) {
              on_bd_rise[idet] = 401;
              e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);
              printf(" ****> Det %d on-board rise time: %d -> %d\n", idet, 465, 401);
              on_bd_change[idet] = 0;
            }
          } else {
            on_bd_change[idet] = 0;
          }
        }
        /* ----- end of special hack for on-board rise times ----- */

        // note special hack for ch 16; no big pulser signals
        if ((e_trapmax > 60 && e_onbd > 50) ||
            (pulser && chan%100 == 16 && e_trapmax > ptInfo.elo[chan])) {
          // note special hack for ch 16; no big pulser signals

          /* histogram difference between the on-board and offline energies */
          s1 = 100.0 * (e_onbd - e_offline) / (double) on_bd_rise[idet];
          de = (int) (s1 + 1000.5);     // onboard - offline difference, 0.01 ADC units
          if (de > 0 && de < 2000) his[2000+chan][de]++;  // 0.01 ADC units
          de = de/20 + 2950;            // units are now 0.2 ADC cts
          if (de < 2001) de = 2001;
          if (de > 3999) de = 3999;
          his[2000+chan][de]++;         // 0.2 ADC units

          /* histogram tmax and t_1 */
          if (tmax > 0 && tmax < 2000) his[1800+chan][tmax]++;       // trap t_max
          t_1 = sig_frac_time(signal, 0.01, 30, tmax+400);
          if (t_1 > 0 && t_1 < 2000)   his[1800+chan][2000+t_1]++;   // 1% signal time

          /* histogram baseline E=0 trap (4-1.5-4) */
          i = (trap_fixed(signal, 20, 400, 100) + 6000*20)/20;  // 0.05-ADC units
          if (i > 4000 && i < 8000) his[1800+chan][i]++;
        }

        /* histogram resting average baseline and RMS */
#ifdef SHORT_BASELINE
        // NOTE: Need only 620 samples of signal baseline for next section
        s1 = s2 = 0;
        for (i=20; i < 600; i++) {
          s1 += signal[i];
          s2 += (int) signal[i] * (int) signal[i];
        }
        s2 = sqrt(600*s2 - s1*s1);                          // rms * 800
        int slope = trap_fixed(signal, 20, 200, 200);       // baseline slope
        if (chan < 100) {                   // HG channel
          s2 /= 180.0;                      // baseline RMS, l 0.3 ADC units
          s1 /= 180.0;                      // baseline mean,l 0.3 ADC units
          slope = (slope + 3000*12)/12;     // baseline slope, 0.06 ADC units
        } else {                            // LG channel
          s2 /= 60.0;                       // baseline RMS, l 0.1 ADC units
          s1 /= 60.0;                       // baseline mean,l 0.1 ADC units
          slope = (slope + 3000*4)/4;       // baseline slope, 0.02 ADC units
        }
#else
        // NOTE: Need at least 820 samples of signal baseline for next section
        s1 = s2 = 0;
        for (i=20; i < 820; i++) {
          s1 += signal[i];
          s2 += (int) signal[i] * (int) signal[i];
        }
        s2 = sqrt(800*s2 - s1*s1);                          // rms * 800
        int slope = trap_fixed(signal, 20, 300, 200);       // baseline slope
        if (chan < 100) {                   // HG channel
          s2 /= 240.0;                      // baseline RMS, l 0.3 ADC units
          s1 /= 240.0;                      // baseline mean,l 0.3 ADC units
          slope = (slope + 3000*18)/18;     // baseline slope, 0.06 ADC units
        } else {                            // LG channel
          s2 /= 80.0;                       // baseline RMS, l 0.1 ADC units
          s1 /= 80.0;                       // baseline mean,l 0.1 ADC units
          slope = (slope + 3000*6)/6;       // baseline slope, 0.02 ADC units
        }
 #endif
        if (s1 > -500 && s1 < 1000) {  // histogram the three values (baseline, RMS, slope)
          int offset = 4000;                                    // low energy, not pulser
          if (e_trapmax >= 50 ||
              (pulser && chan%100 == 16 && e_trapmax > ptInfo.elo[chan]))  // special low-e case for ch 16
            offset = 0;                                           // high energy / pulser
          his[1600+chan][offset + (int) (s1+1000.5)]++;
          if (s2 < 500) his[1600+chan][offset + (int) (s2+2000.5)]++;
          if (slope > 2500 && slope < 3500) his[1600+chan][offset + slope]++;
        }
      }

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
      } else {
        dirty_sig = 0;
      }
      if (GRAN_CUT && granularity) dirty_sig += 32;

#ifndef DO_PSA    // ifdef DO_PSA, use e_adc and e_ctc, further below
      j = e_trapmax;
      if (chan < runInfo->nGe) {
        j = ChData[ievt]->e * Dets[chan].HGcalib[0] * 2.0;
      } else if (chan > 99 && chan < 100 + runInfo->nGe) {
        j = ChData[ievt]->e * Dets[chan-100].LGcalib[0] * 2.0;
      }
      if (dirty_sig) {
        his[400+chan][e_trapmax]++;
        if (j>0 && j<8192) his[800+chan][j]++;
      } else if (!pulser) {
        his[600+chan][e_trapmax]++;
        if (j>0 && j<8192) his[1000+chan][j]++;
      }
#endif

#ifdef DO_PSA
      /* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

      /* This section of signal processing is not requred for somee programs (e.g. noise,
       *      presort, and deadtime) so it is included only by defining the flag DO_PSA
       */

      double e_raw, e_adc, e_ctc, e_lamda, gain;
      float  fsignal[8192], drift, aovere, dcr, lamda;
      int    t0, t90, t100, ebin, a_e_good = 0, a_e_high = 0, dcr_good = 0, lamda_good = 0;
      if (chan < 100) {
        gain = Dets[chan].HGcalib[0];
      } else {
        gain = Dets[chan-100].LGcalib[0];
      }

      /* find t100 and t90*/
      t100 = 700;                 // FIXME? arbitrary 700?
      for (i = t100+1; i < 1500; i++)
        if (signal[t100] < signal[i]) t100 = i;
      /* get mean baseline value */
      int bl = 0;
      for (i=300; i<400; i++) bl += signal[i];
      bl /= 100;
      for (t90 = t100-1; t90 > 500; t90--)
        if ((signal[t90] - bl) <= (signal[t100] - bl)*19/20) break;

      if (t100   > 1300 ||      // important for cleaning, gets rid of pileup
          t90+760 > siglen) {   // dcr trap extends past end of signal
        if (VERBOSE || DEBUG)
          printf(">>> Error: chan %d signal at timestamp %lld is too late!\n", chan, time);
        dirty_sig += 64;  // FIXME; this is too late to be processed
        bad[ievt] = 7;
        // continue;   // FIXME
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
      int   tlo = t100+50, thi = siglen - 10;
      if (thi > tlo + 1500) thi = tlo + 1500;  // FIXME; check performance
      chisq = pz_fitter(fsignal, tlo, thi, chan, &PZI, &lamda1, &frac2, &lamda);
      if (chisq < 0.01 || chisq > 50.0) {       // fit failed, or bad fit chisq
        dirty_sig += 128;  // FIXME; this is too late to be processed
        bad[ievt] = 8;
        // continue;      // FIXME
      }
      lamda = (0.01 / PZI.tau[chan] - lamda) * 1.5e6 + 0.3;  // 0.3 fudge to get mean ~ 0

      /* do (required) PZ correction */
      if (PZ_fcorrect(fsignal+10, siglen-10, chan, &PZI))
        printf(" >>> PZ_fcorrect return error for chan %d\n", chan);

      /* get energy (e_raw and e_ctc) and effective drift time */
      e_ctc = get_CTC_energy(fsignal, siglen, chan, Dets, &PSA,
                             &t0, &e_adc, &e_raw, &drift, CTC.e_dt_slope[chan]);
      if (e_ctc < 0.001) {
        if (VERBOSE) printf("E_ctc = %.1f < 1 eV!\n", e_ctc);
        if (VERBOSE) printf("chan, t0, t100: %d %d %d; e, drift: %.2f %.2f\n",
                            chan, t0, t100, e_raw, drift);
        dirty_sig += 64;  // FIXME; this is too late to be processed
        bad[ievt] = 7;
        // continue;      // FIXME
      }
      ebin = 2.0 * e_ctc + 0.5;
      if (ebin < 0) ebin = 0;
      if (ebin > 8190) ebin = 8190;

      // do lamda-CT correction to get e_lamda
      e_lamda = (e_raw + lamda * CTC.e_lamda_slope[chan] * e_raw/6000.0) * CTC.e_lamda_gain[chan];
      e_lamda *= gain;

      lamda *= 8.0;                   // scaled to roughly match DCR at 2 MeV

      /* find A/E */
      if (e_ctc > 50 && t0 > 600 && t0 < 1300) {
        aovere = float_trap_max_range(fsignal, &tmax, PSA.a_e_rise[chan], 0, t0-20, t0+300);
        aovere *= PSA.a_e_factor[chan] / e_raw;
      } else {
        aovere = 0;
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
      if (siglen < 2450) {
        dcr = float_trap_fixed(fsignal, t90+50, 100, 500) / 25.0;
      } else {
        dcr = float_trap_fixed(fsignal, t90+50, 160, 800) / 40.0;
      }

      /* do drift-time and energy corrections to A/E, DCR, lamda */
      float dtc = drift*3.0*2614.0/e_adc;   // in (x)us, for DT- correction to A/E
      if (chan%100 == 50) dtc *= 3.0;       // special hack for detector 50; has especially strong variation
      float ec  = e_ctc/1000.0 - 2.0;       // in MeV, for E-correction to A/E
      aovere += PSA.ae_dt_slope[chan]    * dtc;
      aovere += PSA.ae_e_slope[chan]     * ec;
      dcr    -= PSA.dcr_dt_slope[chan]   * drift;
      lamda  -= PSA.lamda_dt_slope[chan] * drift;
      if (SUBTR_DCR_MEAN) {  // correct for residual INL; mean dcr and lamda vs. energy
        dcr   -= (float) dcr_mean[chan][(int) e_adc/2] / 10.0;
        lamda -= (float) dcr_mean[chan][4000 + (int) e_adc/2] / 10.0;
      }
      /* ------- end of signal shape processing; now do cuts and histogramming  -------- */

      /* ------- do A/E cut ------- */
      a_e_good = a_e_high = 0;
      // adjust for energy dependence of cut
      // first correct for series noise: assume 2*sigma cut
      // series noise contribution to A/E sigma = BL_RMS * sqrt(2*rise) * factor / E_raw
      if (AOE_CORRECT_NOISE)
        aovere -= (2.0 * PZI.bl_rms[chan] * sqrt(2.0 * (float) PSA.a_e_rise[chan]) *
                   PSA.a_e_factor[chan] * gain/1593.0 * (1.0 - 1593.0/e_ctc));        // 1593 keV = DEP
      // now deal with energy dependence of variation in A/E due to bremsstrahlung etc
      aovere += AOE_CORRECT_EDEP * (PSA.ae_pos[chan] - PSA.ae_cut[chan]) * (1.0 - 1593.0/e_ctc);  // FIXME: Add limit at low e_ctc

      if (aovere >= PSA.ae_cut[chan]) a_e_good = 1;
      if (aovere >= 2.0 * PSA.ae_pos[chan] - PSA.ae_cut[chan]) a_e_high = 1;
      if (aovere > 0 && aovere * 800.0/PSA.ae_pos[chan] < 2000) {
        his[1200+chan][(int) (aovere * 800.0/PSA.ae_pos[chan] + 0.5)]++;
        his[1200+chan][(int) ((aovere-PSA.ae_cut[chan]) * 800.0/PSA.ae_pos[chan] + 2800.5)]++;
      }

      /* ------- do DCR cut ------- */
      dcr_good = 0;
      if (dcr > -100 && e_ctc > 1000) {
        if (dcr < PSA.dcr_lim[chan]) dcr_good = 1;
        if (dcr < 300) {
          his[1400+chan][(int) (500.5  + dcr)]++;
          his[1400+chan][(int) (1000.5 + dcr - PSA.dcr_lim[chan])]++;
          if (!a_e_good) his[1400+chan][(int) (1500.5 + dcr - PSA.dcr_lim[chan])]++;
        }
        if (dcr < 3600)
          his[1400+chan][(int) (2000.5  + dcr/2.0)]++; // coarse scale to look for alpha distribution
      }
      /* ------- do lamda cut ------- */
      lamda_good = 0;
      if (lamda > -100 && e_ctc > 1000) {
        if (lamda < PSA.lamda_lim[chan]) lamda_good = 1;
        if (lamda < 300) {
          his[1400+chan][(int) (4500.5 + lamda)]++;
          his[1400+chan][(int) (5000.5 + lamda - PSA.lamda_lim[chan])]++;
          if (!a_e_good) his[1400+chan][(int) (5500.5 + lamda - PSA.lamda_lim[chan])]++;
        }
        if (lamda < 3600)
          his[1400+chan][(int) (6000.5 + lamda/2.0)]++; // coarse scale to look for alpha distribution
      }
      /* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

      if (0 && e_raw > 2000 && e_raw < 2010)
        printf("$$$$  ch %d e_raw = %.2f dirty = %d\n", chan, e_raw, dirty_sig);
      if (dirty_sig) {  // DIRTY signals
        if (e_adc >= 0 && e_adc < 8190) his[400+chan][(int) (e_adc + 0.5)]++;
        his[800+chan][ebin]++;
        his2[1200+chan][ebin]++;
        his2[99][1200+chan]++;
        if ((dirty_sig&48) == 0) { // exclude granularity and saturated signals
          his2[1400+chan][ebin]++;
          his2[99][1400+chan]++;
        }
        if (EVENTLIST && e_ctc > EVENTLIST && e_ctc < EVENTLISTMAX && chan < 100) {
          fprintf(f_evl, "%3d %7.1f %15lld %d%d%d.%d%d%d%d.%d%d%d%d",
                  chan, e_ctc, time, lamda_good, dcr_good, a_e_good,
                  DBIT(7), DBIT(6), DBIT(5), DBIT(4), DBIT(3), DBIT(2), DBIT(1), DBIT(0));
          if (chan > 99)
            fprintf(f_evl, " %s LG Q0  %8.1f %6.1f %6.1f   %d %3d %3d  %8.2f %6.2f %6.2f\n",
                    Dets[chan-100].StrName, aovere, dcr, lamda, runInfo->runNumber, t90-t0, t100-t90,
                    aovere - PSA.ae_cut[chan], dcr - PSA.dcr_lim[chan], lamda - PSA.lamda_lim[chan]);
          else
            fprintf(f_evl, " %s HG Q0  %8.1f %6.1f %6.1f   %d %3d %3d  %8.2f %6.2f %6.2f\n",
                    Dets[chan].StrName, aovere, dcr, lamda, runInfo->runNumber, t90-t0, t100-t90,
                    aovere - PSA.ae_cut[chan], dcr - PSA.dcr_lim[chan], lamda - PSA.lamda_lim[chan]);
        }

      } else if (!bad[ievt] && !pulser) {     // CLEAN non-pulser signals

        if (e_adc >= 0 && e_adc < 8190) his[600+chan][(int) (e_adc + 0.5)]++;
        his[1000+chan][ebin]++;

        // ----------------------------------------------------------------
        // histogram drift time (t90 - t0)
        if (t90 >= t0 && t90 < t0+450 && e_ctc > 1000 && e_ctc < 3000) {
          his[2000+chan][5000+t90-t0]++;
          if (a_e_good) his[2000+chan][5500+t90-t0]++;
          if (dcr_good) his[2000+chan][6000+t90-t0]++;
          if (!dcr_good) his[2000+chan][6500+t90-t0]++;
          if (a_e_good && !dcr_good) his[2000+chan][7000+t90-t0]++;
          if (a_e_good && dcr_good) his[2000+chan][7500+t90-t0]++;
        }
        // ----------------------------------------------------------------
        /* histogram final results for clean signals */
        his2[chan][ebin]++;            // no cuts
        his2[99][chan]++;
        // if (a_e_high) {             // A/E cut // TEMPORARY
        if (a_e_good) {                // A/E cut
          his2[200+chan][ebin]++;
          his2[99][200+chan]++;
        }
        if (dcr_good) {                // DCR cut
          his2[400+chan][ebin]++;
          his2[99][400+chan]++;
        }
        if (lamda_good) {              // lamda cut
          his2[600+chan][ebin]++;
          his2[99][600+chan]++;
        }
        if (a_e_good && dcr_good && lamda_good) { // all cuts
          his2[800+chan][ebin]++;
          his2[99][800+chan]++;
          int hg_lg_only = 1;
          int jevt;
          for (jevt = 0; jevt < nChData; jevt++) {
            if (chan == ChData[jevt]->chan + 100 ||
                chan == ChData[jevt]->chan - 100) hg_lg_only = 0;
          }
          if (hg_lg_only) {         // all cuts and only the HG or LG ch is hit
            his2[1000+chan][ebin]++;
            his2[99][1000+chan]++;
          }
          // ----------------------------------------------------------------

          if (WRITE_SIGS && chan < runInfo->nGe) {
            // && e_ctc > 2605 && e_ctc < 2625 &&  // TEMPORARY
            //  (aovere > 0.96*PSA.ae_cut[chan] || PSA.ae_cut[chan] < 100)) {  // TEMPORARY
            if (MAKE_2D)
              fprintf(f_2d, "%6.1f, %6.1f, c%2.2d\n", e_ctc, aovere * 800.0/PSA.ae_pos[chan], chan);

            /* write out good 2614-keV signals */
            /* find mean baseline value        */
            int t50;
            bl = 0.0;
            for (i=300; i<700; i++) bl += fsignal[i];
            bl /= 400.0;
            // subtract initial baseline and convert to shorts
            for (i=0; i<siglen; i++) signal[i] = lrintf(fsignal[i] - bl);
            /* find t50 */
            for (t50 = 900; t50 < 1700; t50++)
              if ((signal[t50]) >= e_adc/2) break;
            for (; t50 > 500; t50--)
              if ((signal[t50]) <= e_adc/2) break;
            /* find max value just after rise */
            s1 = 0;
            for (i = t50+21; i < t50+56; i++) s1 += fsignal[i] - bl;
            s1 /= 35.0;
            if (s1 < 0.98*e_adc) continue;  // throw away pileup and slow signals
            // align t50 to t=200
            j = t50 - 200;
            for (i=0; i < 256; i++) signal[i] = signal[i+j];
            for (i=256; i < 512; i++) signal[i] = 0;
            // take derivative
            k = 0;
            for (i=0; i < 4; i++) k += signal[i+5] - signal[i+1];
            signal[256] = k/2;
            for (i=1; i < 247; i++) {
              k += signal[i+8] - 2*signal[i+4] + signal[i];
              signal[i+256] = k/2;
            }
            // add some extra information at the start
            signal[1] = (short) e_ctc;
            signal[2] = chan;
            signal[3] = a_e_good + a_e_high - 1;
            signal[4] = aovere * 800.0/PSA.ae_pos[chan];
            signal[5] = aovere;
            fwrite(signal, 1024, 1, f_sig);
          }
        }

        if (EVENTLIST && e_ctc > EVENTLIST && e_ctc < EVENTLISTMAX && chan < 100) {
          fprintf(f_evl, "%3d %7.1f %15lld %d%d%d.%d%d%d%d.%d%d%d%d",
                  chan, e_ctc, time, lamda_good, dcr_good, a_e_good,
                  DBIT(7), DBIT(6), DBIT(5), DBIT(4), DBIT(3), DBIT(2), DBIT(1), DBIT(0));
          if (chan > 99)
            fprintf(f_evl, " %s LG Q%d  %8.1f %6.1f %6.1f   %d %3d %3d  %8.2f %6.2f %6.2f\n",
                    Dets[chan-100].StrName, a_e_good + 2*dcr_good + 4*lamda_good,
                    aovere, dcr, lamda, runInfo->runNumber, t90-t0, t100-t90,
                    aovere - PSA.ae_cut[chan], dcr - PSA.dcr_lim[chan], lamda - PSA.lamda_lim[chan]);
          else
            fprintf(f_evl, " %s HG Q%d  %8.1f %6.1f %6.1f   %d %3d %3d  %8.2f %6.2f %6.2f\n",
                    Dets[chan].StrName, a_e_good + 2*dcr_good + 4*lamda_good,
                    aovere, dcr, lamda, runInfo->runNumber, t90-t0, t100-t90,
                    aovere - PSA.ae_cut[chan], dcr - PSA.dcr_lim[chan], lamda - PSA.lamda_lim[chan]);
        }
      }
      // ----------------------------------------------------------------
#endif

      /* write out good (or selected) signals */
#ifndef QUIET
      int wsig_len = siglen;
      if (0) {        // TURNED OFF FOR PRESORT
        // optionally subtract initial baseline
        j = (signal[20]+signal[21]+signal[22]+signal[23]+signal[24])/5;
        for (i=0; i<wsig_len; i++) signal[i] -= j;
      }
      if (1) {
        // optionally add some extra information at the start
        signal  -= 8; // 8 extra words (shorts)
        wsig_len += 8;
        for (i=0; i<8; i++) signal[i] = 0;
        signal[1] = (short) (e_onbd/on_bd_rise[idet]);
        signal[2] = (short) e_trapmax;
        signal[4] = (short) ((head2[8]>>11)&0xf)*100;
        signal[5] = (short) dirty_sig;
        signal[6] = (short) chan;
      }

      if (sigout_channum == chan && e_trapmax > 0 &&
          ((sigout_ptag == pulser &&
                        (!pulser || e_trapmax > ptInfo.elo[chan]) &&
                        (sigout_ctag < 0 ||
                         (sigout_ctag == 0 && dirty_sig) ||
                         (sigout_ctag == 1 && !dirty_sig))) ||
           (sigout_ptag < 0 && ((sigout_ctag == 0 && dirty_sig) ||
                                (sigout_ctag == 1 && !dirty_sig)))))
          write_sig(signal, wsig_len, chan, f_out);

      if (sigout_channum == -1 && !pulser && // chan<99 &&
          ((doDC && !dirty_sig && e_trapmax > 2) ||
           (!doDC && e_trapmax > 4))) {
        write_sig(signal, wsig_len, chan, f_out);
      }

      /* write out dirty signals */
      if (sigout_channum >= -1 && doDC &&
          dirty_sig && (dirty_sig&48) == 0 && // exclude granularity and saturated signals
          //e_ctc > 1000 &&
          e_trapmax > 4) {
        write_sig(signal, wsig_len, chan, bf_out);
      }

      // test inl correction
      if (DEBUG && !pulser && e_trapmax >= 10 && chan < 100+runInfo->nGe) {
        float fsignal[8192];
        if (inl_correct(Dets, runInfo, signal, fsignal, siglen, chan)) {
          printf(" >>> inl_correct return error for chan %d\n", chan);
        } else {
          float e1 = float_trap_max(fsignal, &tmax, 401, 200)/401.0;
          printf(" << chan %3d e0, e1, diff = %7.2f %7.2f %5.2f\n",
                 chan, ChData[ievt]->e, e1, ChData[ievt]->e - e1);
        }
      }
#endif
      if (DEBUG) { printf("done...\n"); fflush(stdout); }
    }
    /* ------------- end of Ge detector data processing -------------- */

    /* deal with bad events identified so far */
    if (bad[ievt]) {  // bad event, not just dirty/cleaned
      badevts++;
      badevts_bytype[bad[ievt]]++;
      return 2;
    }

  } // END of loop over channel-events

  goodevts++;
#ifndef QUIET
  if (totevts % 20000 < (totevts-nChData) % 20000) {
    printf(" Processed %d built evts, %d ch-evts, %d bad ch-evts (%d%%)  ",
           builtevts, totevts, badevts, (badevts*100)/totevts);
    for (i=1; i<6; i++) printf(" %d", badevts_bytype[i]);
    printf(" %6d pulser evts\n", pulserevts);
  }
#endif
  if (pulser) return 1;
  return  0;
} /* eventprocess() */

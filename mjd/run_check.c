/* read through channel-events in data file
   look for events that have the right energy to be pulser events
   make sure they are coming at the correct rate
   extract the ch-by-ch offsets of the on-board traps
   suggest threshold values that minimize the positive threshold
     while avoiding small negative thresholds, which result in dead time
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MJDSort.h"

/* ---------- dataId types we want to process ---------- */
char  *GoodDataTypes[] =
  {
   "ORCV830DecoderForEvent",
   "ORCV830DecoderForPolledRead",
   "ORCAEN792DecoderForQdc",
   "ORCAEN792NDecoderForQdc",
   "ORGretina4MWaveformDecoder",
   "ORGretina4AWaveformDecoder",
   "ORRunDecoderForRun",
   "-"
  };
//  "ORMJDPreAmpDecoderForAdc",
//  "ORScriptDecoderForRecord",
//  "ORScriptDecoderForState",
//  "ORiSegHVCardDecoderForHV",

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define VERBOSE   0

int peak_find_small(int *his, int lo, int hi) {  // find highest-energy peak between chs lo and hi
  // extremely simple algorithm, requires zero background counts; could be much improved
  int e, emax = hi;
  for (e=hi-1; e>lo; e--) {
    if (his[e] > his[emax]) emax = e;
    if (his[e] <= 0 && his[emax] > 2) return emax;
  }
  return 0;
}

/*  ------------------------------------------------------------ */

int eventbuild(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int   i, j, k, crate=0, slot=0, board, evlen, ch, chan, idet;
  int   dataId[32], goodDataId[32], idNum, board_type;
  int   dataIdGM=0, dataIdGA=0, dataIdRun=0;
  int   his[200][8192] = {{0}};
  int   e_onbd, e_offline, e_trapmax, de, step;
  unsigned int  head[2], evtdat[20000];
  unsigned short *head2;
  short  *signal, sigin[2048], sigu[2048];
  double s1;
  FILE   *f_out;

  int totevts = 0, out_evts = 0;
  int start_time = 0, tmax; //  end_time = 0;
  long long int time=0, min_time=1, max_time=0;
  int module_lu[NCRATES+1][21]; // lookup table to map VME crate&slot into module IDs
  int det_lu[NBDS][16];         // lookup table to map module&chan into detector IDs
  int chan_lu[NBDS][16];        // lookup table to map module&chan into parameter IDs
  int e_thresh[200];            // Ge energy thresholds for event building
  int on_bd_rise[200], on_bd_flat[200];  // on-board trap rise and flat times (10ns)
  int on_bd_change[200];        // flag for consecutive wrong rise-times for on-board trap
  PTag  ptInfo;                 // pulser energy and delta-time info for tagging

  int    e, dt=0, traptime[200], traparea[200];
  float  pos, area, fwhm, obt_offset[200];
  char   spname[64];

  
  if (runInfo->analysisPass < 0) return 0;

  /* initialize */
  if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;

  /* get pulser period and amplitude information */
  if (!pulser_tag_init(Dets, runInfo, &ptInfo)) {
    printf("\n Sorry, I need the pulser period and amplitude information!\n");
    return -1;
  }

  /* set energy thresholds and trap times for analysis */
  for (i=0; i<200; i++) e_thresh[i] = 50;
  e_thresh[16] = 10;
  e_thresh[116] = 3;
  for (i=0; i<200; i++) {
    on_bd_rise[i] = 401;
    on_bd_flat[i] = 201;
    on_bd_change[i] = 0;
    obt_offset[i] = 0;
  }

  /* identify the data types that we want to decode */
  idNum = runInfo->idNum;
  dataIdGM = runInfo->dataIdGM;
  dataIdGA = runInfo->dataIdGA;
  for (i=0; i<idNum; i++) {
    dataId[i] = runInfo->dataId[i];
    goodDataId[i] = 0;
    for (j=0; GoodDataTypes[j][0]=='O'; j++)
      if (strstr(runInfo->decoder[i], GoodDataTypes[j])) {
        goodDataId[i] = 1;
        if (strstr(runInfo->decoder[i], "ORGretina4MWaveformDecoder"))  dataIdGM    = dataId[i];
        if (strstr(runInfo->decoder[i], "ORGretina4AWaveformDecoder"))  dataIdGA    = dataId[i];
        if (strstr(runInfo->decoder[i], "ORRunDecoderForRun"))          dataIdRun   = dataId[i];
      }
  }
  goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below

  /* start loop over reading events from input file
     ============================================== */

  while (fread(head, sizeof(head), 1, f_in) == 1) {
    board_type = head[0] >> 18;
    evlen = (head[0] & 0x3ffff);

    if (board_type == 0) {  // a new runfile header! file must be corrupt?
      printf("\n >>>> ERROR: DataID = 0; found a second file header??"
             " Ending scan of this file!\n"
             " >>>> head = %8.8x %8.8x  evlen = %d\n", head[0], head[1], evlen);
      break;
    }

    /* see if the event ID matches a known type */
    for (j=0; j<idNum; j++) {
      if (dataId[j] == board_type) break;  // found a good dataId
    }

    /* --------------------------------------------------------- */
    if (j >= idNum) {  // the dataID type is unknown (not found in the file header)
      evlen -= 2;
      if (evlen > 10000) {
        printf("\n >>>> ERROR: Event length too long??\n"
               " >>>> This file is probably corruped, ending scan!\n");
        break;
      }
      while (evlen > 10000) {
        fread(evtdat, 4, 10000, f_in);
        evlen -= 10000;
      }
      if (evlen > 0) fread(evtdat, 4, evlen, f_in);
     continue;
    }
    /* --------------------------------------------------------- */

    /* if we don't want to decode this type of data, just skip forward in the file */
    if (!goodDataId[j]) {
      // fseek(f_in, 4*(evlen-2), SEEK_CUR);
      fread(evtdat, 4, evlen-2, f_in);
      continue;
    }

    /* --------------- ORRunDecoderForRun ---------------- */
    if (board_type == dataIdRun) {
      fread(evtdat, 8, 1, f_in);
      if (head[1] & 0x21) {
        printf("------- START Run %d at %d", evtdat[0], evtdat[1]);
        if (head[1] & 0x4) printf("; remote");
        if (head[1] & 0x2) printf("; quick-start");
        printf(" -------\n");
        if (start_time == 0) start_time = evtdat[1];
      } else if (head[1] & 0x8) {
        if (VERBOSE)
          printf(" -- %8.8x Heartbeat step %d at %d\n", head[1], evtdat[0], evtdat[1]);
      } else if ((head[1] & 0x29) == 0) {
        printf("------- END run %d at %d;  duration = %d s -------\n",
               evtdat[0], evtdat[1], evtdat[1]-start_time);
        // end_time = evtdat[1];
      } else {
        printf("***** ORRunDecoderForRun %8.8x %d %d\n",
               head[1], evtdat[0], evtdat[1]);
      }
      continue;
    }
    /* ------------- end ORRunDecoderForRun end -------------- */

    slot  = (head[1] >> 16) & 0x1f;
    crate = (head[1] >> 21) & 0xf;
    board = -1;
    
    if (crate < 0 || crate > NCRATES ||
        slot  < 0 || slot > 20 ||
        (board = module_lu[crate][slot]) < 0) {
      printf("ERROR: Illegal VME crate or slot number %d %2d; %8.8x %8.8x %d\n",
             crate, slot, head[0], head[1], board);
      if (board >= 0) printf("  --> %s\n", runInfo->decoder[board]);
      //fwrite(head, sizeof(head), 1, f_out_bad);
      if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) break;
      //fwrite(evtdat, sizeof(int), evlen, f_out_bad);
      continue;
    }
    if (VERBOSE && totevts < 50)
      printf("len, cr, slot, board: %d %d %d %d\n", evlen, crate, slot, board);

    /* ========== read in the rest of the event data ========== */
    if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
      printf("  No more data...\n");
      break;
    }

    totevts++;
    if (board_type != dataIdGM && board_type != dataIdGA) continue;   // (not) GRETINA digitizer

    /* --------------- Gretina4M digitizer ---------------- */

    ch = (evtdat[1] & 0xf);   // for Gretina4M data
    if (board < 0 || ch > 9) continue;
    chan = chan_lu[board][ch];
    idet = det_lu[board][ch];
    if (idet >= runInfo->nGe) continue;
    /* extract 48-bit time stamp */
    time = (evtdat[3] & 0xffff);
    time = time << 32 | evtdat[2];
    /* if time is less than 0, then this is probably a bad event */
    if (time < 0) continue;
    if (min_time == 1 || min_time > time) min_time = time;
    if (max_time < time) max_time = time;

    /* ------------------------------------------------------------------ */
    /* process the event */
    signal = (short *) evtdat + 28;
    if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
        2020 == decompress_signal((unsigned short *)signal, sigu, 2*(evlen - 2) - 28)) {
      signal = sigu;
      evlen = 1026;
    }
    head2  = (unsigned short *) evtdat;
    e_onbd = (((head2[8] & 0x1ff) << 16) | head2[7]);  // on-board energy from event header

    /* --- look for transition to 4x presumming ------- */
    if (chan <= 99) {
      k = Dets[chan].HGPrerecnt + Dets[chan].HGPostrecnt;  // expected location of presumming start
    } else {
      k = Dets[chan-100].LGPrerecnt + Dets[chan-100].LGPostrecnt;
    }
    if (k > 10 && k < 2008) { // correct for presumming
      // FIXME: replace hard-coded factor of 4 with data from header
      step = k;
      j = (signal[step-5] + signal[step-4] + signal[step-3] + signal[step-2])/4;
      if (j > 10 || j < -10) {
        /* --- Find transition to 4x presumming ------- */
        for (step=k; step<k+2; step++) {
          if (j > 0 && signal[step] > 2*j) break;
          if (j < 0 && signal[step] < 2*j) break;
        }
        if (step == k+2) {
          step = k+1;
          //printf("Hmmm... step not found in %d - %d in chan %3d, counts %d\n",
          //       step-1, step, chan, j);
        }
      }

      /* --- uncompress --- */
      /* --- make a copy of the compressed part --- */
      for (i=0; i<(2018-step)/4; i++) sigin[i] = signal[i+step];
      /* --- and uncompress it -------------------- */
      for (i=0; i<(2018-step)/4; i++) {         // NOTE: This is wrong for signal < 0
        for (j=0; j<4; j++)
          signal[step + 4*i + j] = sigin[i]/4;  // distribute sums over 4 bins each
        for (j=0; j<sigin[i]%4; j++)
          signal[step + 4*i + j]++;             // and make sure sum is right
      }
    }
    /* --- Done with 4x presumming ------- */

    e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);
    e_trapmax = e_offline/on_bd_rise[idet];
    if (e_offline < e_thresh[chan] * on_bd_rise[idet]) continue;  // sub-energy-threshold; skip event
    out_evts++;

    /* ----- check for bad trap rise time; 401 -> 465 or vice versa ----- */
    if (1 && e_trapmax > 200 && e_trapmax < 6000) {
      s1 = (float) e_onbd / (float) e_offline;
      if (on_bd_rise[idet] == 401 && s1 > 1.13 && s1 < 1.17) {
        // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
        if (on_bd_change[idet]++ > 1) {
          on_bd_rise[idet] = 465;
          e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);
          printf(" ****> Det %d on-board rise time: %d -> %d\n", idet, 401, 465);
          on_bd_change[idet] = 0;
        }
      } else if (on_bd_rise[idet] == 465 && 1.0/s1 > 1.13 && 1.0/s1 < 1.17) {
        // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
        if (on_bd_change[idet]++ > 1) {
          on_bd_rise[idet] = 401;
          e_offline = trap_max(signal, &tmax, on_bd_rise[idet], on_bd_flat[idet]);
          printf(" ****> Det %d on-board rise time: %d -> %d\n", idet, 465, 401);
          on_bd_change[idet] = 0;
        }
      } else {
        on_bd_change[idet] = 0;
      }
    }
    s1 = 100.0 * (e_onbd - e_offline) / (double) on_bd_rise[idet];
    de = (int) (s1 + 1000.5);    // onboard - offline difference, 0.01 ADC units
    if (de > 0 && de < 2000) his[chan][de]++;  // 0.01 ADC units
    de = de/20 + 3950;           // units are now 0.2 ADC cts
    if (de < 2000 || de > 6000) de = 6000;
    his[chan][de]++;             // 0.2 ADC units

    /* now look to see if the energy of the event matches a pulser event
       first make sure we're using a trap rise time of 401 samples
     */
    if (on_bd_rise[idet] != 401) {
      e_offline = trap_max(signal, &tmax, 401, 201);
      e_trapmax = e_offline/401;
    }
    /* now compare with pulser amplitude limits */
    if (e_trapmax >= ptInfo.elo[chan] &&
        e_trapmax  <= ptInfo.ehi[chan] &&
        ptInfo.pdt[chan] > 0) {
      // signal meets pulser energy gate
      his[chan][6000+tmax]++;    // trap-max time plus 6000
    }

    // if (out_evts > 500) break;

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  dt = (max_time - min_time)/100000000;  // GRETINA card run time in s
  printf("\n %d events in, %d events out, digitizer time range = %d s (%d min) (\n",
         totevts, out_evts,  dt, (dt+30)/60);

  if (totevts < 200) return 0; // too few events to process

  /*  ======  extract on-board trap offsets and thresholds, and tmax  ======  */
  printf("\n Here an offset of -99.99 indicates insufficient valid data!\n"
         "                               On-bd            True             "
         "               On-bd            True\n"
         "    Detector    rise  HGChan   offset        threshold  Cnts Tmax"
         "     LGChan    offset        threshold  Cnts Tmax\n");
  for (chan=0; chan < 100+runInfo->nGe; chan++) {
    obt_offset[chan] = -99.99;
    traptime[chan] = traparea[chan] = 0;
    if (chan >= runInfo->nGe && chan < 100) continue;
    if (chan < runInfo->nGe && !Dets[chan].HGChEnabled) continue;
    if (chan > 99 && chan < 100+runInfo->nGe && !Dets[chan-100].LGChEnabled) continue;

    /* look for tmax peak */
    if ((e = peak_find_small(his[chan], 6001, 8000))) { // peak found
      if ((pos = autopeak(his[chan], e, &area, &fwhm)) != 0) {
        traptime[chan] = pos - 6000;
        traparea[chan] = area;
      }
    }
    /* look for peak around bin 1000 (0.01 ADC) */
    if ((e = peak_find_small(his[chan], 1, 1999))) { // peak found
      if ((pos = autopeak(his[chan], e, &area, &fwhm)) == 0) pos = e;
      pos = 0.01 * (pos - 1000.0);
    } else {  /* look for peak above bin 4000 (0.1 ADC)*/
      if (!(e = peak_find_small(his[chan], 2001, 7999)))  continue; // no peak found
      if ((pos = autopeak(his[chan], e, &area, &fwhm)) == 0) pos = e;
      pos = 0.2 * (pos - 4000.0);
    }
    obt_offset[chan] = pos;
    if (chan%100 == 16 && traparea[chan] < area) traparea[chan] = area;
  }
  for (chan=0; chan < runInfo->nGe; chan++) {
    if (!Dets[chan].HGChEnabled) continue;
    printf("%8s %s %4d %6d", Dets[chan].DetName, Dets[chan].StrName, on_bd_rise[chan], chan);
    printf("%9.2f %7.2f ->%7.2f %4d %4d",
           obt_offset[chan],
           Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan],
           Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan],
           traparea[chan], traptime[chan]);
    printf("%10d %9.2f  %7.2f ->%7.2f %4d %4d\n",
           chan+100, obt_offset[chan+100],
           Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan],
           Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan+100],
           traparea[chan+100], traptime[chan+100]);
  }

  printf("\nChannels with low pulser rate or bad timimg:\n"
         "   tMax should be between 350 and 550\n"
         "cr,sl,ch DetID  Rate(%%)  tMax   DetPos      Notes\n");
  for (i=0; i < runInfo->nGe; i++) {
    if (Dets[i].HGChEnabled) {
      chan = i;
      if (traparea[chan] < dt/16 ||
          traptime[chan] < 350 || traptime[chan] > 600) {
        printf("%3d,%2.2d,%d %6d %6.0f %7d  %8s HG",
               Dets[i].crate, Dets[i].slot, Dets[i].chanHi, chan,
               838.8 * (float) traparea[chan] / (float) dt, traptime[chan], Dets[i].StrName);
        if (i==16) printf(" [bad tMax expected]");
        if (i==42) printf(" [no pulser expected]");
        printf("\n");
      }
    }
    if (Dets[i].LGChEnabled) {
      chan = i+100;
      if (traparea[chan] < dt/16 ||
          traptime[chan] < 350 || traptime[chan] > 600) {
        printf("%3d,%2.2d,%d %6d %6.0f %7d  %8s LG",
               Dets[i].crate, Dets[i].slot, Dets[i].chanLo, chan,
               838.8 * (float) traparea[chan] / (float) dt, traptime[chan], Dets[i].StrName);
        if (i==16) printf(" [bad tMax expected]");
        if (i==42) printf(" [no pulser expected]");
        printf("\n");
      }
    }
  }
  

  /* ===== make a file listing channels for which the threshold finder should be re-run ===== */
  if (!(f_out = fopen("TFChannels.txt", "w"))) return 0;
  j = 0;
  for (i=0; i < runInfo->nGe; i++) {
    if (i == 42) continue;  // no pulser for this channel
    if (Dets[i].HGChEnabled &&
        obt_offset[i] < -0.3 * Dets[i].HGTrapThreshold / (float) on_bd_rise[i]) {
      if (obt_offset[i] < -99) {
        fprintf(f_out, "%d,%d,%d    %s HG; current thresh = %d;"
                " No data; raise thresh to %d and repeat\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanHi,
                Dets[i].StrName, Dets[i].HGTrapThreshold, Dets[i].HGTrapThreshold+200);
      } else {
        k = 0.8 * Dets[i].HGTrapThreshold - 1.8 * obt_offset[i] * (float) on_bd_rise[i];
        fprintf(f_out, "%d,%d,%d    %s HG; current thresh = %d; suggested thresh = %d\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanHi,
                Dets[i].StrName, Dets[i].HGTrapThreshold, k);
      }
      j++;
    }
    if (Dets[i].LGChEnabled &&
        obt_offset[i+100] < -0.3 * Dets[i].LGTrapThreshold / (float) on_bd_rise[i]) {
      if (obt_offset[i+100] < -99) {
        fprintf(f_out, "%d,%d,%d    %s LG; current thresh = %d;"
                " No data; raise thresh to %d and repeat\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanLo,
                Dets[i].StrName, Dets[i].LGTrapThreshold, Dets[i].LGTrapThreshold+200);
      } else {
        k = 0.8 * Dets[i].LGTrapThreshold - 1.8 * obt_offset[i+100] * (float) on_bd_rise[i];
        fprintf(f_out, "%d,%d,%d    %s LG; current thresh = %d; suggested thresh = %d\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanLo,
                Dets[i].StrName, Dets[i].LGTrapThreshold, k);
      }
      j++;
    }
  }
  fclose(f_out);

  if (j)
    printf("\n %d channels require the Threshold Finder; see TFChannels.txt\n", j);
  else
    printf("\n No channels require the Threshold Finder.\n");
  
  /* if needed, write out diagnostic histograms */
  if (0) {
    f_out = fopen("rc.rms", "w");
    for (i=0; i<200; i++) {
      sprintf(spname, "%d run_check of run %d, chan %d", i, runInfo->runNumber, i);
      write_his(his[i], 8192, i, spname, f_out);
    }
    fclose(f_out);
  }

  return  out_evts;
} /* eventbuild() */

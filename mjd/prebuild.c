/*
   MJD prebuilder for runs where one GRETINA card has lost sync
   David Radford Feb 2017

   Run this program twice... once to find the problem, then again to fix it.

   Type 1 run: Find time offset for boards out-of-sync
      read through channel-events in data file
      look for events that have the right energy to be pulser events
      look for other boards with the same CC pulser and check time differences
      if there is a significant difference then one of them has lost sync,
          so save the time difference to run?????.pb
   Type 2 run: Use the time differences from run?????.pb or default.pb to fix the
               times in the file
      read through channel-events in data file and re-write to DS* output file
      for all places where the bad board ID occurs,
          change the time in the event before rewriting
      can use the same pb_dt.input file for multiple runs

   Build with presort.c and up_util.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MJDSort.h"

/* ---------- dataId types we want to process ---------- */
char  *GoodDataTypes[] =
  {
    //   "ORCV830DecoderForEvent",
    //   "ORCV830DecoderForPolledRead",
    //   "ORCAEN792DecoderForQdc",
    //   "ORCAEN792NDecoderForQdc",
    "ORGretina4MWaveformDecoder",
    "ORGretina4AWaveformDecoder",
    //   "ORRunDecoderForRun",
    "-"
  };
//  "ORMJDPreAmpDecoderForAdc",
//  "ORScriptDecoderForRecord",
//  "ORScriptDecoderForState",
//  "ORiSegHVCardDecoderForHV",

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define VERBOSE   0

/*  ------------------------------------------------------------ */

int eventprescan(FILE *f_in, FILE *ps_f_out, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int   i, j, k, crate=0, slot=0, board, evlen, ch, chan, idet;
  int   dataId[32], goodDataId[32], idNum, board_type;
  int   dataIdGM = 0, dataIdGA = 0;
  int   his[200][8192] = {{0}};
  int   e_offline, e_trapmax, step;
  unsigned int  head[2], evtdat[20000];
  short  *signal, sigin[2048];
  FILE   *f_out=0, *f;

  int totevts = 0, out_evts = 0;
  int tmax;
  long long int time=0, min_time=1, max_time=0;
  int module_lu[NCRATES+1][21]; // lookup table to map VME crate&slot into module IDs
  int det_lu[NBDS][16];         // lookup table to map module&chan into detector IDs
  int chan_lu[NBDS][16];        // lookup table to map module&chan into parameter IDs
  int e_thresh[200];            // Ge energy thresholds for event building
  int on_bd_rise[200], on_bd_flat[200];  // on-board trap rise and flat times (10ns)
  PTag  ptInfo;                 // pulser energy and delta-time info for tagging

  int    lastboard = -1, lastCC = -1, no_time_changes = 1;
  int    pmod[NBDS+1] = {0}, ndt[NBDS] = {0};
  long long int dtmod[NBDS][100] = {{0}}, lasttime=0, dt=0, t_offset[NBDS] = {0};
  char   spname[64], line[128], fname[256];
  double dts;

  
  if (runInfo->analysisPass < 0) return 0;

  /* initialize */
  if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;

  /* get pulser period and amplitude information */
  if (!pulser_tag_init(Dets, runInfo, &ptInfo)) {
    printf("\n Sorry, I need the pulser period and amplitude information!\n");
    return -1;
  }
  /* if a prebuild time offset file exists, read the time offsets from it */
  sprintf(fname, "run%d.pb", runInfo->runNumber);
  if ((f = fopen(fname, "r")) ||
      (f = fopen("default.pb", "r"))) {
    printf("\n");
    while (fgets(line, sizeof(line), f) &&
           sscanf(line, "%d %lld", &j, &dt) == 2) {
      t_offset[j] = dt;
      printf(" >>>  From %s or default.pb: board %d prebuild time offset = %lld\n",
             fname, j, dt);
      no_time_changes = 0;
    }
    printf("\n");
    fclose(f);
  }    

  /* set energy thresholds and trap times for analysis */
  for (i=0; i<200; i++) e_thresh[i] = 50;
  e_thresh[16] = 10;
  e_thresh[116] = 3;
  for (i=0; i<200; i++) {
    on_bd_rise[i] = 401;
    on_bd_flat[i] = 201;
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
        if (strstr(runInfo->decoder[i], "ORGretina4MWaveformDecoder"))  dataIdGM   = dataId[i];
        if (strstr(runInfo->decoder[i], "ORGretina4AWaveformDecoder"))  dataIdGA   = dataId[i];
      }
  }
  goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below

  /* start loop over reading events from input file
     ============================================== */

  while (fread(head, sizeof(head), 1, f_in) == 1) {
    fwrite(head, sizeof(head), 1, ps_f_out);  // copy input data to output file
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
        fwrite(evtdat, 4, 10000, ps_f_out);  // copy input data to output file
        evlen -= 10000;
      }
      if (evlen > 0) {
        fread(evtdat, 4, evlen, f_in);
        fwrite(evtdat, 4, evlen, ps_f_out);  // copy input data to output file
      }
     continue;
    }
    /* --------------------------------------------------------- */

    /* if we don't want to decode this type of data, just skip forward in the file */
    if (!goodDataId[j]) {
      fread(evtdat, 4, evlen-2, f_in);
      fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file
      continue;
    }

    slot  = (head[1] >> 16) & 0x1f;
    crate = (head[1] >> 21) & 0xf;
    board = -1;
    
    if (crate < 0 || crate > NCRATES ||
        slot  < 0 || slot > 20 ||
        (board = module_lu[crate][slot]) < 0) {
      printf("ERROR: Illegal VME crate or slot number %d %2d; %8.8x %8.8x %d\n",
             crate, slot, head[0], head[1], board);
      if (board >= 0) printf("  --> %s\n", runInfo->decoder[board]);
      if (fread(evtdat, 4, evlen-2, f_in) != evlen-2) break;
      fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file
      continue;
    }
    if (VERBOSE && totevts < 50)
      printf("len, cr, slot, board: %d %d %d %d\n", evlen, crate, slot, board);

    /* ========== read in the rest of the event data ========== */
    if (fread(evtdat, 4, evlen-2, f_in) != evlen-2) {
      printf("  No more data...\n");
      break;
    }

    totevts++;
    if (board_type != dataIdGM && board_type != dataIdGA) {
      fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file
      continue;   // (not) GRETINA digitizer
    }

    /* --------------- Gretina4M digitizer ---------------- */

    ch = (evtdat[1] & 0xf);   // for Gretina4M data
    if (board < 0 || ch > 9) continue;
    chan = chan_lu[board][ch];
    idet = det_lu[board][ch];
    if (idet >= runInfo->nGe || idet == 16) {
      fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file
      continue;
    }
    /* extract 48-bit time stamp */
    time = (evtdat[3] & 0xffff);
    time = time << 32 | evtdat[2];
    if (t_offset[board]) {
      time += t_offset[board];
      evtdat[3] = (evtdat[3] & 0xffff0000) | ((time >> 32) & 0xffff);
      evtdat[2] = time & 0xffffffff;
      dt = (evtdat[3] & 0xffff);
      dt = dt << 32 | evtdat[2];
      if (dt != time) printf("\n\nERROR! dt = %lld time = %lld\n\n", dt, time);
    }
    fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file

    /* if time is less than 0, then this is probably a bad event */
    if (time < 0) continue;
    if (min_time == 1 || min_time > time) min_time = time;
    if (max_time < time) max_time = time;

    /* ------------------------------------------------------------------ */
    /* process the event */
    signal = (short *) evtdat + 28;
    if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
        2020 == decompress_signal((unsigned short *)signal, sigin, 2*(evlen - 2) - 28)) {
      memcpy(signal, sigin, 4040);
      evlen = 1026;
    }

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

    /* look to see if the energy of the event matches a pulser event */
    e_offline = trap_max(signal, &tmax, 401, 201);
    e_trapmax = e_offline/401;
    if (e_trapmax >= ptInfo.elo[chan] &&
        e_trapmax  <= ptInfo.ehi[chan] &&
        ptInfo.pdt[chan] > 0) {
      // signal meets pulser energy gate
      if (pmod[NBDS] < 8000 &&
          (pmod[NBDS] == 0 || pmod[NBDS] != pmod[board])) {
        pmod[board] = ++pmod[NBDS];
        his[board][pmod[NBDS]] = 10 + Dets[idet].CCnum;
        his[NBDS][pmod[NBDS]]  = 10 + 20*Dets[idet].CCnum + board;
      }
      if (board < 20 &&
          lastCC == Dets[idet].CCnum &&
          lastboard != board) {
        // previous observed pulser signal was in a different board but on the same CC
        dt = time - lasttime;
        if (ndt[board] < 99) dtmod[board][ndt[board]++] = dt;
        if (ndt[lastboard] < 99) dtmod[lastboard][ndt[lastboard]++] = -dt;
        k = 4000 + dt;
        if (k > 0 && k < 8000) his[board+40][k]++;
        k = 4000 - dt;
        if (k > 0 && k < 8000) his[lastboard+40][k]++;
        k = 4000 + dt/2000;
        if (k > 0 && k < 8000) his[board+60][k]++;
        k = 4000 - dt/2000;
        if (k > 0 && k < 8000) his[lastboard+60][k]++;
        k = 4000 + dt/4000000;
        if (k > 0 && k < 8000) his[board+80][k]++;
        k = 4000 - dt/4000000;
        if (k > 0 && k < 8000) his[lastboard+80][k]++;
        k = 4000 + dt/8000000000;
        if (k > 0 && k < 8000) his[board+100][k]++;
        k = 4000 - dt/8000000000;
        //          6308959166189
        if (k > 0 && k < 8000) his[lastboard+100][k]++;
      } else if (board < 16 && lastboard < 16 &&
                 lastCC != Dets[idet].CCnum &&
                 lastboard != board) {
        // previous observed pulser signal was in a different board and a different CC
        dt = time - lasttime;
        k = 4000 + dt;
        if (k > 0 && k < 8000) his[board+120][k]++;
        k = 4000 + dt/2000;
        if (k > 0 && k < 8000) his[board+140][k]++;
        k = 4000 + dt/4000000;
        if (k > 0 && k < 8000) his[board+160][k]++;
        k = 4000 + dt/8000000000;
        if (k > 0 && k < 8000) his[board+180][k]++;
      }
    }
    lastboard = board;
    lastCC = Dets[idet].CCnum;
    lasttime = time;

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */
  printf("Ended reading at file position %ld\n\n", ftell(f_in));

  k = 0;
  sprintf(fname, "run%d.pb", runInfo->runNumber);
  for (i=0; i<16; i++) {
    if (ndt[i] > 8) {
      dts = 0;
      for (j=4; j<ndt[i]; j++) dts += dtmod[i][j];
      dts /= (ndt[i] - 4);
      printf("Board %2d dt = %lld %lld %lld %lld %lld %lld; mean = %.0lf\n",
             i, dtmod[i][4], dtmod[i][5], dtmod[i][6], dtmod[i][7], dtmod[i][8],
             dtmod[i][ndt[i]-1], dts);
      if (dts < -4000) {
        if (k++ == 0) {
          f_out = fopen(fname, "w");
          for (j=0; j<NBDS; j++)
            if (t_offset[j]) fprintf(f_out, "%2d %lld\n", j, t_offset[j]);
        }
        fprintf(f_out, "%2d %.0f\n", i, -dts);
      }
    }
  }
  if (k) {
    fclose(f_out);
    printf("\n\n Created new version of time-offset file %s\n"
           " ***************** %d new case(s) found with lost sync *****************\n"
           " You should re-run this code to fix it and make a new DS%d file\n\n",
           fname, k, runInfo->runNumber);
  } else if (no_time_changes) {   // no lost lock found
    printf("\n\n --------------- No lost sync found ---------------\n"
           " You should probably delete the output file DS%d\n\n",
           runInfo->runNumber);
  }

  /* if needed, write out diagnostic histograms */
  if (0) {
    f_out = fopen("pb.rms", "w");
    for (i=0; i<200; i++) {
      sprintf(spname, "%d run_check of run %d, chan %d", i, runInfo->runNumber, i);
      write_his(his[i], 8192, i, spname, f_out);
    }
    fclose(f_out);
  }

  return  out_evts;
} /* eventprescan() */

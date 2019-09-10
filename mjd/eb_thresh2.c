/*
   eventbuild.c

   event builder code for MJD
   - up to 16 digitizer modules for the Ge detectors, plus
   - QDCs etc for the veto

   David Radford   Nov 2016
*/
/* read through channel-events in data file and build global events
   which are then passed out to eventprocess() for analysis.
   returns: -1 on error
            otherwise the actual number of channel-events processed
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MJDSort.h"

#define  INPUT_THRESH_FNAME "thresholds_for_DAQ.txt"

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
  int   dataIdGM=0, dataIdRun=0;
  int   his[200][8192] = {{0}};
  int   e_onbd, e_offline, e_trapmax, de, step;
  unsigned int  head[2], evtdat[20000];
  unsigned short *head2;
  short  *signal, sigin[2048], sigu[2048];
  double s1;
  FILE   *f_in2, *f_out;

  int totevts = 0, out_evts = 0;
  int start_time = 0;  // end_time = 0;
  long long int time=0;
  int module_lu[NCRATES+1][21]; // lookup table to map VME crate&slot into module IDs
  int det_lu[NBDS][16];         // lookup table to map module&chan into detector IDs
  int chan_lu[NBDS][16];        // lookup table to map module&chan into parameter IDs
  int e_thresh[200];            // Ge energy thresholds for event building
  int on_bd_rise[200], on_bd_flat[200];  // on-board trap rise and flat times (10ns)
  int on_bd_change[200];        // flag for consecutive wrong rise-times for on-board trap

  int    e, repeat=0, load_thresholds = 0, argn;
  float  pos, area, fwhm, t, obt_offset[200], daq_thresh[200], suggested_thresh[200], tfactor=1.0;
  char   line[256], command[40000], *c, fname[256];


  if (runInfo->analysisPass < 0) return 0;

  /* initialize */
  if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
  strncpy(fname, INPUT_THRESH_FNAME, sizeof(fname));

  /* look at command line arguments to see if threshold file and/or factor are specified */
  for (argn=1; argn < runInfo->argc-1; argn++) {
    if (strstr(runInfo->argv[argn], "-t"))      // specifies input thresholds file
      strncpy(fname, runInfo->argv[++argn], sizeof(fname));
    if (strstr(runInfo->argv[argn], "-f")) {    // specifies input thresholds factor
      if (sscanf(runInfo->argv[++argn], "%f", &tfactor) != 1) {
        fprintf(stderr, "\nERROR:  Failed to read threshold factor from arguments %s %s\n",
                runInfo->argv[argn-1], runInfo->argv[argn]);
        return 0;
      }
    }
  }

  /* read in requested absolute thresholds */
  for (i=0; i<200; i++) daq_thresh[i] = suggested_thresh[i] = 3.0f;
  if (!(f_in2 = fopen(fname, "r"))) {
    fprintf(stderr, "\nERROR: Cannot open input threshold file %s\n  Aborting!\n\n", fname);
    return -1;
  }
  while (fgets(line, sizeof(line), f_in2)) {
    if (line[0] == '#' || strlen(line) < 4) continue;
    if (sscanf(line, "%d %f", &i, &t) != 2 || i < 0 || i > 199 || i%100 > 57) {
      fprintf(stderr," ERROR: reading of input threshold  failed!\n line: %s\n", line);
      fclose(f_in2);
      return -1;
    }
    daq_thresh[i] = t*tfactor;
    if (i < 100) daq_thresh[i+100] = t*tfactor;
  }
  fclose(f_in2);
  printf("Thresholds read from file %s\n", fname);

  /* set separate offline energy thresholds and trap times for analysis here */
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
  for (i=0; i<idNum; i++) {
    dataId[i] = runInfo->dataId[i];
    goodDataId[i] = 0;
    for (j=0; GoodDataTypes[j][0]=='O'; j++)
      if (strstr(runInfo->decoder[i], GoodDataTypes[j])) {
        goodDataId[i] = 1;
        if (strstr(runInfo->decoder[i], "ORGretina4MWaveformDecoder"))  dataIdGM    = dataId[i];
        if (strstr(runInfo->decoder[i], "ORRunDecoderForRun"))          dataIdRun   = dataId[i];
      }
  }
  goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below

  /* start loop over reading events from input file
     ============================================== */

  while (1) {

    if (fread(head, sizeof(head), 1, f_in) != 1) break;
    board_type = head[0] >> 18;
    evlen = (head[0] & 0x3ffff);

    if (board_type == 0) {  // a new runfile header! file must be corrupt?
      fprintf(stderr, "\n >>>> ERROR: DataID = 0; found a second file header??"
              " Ending scan of this file!\n"
              " >>>> head = %8.8x %8.8x  evlen = %d\n", head[0], head[1], evlen);
      break;
    }

    /* see if the event ID matches a known type */
    for (j=0; j<idNum; j++) {
      if (dataId[j] == board_type) break;  // found a good dataId
    }

    /* --------------------------------------------------------- */
    if (j >= idNum) {  // the dataID type is not found in the file header
      evlen -= 2;
      if (evlen > 10000) {
        fprintf(stderr, "\n >>>> ERROR: Event length too long??\n"
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
      fread(evtdat, 4, evlen-2, f_in); // fread lets us use stdin as raw data input
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
    if (board_type != dataIdGM) continue;   // (not) GRETINA digitizer

    /* --------------- Gretina4M digitizer ---------------- */
    /* if time is less than 0, then this is probably a bad event */
    if (time < 0) continue;

    ch = (evtdat[1] & 0xf);   // for Gretina4M data
    if (board < 0 || ch > 9) continue;
    chan = chan_lu[board][ch];
    idet = det_lu[board][ch];
    if (idet >= runInfo->nGe) continue;

    /* ------------------------------------------------------------------ */
    /* process the event */
    signal = (short *) evtdat + 28;
    if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
        2020 == decompress_signal((unsigned short *)signal, sigu, 2*(evlen - 2) - 28)) {
      signal = sigu;
      evlen = 1026;
    }
    head2  = (unsigned short *) evtdat;
    e_onbd = (((head2[8] & 0x1ff) << 16) | head2[7]);

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

    e_offline = trap_max(signal, &j, on_bd_rise[idet], on_bd_flat[idet]);
    e_trapmax = e_offline/on_bd_rise[idet];
    if (e_offline < e_thresh[chan] * on_bd_rise[idet]) continue;  // sub-threshold
    out_evts++;

    /* ----- check for bad trap rise time in the onboard value; 401 -> 465 or vice versa ----- */
    if (1 && e_trapmax > 200 && e_trapmax < 6000) {
      s1 = (float) e_onbd / (float) e_offline;
      if (on_bd_rise[idet] == 401 && s1 > 1.13 && s1 < 1.17) {
        // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
        if (on_bd_change[idet]++ > 1) {
          on_bd_rise[idet] = 465;
          e_offline = trap_max(signal, &j, on_bd_rise[idet], on_bd_flat[idet]);
          printf(" ****> Det %d on-board rise time: %d -> %d\n", idet, 401, 465);
          on_bd_change[idet] = 0;
        }
      } else if (on_bd_rise[idet] == 465 && 1.0/s1 > 1.13 && 1.0/s1 < 1.17) {
        // on_bd_change[] counts number of consecutive pulses that indicate wrong rise time
        if (on_bd_change[idet]++ > 1) {
          on_bd_rise[idet] = 401;
          e_offline = trap_max(signal, &j, on_bd_rise[idet], on_bd_flat[idet]);
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
    if (de < 2000 || de > 8000) de = 8001;
    his[chan][de]++;             // 0.2 ADC units

    // if (out_evts > 500) break;

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  printf("\n %d events in, %d events out\n", totevts, out_evts);

  if (totevts < 200) return 0; // too few events to process

  /*  ======  extract on-board trap offsets and thresholds, and tmax  ======  */
  for (chan=0; chan < 100+runInfo->nGe; chan++) {
    obt_offset[chan] = -99.99;
    if (chan >= runInfo->nGe && chan < 100) continue;
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
  }

  /* report results, including new suggested threshold */
  printf("\n Here an offset of -99.99 indicates insufficient valid data!\n"
         "                               On-bd            True  "
         "                               On-bd            True\n"
         "    Detector    rise  HGChan   offset        threshold     req adjusted"
         "    LGChan    offset        threshold     req adjusted\n");
  for (chan=0; chan < runInfo->nGe; chan++) {
    if (!Dets[chan].HGChEnabled) continue;

    suggested_thresh[chan] = daq_thresh[chan] + obt_offset[chan];
    if (suggested_thresh[chan] < 0.8f * daq_thresh[chan] - obt_offset[chan])
      suggested_thresh[chan] = 0.8f * daq_thresh[chan] - obt_offset[chan];
    suggested_thresh[chan+100] = daq_thresh[chan+100] + obt_offset[chan+100];
    if (suggested_thresh[chan+100] < 0.8f * daq_thresh[chan+100] - obt_offset[chan+100])
      suggested_thresh[chan+100] = 0.8f * daq_thresh[chan+100] - obt_offset[chan+100];

    printf("%8s %s %4d %6d", Dets[chan].DetName, Dets[chan].StrName, on_bd_rise[chan], chan);
    printf("%9.2f %7.2f ->%7.2f %7.2f %7.2f",
           obt_offset[chan],
           Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan],
           Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan],
           daq_thresh[chan], suggested_thresh[chan]);
    printf("%10d %9.2f %7.2f ->%7.2f %7.2f %7.2f\n",
           chan+100, obt_offset[chan+100],
           Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan],
           Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan+100],
           daq_thresh[chan+100], suggested_thresh[chan+100]);
  }

  /* ===== make a file listing the thresholds to be loaded ===== */
  if (!(f_out = fopen("TChannels.txt", "w"))) return 0;
  load_thresholds = 0;
  for (i=0; i < runInfo->nGe; i++) {
    if (i == 42) continue;  // no pulser for this channel
    if (Dets[i].HGChEnabled) {
      if (obt_offset[i] < -99) {
        fprintf(f_out, "%d,%d,%d, %d ;   %s HG; current thresh = %d;"
                " No data; raise thresh to %d and repeat\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanHi,
                (int) (suggested_thresh[i] * (float) on_bd_rise[i] + 0.5),
                Dets[i].StrName, Dets[i].HGTrapThreshold, Dets[i].HGTrapThreshold+200);
      } else {
        fprintf(f_out, "%d,%d,%d, %d ;   %s HG; current thresh = %d\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanHi,
                (int) (suggested_thresh[i] * (float) on_bd_rise[i] + 0.5),
                Dets[i].StrName, Dets[i].HGTrapThreshold);
      }
      load_thresholds++;
    }
    if (Dets[i].LGChEnabled) {
      if (obt_offset[i+100] < -99) {
        fprintf(f_out, "%d,%d,%d %d ;    %s LG; current thresh = %d;"
                " No data; raise thresh to %d and repeat\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanLo,
                (int) (suggested_thresh[i+100] * (float) on_bd_rise[i] + 0.5),
                Dets[i].StrName, Dets[i].LGTrapThreshold, Dets[i].LGTrapThreshold+200);
      } else {
        fprintf(f_out, "%d,%d,%d, %d ;   %s LG; current thresh = %d\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanLo,
                (int) (suggested_thresh[i+100] * (float) on_bd_rise[i] + 0.5),
                Dets[i].StrName, Dets[i].LGTrapThreshold);
      }
      load_thresholds++;
    }
  }
  fclose(f_out);

  /* ===== save current onboard thresholds in a new file; thresholds_from_DAQ.save" ===== */
  if (!(f_out = fopen("thresholds_from_DAQ.save", "w"))) return 0;
  fprintf(f_out, "# ID ADCthresh rise\n");
  for (i=0; i < 100+runInfo->nGe; i++) {
    if (i > runInfo->nGe && i < 100) continue;
    if (i < 100 && Dets[i].HGChEnabled && obt_offset[i] > -99)
      fprintf(f_out, "%3d %6.2f %10d   %s\n",
              i, Dets[i].HGTrapThreshold / (float) on_bd_rise[i] - obt_offset[i],
              on_bd_rise[i], Dets[i].StrName);
    else if (i > 99 && Dets[i-100].LGChEnabled && obt_offset[i] > -99)
      fprintf(f_out, "%3d %6.2f %10d   %s\n",
              i, Dets[i-100].LGTrapThreshold/ (float) on_bd_rise[i-100] - obt_offset[i],
              on_bd_rise[i-100], Dets[i-100].StrName);
  }
  fclose(f_out);

  if (!load_thresholds)
    printf("\n No channels require the Threshold Finder.\n");
  else {
    printf("\n %d channels require the Threshold Finder; see TFChannels.txt\n", load_thresholds);

    /* ----------  create command string to load new thresholds into Orca ---------- */
    repeat = 0;                    // do we need to repeat this with the new thresholds?
    bzero(command, sizeof(command));
    c = command;
    sprintf(c, "[RunControl stopRun];\n"); c = command+strlen(command);
    for (k = 0; k < runInfo->nGD; k++) {  // loop through all digitizer modules
      load_thresholds = 0;                              // count number of thresholds to modify
      for (i=0; i < runInfo->nGe; i++) {
        if (i == 42) continue;            // no pulser for channel 42
        if (Dets[i].crate != runInfo->GDcrate[k] ||
            Dets[i].slot != runInfo->GDslot[k]) continue;   // detector is not on this digitizer

        if (Dets[i].HGChEnabled) {
          if (obt_offset[i] < -99) {
            sprintf(c, "[ORGretina4MModel,%d,%d setTrapThreshold:%d withValue:%d];\n",
                    Dets[i].crate, Dets[i].slot, Dets[i].chanHi, Dets[i].HGTrapThreshold + 200);
            c = command+strlen(command);
            load_thresholds++;
            repeat = 1;
          } else {
            sprintf(c, "[ORGretina4MModel,%d,%d setTrapThreshold:%d withValue:%d];\n",
                    Dets[i].crate, Dets[i].slot, Dets[i].chanHi,
                    (int) (suggested_thresh[i] * (float) on_bd_rise[i] + 0.5));
            c = command+strlen(command);
            load_thresholds++;
          }
        }
        if (Dets[i].LGChEnabled) {
          if (obt_offset[i+100] < -99) {
            sprintf(c, "[ORGretina4MModel,%d,%d setTrapThreshold:%d withValue:%d];\n",
                    Dets[i].crate, Dets[i].slot, Dets[i].chanLo, Dets[i].LGTrapThreshold + 200);
            c = command+strlen(command);
            load_thresholds++;
            repeat = 1;
          } else {
            sprintf(c, "[ORGretina4MModel,%d,%d setTrapThreshold:%d withValue:%d];\n",
                    Dets[i].crate, Dets[i].slot, Dets[i].chanLo,
                    (int) (suggested_thresh[i+100] * (float) on_bd_rise[i] + 0.5));
            c = command+strlen(command);
            load_thresholds++;
          }
        }
      }   // for (i=0; i < runInfo->nGe; i++) {
      if (load_thresholds > 0) {  // modified some thresholds in this card
        sprintf(c, "[ORGretina4MModel,%d,%d loadThresholds];\n",
                runInfo->GDcrate[k], runInfo->GDslot[k]);
        c = command+strlen(command);
      }
    }   // for (k = 0; k < runInfo->nGD; k++) {
    sprintf(c, "[RunControl startRun];\n"); c = command+strlen(command);

    printf("\n -------------------------\n%s", command);
  }

  if (repeat)
    fprintf(stderr,
            "\n ****** NOTE: There were enabled channels with no data! ******\n"
            "              Suggest repeat running of %s\n", runInfo->argv[0]);

  return  out_evts;
} /* eventbuild() */

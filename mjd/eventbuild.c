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
/*
 *  NOTE: define compile-time option PRESORT for version
 *         where we are making presorted data subsets
 *         (was presortbuild.c)
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

#define VERBOSE 0
#define DEBUG   0
#define DO_NOT_COMPRESS 0
// #define DETSEL (chan == 29 || chan == 30)

/* set BUILD_VETO to 1 to include veto data in the built events sent to eventprocess */
#define BUILD_VETO 0
#define NEVTS  50000               /* no. of events to use in each board's event buffer
                                      (should only need ~ 50 in low-rate situations) */
#define RESOLVING_TIME        400  /* up to 4 us delay allowed within one built event */ 
#define BUFFER_TIME (20*100000000) /* up to 20s delay allowed before forcing event building */
/*  ------------------------------------------------------------ */

int flush_buffers(MJDetInfo *Dets, MJRunInfo *runInfo, BdEvent *modBuf[NBDS], int nModBuf,
                  int nevts[NBDS], int iptr[NBDS], int order[NBDS], int *nBdsAvail,
#ifdef PRESORT
                  int *out_evts, int *built_evts, FILE *ps_f_out,
#else
                  int *out_evts, int *built_evts,
#endif
                  int *badevts, int his[16][8192], long long int lastBuiltTime) {

  int     i, j, k, n, nChData;
  BdEvent *ChData[200];        // pointers for constructing a built-event
  long long int dt, min_time;
#ifdef PRESORT
  unsigned short sigc[2048];
  short         *signal;
#endif

  //printf("fb_in\n"); fflush(stdout);
  /* nBdsAvail is the number of data channels for which
     we have at least one stored event */
  while (*nBdsAvail > 0) {
    n = 0;
    for (i=0; i<nModBuf; i++) {
      if (nevts[i] <= 0) continue;
      if (n == 0 ||
          modBuf[order[n-1]][iptr[order[n-1]]].time <= modBuf[i][iptr[i]].time) {
        order[n++] = i;
      } else {
        for (j=0; j<n; j++) {
          if (modBuf[order[j]][iptr[order[j]]].time > modBuf[i][iptr[i]].time) {
            for (k=n; k>j; k--) order[k] = order[k-1];
            order[j] = i;
            n++;
            break;
          }
        }
      }
    }
    if (n != *nBdsAvail) {
      printf("\nAckk! In flush_buffers(), n = %d != nBdsAvail = %d !\n\n",
             n, *nBdsAvail);
      return -1;
    }
    min_time = modBuf[order[0]][iptr[order[0]]].time;  // earliest time in all module buffers

    nChData = 0;
    for (i=0; i<n; i++) {   // assemble and process board data
      j = order[i];
      while (//modBuf[j][iptr[j]].time > min_time - 20 &&  // FIXME
             modBuf[j][iptr[j]].time - min_time < RESOLVING_TIME) {
        // extend event time by later times? if so, uncomment next line...
        // if (min_time < modBuf[j][iptr[j]].time) min_time = modBuf[j][iptr[j]].time;
        if (modBuf[j][iptr[j]].time - min_time < -1000000000) { // 10 s; clocks must have been reset?
          printf("Error: Large jump back in timestamps. Clocks reset?\n");
          break;
        }
        his[3][210+modBuf[j][iptr[j]].crate]++;
        if (nChData >= 200) {
          printf("ERROR: More than 200 sub-events in built event!\n\n");
          printf("   min_time = %12lld\n\n Ordering:\n", min_time);
          for (k=0; k<n; k++)
            printf("%12lld\n", modBuf[order[k]][iptr[order[k]]].time);
          printf("\n");
          for (k=0; k<200; k+=4)
            printf("%12lld %12lld %12lld %12lld\n", ChData[k]->time,
                   ChData[k+1]->time, ChData[k+2]->time, ChData[k+3]->time);
          //return -1;
          break;
        }
        (*out_evts)++;
        if (DEBUG) {
          printf("%8d out; i,j,time,mintime: %2d %2d %lld %lld\n",
                 *out_evts, i, j, modBuf[j][iptr[j]].time, min_time);
          fflush(stdout);
        }
        ChData[nChData++] = &modBuf[j][iptr[j]];  // insert channel-event into built-event
        if (++iptr[j] >= NEVTS) iptr[j] -= NEVTS; // and remove from module buffer
        if (--nevts[j] <= 0) {
          (*nBdsAvail)--;
          break;
        }
      }
    }
    /* now analyze/process the built-event, and write to skim file for PRESORT */
    /* first histogram time since last built event */
    dt = (ChData[0]->time - lastBuiltTime) / 100;  // us
    lastBuiltTime = ChData[0]->time;
    if (dt > -1000 && dt < 7000) his[11][1000+dt]++;
    dt /= 1000; // ms
    if (dt > -1000 && dt < 7000) his[12][1000+dt]++;

    /* eventprocess return values:
       -1 on error,
       1 on pulser event,
       2 on bad event (bad = cannot be processed; not dirty)
       0 otherwise (good event, can be clean or dirty)
    */
    (*built_evts)++;
    //printf("ep_in\n"); fflush(stdout);
    if ((k = eventprocess(Dets, runInfo, nChData, ChData)) < 0) return k;
    //printf("ep_out\n"); fflush(stdout);
#ifdef PRESORT
    if (k < 1) {    // good event, not a pulser; can be dirty
      for (i = 0; i < nChData; i++) {
        if (!DO_NOT_COMPRESS &&
            (ChData[i]->orca_type == runInfo->dataIdGM ||
             ChData[i]->orca_type == runInfo->dataIdGA) &&   // GRETINA digitizer
            ChData[i]->evlen == 1026) {  // uncompressed; compress the signal
          signal = (short *) ChData[i]->evbuf + 32;
          int m = compress_signal(signal, sigc, 2020);
          memcpy(signal, sigc, m*2);
          ChData[i]->evlen = 16 + m/2;
          ChData[i]->evbuf[0] = (ChData[i]->evbuf[0] & 0xfffc0000) + ChData[i]->evlen;
        }
        fwrite(ChData[i]->evbuf, sizeof(int), ChData[i]->evlen, ps_f_out);
      }
    }
#endif
    if (k > 1) (*badevts)++;  // CHECKME
  }
  printf("Flushing finished...\n"); // fflush(stdout);

  return 0;
}

/*  ------------------------------------------------------------ */

#ifdef PRESORT
int eventprescan(FILE *f_in, FILE *ps_f_out, MJDetInfo *Dets, MJRunInfo *runInfo) {
#else
int eventbuild(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo) {
#endif

  int    i, j, k, n, crate=0, slot=0, board, evlen, ch, chan=0, nmod, ee, e_onbd;
  unsigned int  head[2], evtdat[20000];
  static int    dataId[32], goodDataId[32], idNum, board_type;
  static int    dataIdGM=0, dataIdGA=0, dataIdRun=0;
  static int    dataIdC830e=0, dataIdC830p=0, dataIdC792=0, dataIdC792N=0;
  static int    totevts, badevts, subthreshevts, out_evts, veto_evts, built_evts, recordID;

  long long int time=0, time2=0, time3, min_time = -1, prev_ch_time=0;
  int           min_time_board=-1, istep1, istep2, step, step2, siglen;
  static int    ntOutOfOrder = 0, oldCrate=0;
  static int    veto_count, veto_count2, sub_veto_count, veto_error_count=0;
  static int    start_time, start_time_veto, start_time_Ge;
  static int    end_time, end_time_veto, end_time_Ge;
  static long long int dt, dtsum, lastBuiltTime=0, oldTime=0;
  // static long long int lastChTime[NBDS][10];

  static BdEvent  *modBuf[NBDS];  // one buffer for each module, holds up to NEVT chan-events
  static short    *sigBuf[NBDS][NEVTS];  // signal buffers corresponding to modBuf
  static int      nModBuf = 0;
  static int      nevts[NBDS], iptr[NBDS], order[NBDS], nBdsAvail = 0; // book-keeping for modBuf[]
  static int      presum[200] = {0};     // expected presum step (or zero for no presum)
  short    ucsig[32*2048];      // temporary storage for expanded signal (up to a factor of 32)
  BdEvent  *ChData[200];        // pointers for constructing a built-event
  int      nChData;             // number of channel-events in the built-event

  static int module_lu[NCRATES+1][21]; // lookup table to map VME crate&slot into module IDs
  static int det_lu[NBDS][16];         // lookup table to map module&chan into detector IDs
  static int chan_lu[NBDS][16];        // lookup table to map module&chan into parameter IDs
  static int e_thresh[200];            // Ge energy thresholds for event building
  static int his[16][8192] = {{0}};    // histograms for diagnostics/development

  static int first = 1, good_run = 1, good_run_list[16000], n_good_run_list = 0;
  static int first_runNumber = 0, current_runNumber = 0, run_evts_count = 0, run_gretina_evts_count = 0;
  char     *c, spname[64], line[256], fname[256];
  short    *signal, mkrs[8192];
  FILE     *f_in2, *f_out;
  static int           last_ch_energy[200] = {0};
  static long long int last_ch_time[200] = {0}, last_bd_time[20] = {0};
  static long long int t_offset[200] = {0}, last_gretina_time = -1;
  int  cts=0, erro=0, chan_k;


#ifdef PRESORT
  if (runInfo->analysisPass < 0) {
    /* end-of-run cleanup */
    if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                      &nBdsAvail, &out_evts, &built_evts, ps_f_out,
                      &badevts, his, lastBuiltTime) < 0)
      return -1;

    /* change filename stored in runInfo for the benefit of ep_finalize and plot titles, etc */
    /* sprintf(runInfo->filename, "DS%d (runs %d to %d)",
       first_runNumber, first_runNumber, runInfo->runNumber); */
    printf(" Flushing eventprocess...\n");
    sprintf(runInfo->filename, "DS%d", first_runNumber);
    if ((k = eventprocess(Dets, runInfo, -99, ChData)) < 0) return k;
    printf("\n %d events in, %d events out, %d events built, %d * 3 veto events\n",
           totevts, out_evts, built_evts, veto_evts);
    printf(" %d bad events, %d sub-threshold events\n", badevts, subthreshevts);
    return  out_evts;
  }
#else
  /* some versions of the event builder allow for a special call
     to trigger end-of-run cleanup; this isn't one of them, so ignore */
  if (runInfo->analysisPass < 0) return 0;
  //return eventprocess(Dets, runInfo, -99, ChData);
#endif
  
  /* initialize */
  current_runNumber = runInfo->runNumber;

  if (DS0) {
    runInfo->nCC = 8;
    Dets[ 0].CCnum = 1;
    Dets[ 1].CCnum = 1;
    Dets[ 2].CCnum = 1;
    Dets[ 3].CCnum = 2;
    Dets[ 4].CCnum = 0;
    Dets[ 5].CCnum = 3;
    Dets[ 6].CCnum = 4;
    Dets[ 7].CCnum = 1;
    Dets[ 8].CCnum = 0;
    Dets[ 9].CCnum = 0;
    Dets[10].CCnum = 3;
    Dets[11].CCnum = 3;
    Dets[12].CCnum = 3;
    Dets[13].CCnum = 4;
    Dets[14].CCnum = 0;
    Dets[15].CCnum = 4;
    Dets[16].CCnum = 4;
    Dets[17].CCnum = 2;
    Dets[18].CCnum = 2;
    Dets[19].CCnum = 2;
    Dets[20].CCnum = 1;
    Dets[21].CCnum = 1;
    Dets[22].CCnum = 1;
    Dets[23].CCnum = 5;
    Dets[24].CCnum = 0;
    Dets[25].CCnum = 2;
    Dets[26].CCnum = 2;
    Dets[27].CCnum = 2;
    Dets[28].CCnum = 0;
  }

  /* initialize some values */
  if ((nmod = ep_init(Dets, runInfo, module_lu, det_lu, chan_lu)) < 0) return -1;
#ifdef PRESORT
  if (nModBuf > 0 && nmod+1 != nModBuf) {
    fprintf(stderr, "ERROR: Number of VME boards has changed from %d to %d!\n\n", nModBuf, nmod+1);
    return -1;
  }
#else
  if (nModBuf > 0 && nmod != nModBuf) {
    fprintf(stderr, "ERROR: Number of VME boards has changed from %d to %d!\n\n", nModBuf, nmod);
    return -1;
  }
#endif

  /* look for flag indicating that we shoudl preselct only some (good) runs from a presorted subset */
  if (first) {
    for (i=1; i<runInfo->argc; i++) {
      if (runInfo->argv[i][0] == '-' &&
          strstr(runInfo->argv[i], "g")) {  // -g flag defined in argument list; use good_runs.input
        if (!(f_in2 = fopen("good_runs.input", "r"))) {
          printf("ERROR: -g flag set but file good_runs.input does not exist!\n\n");
          exit(-1);
        }
        n_good_run_list = 0;
        while (fgets(line, sizeof(line), f_in2) &&
               (n = sscanf(line, "%d %d", &j, &k)) > 0) {
          if (n > 1) {
            while (j <= k) good_run_list[n_good_run_list++] = j++;
          } else {
            good_run_list[n_good_run_list++] = j;
          }
        }
        fclose(f_in2);
        printf("NOTE: %d good run numbetrs read from file good_runs.input\n", n_good_run_list);
        break;
      }
    }
  }
  
  if (first || runInfo->analysisPass > 0) {
    first = 0;
    runInfo->firstRunNumber = first_runNumber = current_runNumber = runInfo->runNumber;
    i = first_runNumber; // to silence warning when PRESORT is not defined
    totevts = badevts = subthreshevts = out_evts = veto_evts = built_evts = dtsum = recordID = 0;
    min_time = veto_count = veto_count2 = start_time_veto = start_time_Ge = -1;
    start_time =  end_time = end_time_veto = end_time_Ge = 0;
    sub_veto_count = 3;

    for (i=0; i<NBDS; i++) {
      nevts[i] = iptr[i] = 0;
      // for (j=0; j<10; j++) lastChTime[i][j]=0;
    }

    /* if a prebuild time offset file exists, read the time offsets from it */
    sprintf(fname, "run%d.pb", runInfo->runNumber);
    if ((f_in2 = fopen(fname, "r")) ||
        (f_in2 = fopen("default.pb", "r"))) {
      printf("\n");
      while (fgets(line, sizeof(line), f_in2) &&
           sscanf(line, "%d %lld", &j, &dt) == 2) {
        t_offset[j] = dt;
        printf(" >>>  From %s or default.pb: board %d prebuild time offset = %lld\n",
               fname, j, dt);
      }
      printf("\n");
      fclose(f_in2);
    }

    /* read energy thresholds from file thresholds.input*/
    for (i=0; i<200; i++) e_thresh[i] = 0;
    if ((f_in2 = fopen(E_THRESHOLDS_FILENAME, "r"))) {
      int m = 0;
      while (fgets(line, sizeof(line), f_in2)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d %d", &j, &k, &n) != 3) break;
        if (j > k || j < 0 || k > 99) break;
        for (i=j; i<=k; i++) {
          if (i >=0 && i < 100) {
            e_thresh[i] = n;                // high-gain threshold
            e_thresh[100+i] = n/3.4 + 4.0;  // low-gain threshold
            m++;
          }
        }
      }
      fclose(f_in2);
      printf(" %d energy thresholds defined from file %s\n", 2*m, E_THRESHOLDS_FILENAME);
    } else {
      printf("No thresholds defined for event builder, using zeroes\n");
    }

    if (nModBuf == 0) {  // only do this section if this is a first pass through the data
      /* malloc module data buffers for event building */
      for (i=0; i<runInfo->nGD; i++) {
        if (VERBOSE)
          printf("%2d    GeDig %2d of %2d crate %d slot %d\n",
                 nModBuf, i, runInfo->nGD, runInfo->GDcrate[i], runInfo->GDslot[i]);
        /* malloc space to the corresponding module event buffer */
        if (!(modBuf[nModBuf++] = malloc(NEVTS * sizeof(BdEvent)))) {
          printf("Malloc failed for modBuf[%d][%d]\n\n", nModBuf-1, NEVTS);
          return -1;
        }
        /* decide if there is presumming; if so, then we also need to
           malloc space for sigBuf, to hold expanded signals */
        k = 2008;
        chan_k = -1;
        for (ch=0; ch<10; ch++) {
          chan = chan_lu[module_lu[runInfo->GDcrate[i]][runInfo->GDslot[i]]][ch];
          if (chan >= 0 && chan < runInfo->nGe) {
            if (Dets[chan].type == 0 &&
                (Dets[chan].HGPostrecnt + Dets[chan].HGPrerecnt) < 2008) {
              presum[chan] = Dets[chan].HGPostrecnt + Dets[chan].HGPrerecnt;
              if (k > presum[chan]) {
                k = presum[chan];
                chan_k = chan;
              }
            }
            if (Dets[chan].type == 1 && Dets[chan].HGPreSumEnabled) {
              presum[chan] = 1 << Dets[chan].HGdecimationFactor;
              k = 1111;
              chan_k = chan;
            }
            printf("chan, type, presum = %d %d %d\n", chan, Dets[chan].type, presum[chan]);
          }
          if (chan >= 100 && chan < 100 + runInfo->nGe) {
            if (Dets[chan-100].type == 0 &&
                (Dets[chan-100].HGPostrecnt + Dets[chan-100].HGPrerecnt) < 2008) {
              presum[chan] = Dets[chan-100].HGPostrecnt + Dets[chan-100].HGPrerecnt;
              if (k > presum[chan]) {
                k = presum[chan];
                chan_k = chan;
              }
            }
            if (Dets[chan-100].type == 1 && Dets[chan-100].LGPreSumEnabled) {
              presum[chan] = 1 << Dets[chan-100].LGdecimationFactor;
              k = 1111;
              chan_k = chan;
            }
          }
        }
        if (k < 2008 && k > 1000) {  // yes, at least one channel on this card has presumming enabled
          printf("   Presumming: i=%2d; crate %d slot %2d has PS start at %d\n",
                 i, runInfo->GDcrate[i], runInfo->GDslot[i], k);
          if (k == 1111) { // GRETINA4A type, x16 expansion
            k = presum[chan_k]*2048;   // TEST TEMPORARY
          } else {
            k = 2048 + (2018-k)*3;  // max expanded signal length
          }
          for (j=0; j<NEVTS; j++) {
            if (!(sigBuf[nModBuf-1][j] = malloc(k * sizeof(short)))) {
              printf("Malloc failed for sigBuf[%d][%d]\n\n", nModBuf-1, j);
              return -1;
            }
          }
        } else if (DEBUG && chan_k >= 0) {
          printf("No presumming: i=%2d; crate %d slot %2d has k = %d\n",
                 i, runInfo->GDcrate[i], runInfo->GDslot[i], k);
          printf("HGPostrecnt, HGPrerecnt: %d %d\n",
                 Dets[chan_k].HGPostrecnt, Dets[chan_k].HGPrerecnt); fflush(stdout);
        }
      }

      for (i=0; i<runInfo->nV; i++) {
        if (VERBOSE)
          printf("%2d  veto card %2d of %2d crate %d slot %d\n",
                 nModBuf, i, runInfo->nV, runInfo->Vcrate[i], runInfo->Vslot[i]);
        /* malloc space to the corresponding module event buffer */
        if (!(modBuf[nModBuf++] = malloc(NEVTS * sizeof(BdEvent)))) {
          fprintf(stderr, "Malloc failed for Veto modBuf[%d][%d]\n\n", nModBuf-1, NEVTS);
          return -1;
        }
      }
      if (nModBuf > NBDS) {
        fprintf(stderr, "ERROR: Too many VME boards (%d)!\n\n", nModBuf);
        return -1;
      } else if (nModBuf != nmod) {
        fprintf(stderr, "ERROR: Number of VME boards (%d) != result of ep_init (%d)\n\n",
               nModBuf, nmod);
        return -1;
      }
    }
#ifdef PRESORT
    /* malloc space to the rundecoder event buffer */
    if (!(modBuf[nModBuf++] = malloc(NEVTS * sizeof(BdEvent)))) {
      printf("Malloc failed for Run modBuf[%d][%d]\n\n", nModBuf-1, NEVTS);
      return -1;
    }
    if (nModBuf > NBDS) {
      printf("ERROR: Too many VME boards (%d)!\n\n", nModBuf);
      return -1;
    }
#endif

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
          if (strstr(runInfo->decoder[i], "ORCV830DecoderForEvent"))      dataIdC830e = dataId[i];
          if (strstr(runInfo->decoder[i], "ORCV830DecoderForPolledRead")) dataIdC830p = dataId[i];
          if (strstr(runInfo->decoder[i], "ORCAEN792DecoderForQdc"))      dataIdC792  = dataId[i];
          if (strstr(runInfo->decoder[i], "ORCAEN792NDecoderForQdc"))     dataIdC792N = dataId[i];
        }
    }
    goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below
  }

  /* start loop over reading events from input file
     ============================================== */

  while (1) {
    oldCrate = crate;

    if (fread(head, sizeof(head), 1, f_in) != 1) break;
    recordID++;
    run_evts_count++;
    board_type = head[0] >> 18;
    evlen = (head[0] & 0x3ffff);
    ee = 0;

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
    if (j >= idNum) {  // the dataID type is not found in the file header
#ifndef QUIET
      printf("****  Unrecognized event type %d  %8.8x %8.8x  evlen = %d\n",
             board_type, head[0], head[1], evlen);
#endif
      n = evlen-2;
      if (n > 32) n = 32;
      fread(evtdat, 4, n, f_in);
      if (VERBOSE) {
        for (k=0; k<n; k++) {
          printf("   %8.8x %12d", evtdat[k], evtdat[k]);
          if (k%5 == 4) printf("\n");
        }
        c = (char *) evtdat;
        for (k=0; k<4*n; k++) {
          if (k%80 == 0) printf("\n");
          printf("%c", c[k]);
        }
        printf("\n\n");
      }
      evlen -= 2+n;
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

    if (DEBUG && totevts < 50)
      printf("Header: %8.8x %8.8x  Type: %3d Length:%4d  j = %2d -> %s\n",
             head[0], head[1], board_type, evlen, j, runInfo->decoder[j]);

    /* if we don't want to decode this type of data, just skip forward in the file */
    if (!goodDataId[j]) {
      fseek(f_in, 4*(evlen-2), SEEK_CUR);
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
        run_evts_count = 0;
        run_gretina_evts_count = 0;
        current_runNumber = evtdat[0];
        ntOutOfOrder = 0;

        
        if (n_good_run_list > 0) {
          good_run = 0;
          for (i=0; i < n_good_run_list; i++) {
            if (current_runNumber == good_run_list[i]) {
              good_run = 1;
              break;
            }
          }
        }
      } else if (head[1] & 0x8) {
        if (VERBOSE)
          printf(" -- %8.8x Heartbeat step %d at %d\n", head[1], evtdat[0], evtdat[1]);
      } else if ((head[1] & 0x29) == 0) {
        printf("------- END run %d at %d;  duration = %d s;  %d events -------\n",
               evtdat[0], evtdat[1], evtdat[1]-start_time, totevts);
        end_time = evtdat[1];
#ifdef PRESORT
        if (0) { // flush buffers at end of each run?
          printf("Flushing buffers...\n");
          if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                            &nBdsAvail, &out_evts, &built_evts, ps_f_out,
                            &badevts, his, lastBuiltTime) < 0) return -1;
          min_time = -1;
          min_time_board = -1;
        }
#endif
      } else {
        printf("***** ORRunDecoderForRun %8.8x %d %d\n",
               head[1], evtdat[0], evtdat[1]);
      }
      board = nModBuf-1;
      crate = 0;
      slot = 0;
      // time = min_time;
      time = last_gretina_time + 2000; // 20 us offset to avoid events being cut as part of pulser events
      /* ------------- end ORRunDecoderForRun end -------------- */

    } else {   // all other event types
      slot  = (head[1] >> 16) & 0x1f;
      crate = (head[1] >> 21) & 0xf;
      board = -1;

      if (crate < 0 || crate > NCRATES ||
          slot  < 0 || slot > 20 ||
          (board = module_lu[crate][slot]) < 0) {
        if (!DS0) {
          printf("ERROR: Illegal VME crate or slot number %d %2d; %8.8x %8.8x %d\n",
                 crate, slot, head[0], head[1], board);
          if (board >= 0) printf("  --> %s\n", runInfo->decoder[board]);
        }
        if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) break;
        continue;
      }
      if (VERBOSE && totevts < 50)
        printf("len, cr, slot, board: %d %d %d %d\n", evlen, crate, slot, board);
      his[3][1000 + board]++;

      /* ========== read in the rest of the eventdata ========== */
      if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
        printf("  No more data...\n");
        break;
      }
      if (!good_run) continue;
      totevts++;
#ifndef QUIET
      if (totevts % 100000 == 0)
        printf(" %8d evts in, %d out, %d built\n", totevts, out_evts, built_evts); fflush(stdout);
#endif
    }

    /* --------------- Gretina4M or Gretina4A digitizer ---------------- */
    if (board_type == dataIdGM || board_type == dataIdGA) {   // GRETINA digitizer
      // FIXME: Check for 0xaaaaaaaa, other integrity checks
      /* extract 48-bit time stamp */
      if (time > 0) oldTime  = time;
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

      his[3][200]++;
      /* if time is less than 0, then this is probably a bad event */
      if (time < 0) {
        his[3][201]++;
        continue;
      }

      ch = (evtdat[1] & 0xf);
      if (board >= 0 && ch < 10) {
        chan = chan_lu[board][ch];
#ifdef DETSEL
        if (!DETSEL) continue;
#endif
        signal = (short *) evtdat + 28;

        /* signal handling specific to Gretina4M data */
        if (board_type == dataIdGM) {
          if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
              2020 == decompress_signal((unsigned short *)signal, ucsig, 2*(evlen - 2) - 28)) {
            memcpy(signal, ucsig, 4040);
            evlen = 1026;
            head[0] = (head[0] & 0xfffc0000) + evlen;
          }
          siglen = 2018;
          erro = 0;
          if (chan >= 0 && chan < 200 && presum[chan] > 0) {   // need to correct for presumming
            // FIXME: replace hard-coded factor of 4 with data from header
            step = presum[chan];
            cts = (signal[step-5] + signal[step-4] + signal[step-3] + signal[step-2])/4;
            if (cts > 10 || cts < -10) {
              /* --- Find transition to 4x presumming ------- */
              for (step=presum[chan]; step<presum[chan]+2; step++) {
                if (cts > 0 && signal[step] > 2*cts) break;
                if (cts < 0 && signal[step] < 2*cts) break;
              }
              if (step == presum[chan]+2) {
                step = presum[chan]+1;
                erro = 1;
              } else if (DEBUG) {
                printf("step found at %d in chan %3d\n", step, chan);
              }
            }
            /* --- uncompress presumming --- */
            for (i=0; i<step; i++) ucsig[i+10] = signal[i];
            for (i=step; i<siglen; i++) {                     // NOTE: This is wrong for signal < 0
              for (j=0; j<4; j++)
                ucsig[step + 10 + 4*(i-step) + j] = signal[i]/4; // distribute sums over 4 bins each
              for (j=0; j<signal[i]%4; j++)
                ucsig[step + 10 + 4*(i-step) + j]++;             // and make sure sum is right
            }
            siglen += 3*(siglen-step);
            signal = &ucsig[10];
            /* --- Done with 4x presumming ------- */
            if (DEBUG)
              printf("Corrected chan %d for presumming at step %d, siglen = %d\n",
                     chan, step, siglen);
          }
        }

        /* signal handling specific to Gretina4A data */
        else if (board_type == dataIdGA) {
          prev_ch_time = evtdat[5];
          prev_ch_time = (prev_ch_time << 16) + (evtdat[4] >> 16);
          e_onbd = (((evtdat[9] & 0xffff) << 8) + (evtdat[8] >> 24) -
                    (evtdat[8] & 0xffffff));
          siglen = 2*evlen - 32;

          if ((k = presum[chan]) > 1) {  // presumming with factor k = presum[chan] // TEMPORARY 0
            istep1 = step = step2 = 0;
            istep2 = siglen;
            // printf("siglen, istep1, istep2 = %d %d %d\n", siglen, istep1, istep2); fflush(stdout);
            for (i=0; i<siglen; i++) {
              mkrs[i] = (unsigned short) signal[i] >> 14;
              if (0 && mkrs[i]) printf("mkr%d at %3d ; ", mkrs[i], i);
              if (step == 0 && mkrs[i] == 1) step = i;
              if (step2 == 0 && istep1 > 0 && mkrs[i] == 1) step2 = i;
              if (mkrs[i] > 1) {
                if (istep1 == 0) istep1 = i;
                else istep2 = i;
              }
              signal[i] = (signal[i] & 0x3fff) - 8192; // + (mkrs[i]/2) * 200; // TEMPORARY
            }
            // signal[0] = istep1; signal[1] = istep2;
            if (0 && chan < 57) {printf("chan, siglen, istep1, istep2 = %3d %4d %4d %4d\n",
                                        chan, siglen, istep1, istep2); fflush(stdout);}

            if (0 && istep1 > 0) { // correct for presumming
              for (i=j=0; i<istep1; i++, j+=k)
                for (n=0; n<k; n++) ucsig[j+n] = signal[i];
              for (; i<istep2; i++, j++)
                ucsig[j] = signal[i];
              i++;  // skip the first one of the presummed samples; it's redundant?
              for (; i<siglen; i++, j+=k)
              // for (; i<1200; i++, j+=k)
                for (n=0; n<k; n++) ucsig[j+n] = signal[i];
              siglen = j;
              signal = ucsig;
            }
          } else {                 // no presumming
            for (i=0; i<siglen; i++) signal[i] = (signal[i] & 0x3fff) - 8192;
          }
        }

        ee = trap_max(signal, &j, 401, 200);
#ifdef DETSEL
        // if (ee/401 < 4970 || ee/401 > 5700) continue;
#endif
        if (e_thresh[chan] > 0 && ee < e_thresh[chan] * 401) {
          subthreshevts++;
          continue;                         // discard as a bad event
        }
        if (board_type == dataIdGM && erro) {
          printf("Hmmm... step not found in %d - %d in chan %3d, counts %5d, ee = %.1lf\n",
                 step-1, step, chan, cts, (double) ee/401.0);
        }

        /* if necessary, save uncompressed signal in sigBuf */
        j = nevts[board] + iptr[board];
        if (j >= NEVTS) j -= NEVTS;
        if (siglen > 2018) {
          // printf(" j=%d\n", j); fflush(stdout);
          if (DEBUG) printf("uncompressed signal length = %d, board = %d, ptr = %d\n",
                            siglen, board, j); fflush(stdout);
          memcpy(sigBuf[board][j], ucsig, (siglen+10)*sizeof(short));
          signal = sigBuf[board][j]+10;
        } else {
          signal = (short *) modBuf[board][j].evbuf + 32;
        }
        modBuf[board][j].siglen = siglen;
        modBuf[board][j].sig = signal;
      }

      /* if we get here, then it's an okay event; time > 0 and energy > threshold */
      run_gretina_evts_count++;
      if (DEBUG && run_gretina_evts_count < 10) {
        printf("$> GRETINA evt %2d in run %d, chan %3d;  time %lld ms;  last_gretina_time %lld ms\n",
               run_gretina_evts_count, current_runNumber, chan, time/100000, last_gretina_time/100000);
      }

      /* if min_time < 0 or last_gretina_time < 0 then this is the first event with a valid time */
      if (min_time < 0) min_time = time;
      if (last_gretina_time < 0) last_gretina_time = time;

      if (run_gretina_evts_count < 2 &&            // we just started a new run
          last_gretina_time > time + 6000000000) { // and time *decreased* by dt > 1 minute
        /* this is a new run with reset time stamps
           flush build buffers */
        printf(" >>>>  Time stamps seem to have been reset! Flushing buffers!\n");
        if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
#ifdef PRESORT
                          &nBdsAvail, &out_evts, &built_evts, ps_f_out,
#else
                          &nBdsAvail, &out_evts, &built_evts,
#endif
                          &badevts, his, lastBuiltTime) < 0)
          return -1;
        min_time = oldTime = time;
        min_time_board = board;
        for (i=0; i<200; i++) last_ch_time[i] = 0;
        for (i=0; i<20; i++) last_bd_time[i] = 0;
      } else if (min_time > time+10) {  // small-ish time difference; has a board lost sync?
        if (out_evts > 0) {
          if (++ntOutOfOrder <= 10) {
            if ((min_time-time)/100000 < 10)
              printf("Hmmm, times are out of order... (%lld us < %lld us; d = %4lld us)"
                     " crate, slot = %d %2d (board %2d) min_time board = %d   %d evts into run %d\n",
                     time/100, min_time/100, (min_time-time)/100,
                     crate, slot, board, min_time_board,
                     run_evts_count, current_runNumber);
            else if ((min_time-time)/100000 < 10000)
              printf("Hmmm, times are out of order... (%lld ms < %lld ms; d = %5lld ms)"
                     " crate, slot = %d %2d (board %2d) min_time board = %d   %d evts into run %d\n",
                     time/100000, min_time/100000, (min_time-time)/100000,
                     crate, slot, board, min_time_board,
                     run_evts_count, current_runNumber);
            else 
              printf("Hmmm, times are out of order... (%lld  s < %lld  s; d = %6lld  s)"
                     " crate, slot = %d %2d (board %2d) min_time board = %d   %d evts into run %d\n",
                     time/100000000, min_time/100000000, (min_time-time)/100000000,
                     crate, slot, board, min_time_board,
                     run_evts_count, current_runNumber);
            // for (i=0; i<nModBuf; i++) if (nevts[i] > 0) printf("%4d", i);
            // printf("\n");
          }
          if (ntOutOfOrder == 10)
            printf(" Supressing further out-of-order messages. This run requires prebuilding!\n");
        }
        min_time = time;
        min_time_board = board;
      }
      last_gretina_time = time;
      dt = (time - min_time)/10000000;   // in 0.1 seconds
      if (dt > 3999) dt = 3999;
      his[3][dt+4000]++;
      if (start_time_Ge < 0) start_time_Ge = time / 100000000;

      /* add run number, recordID, and previous E and time, to Ge event */
      if (evtdat[8] != current_runNumber) {
        evtdat[8]  = current_runNumber;
        evtdat[9]  = recordID;
        evtdat[10] = last_ch_energy[chan];
        evtdat[11] = (time - last_ch_time[chan])/100; // us
      }
      if (last_bd_time[board] > time+800 && last_ch_time[chan] > time)
        printf(">> ERROR: board chan time previous_time = %2d %3d %12lld %12lld -> %14lld ms\n",
               board, chan, time/100000, last_bd_time[board]/100000,
               (time-last_bd_time[board])/100000);
      last_ch_energy[chan] = ee/401;
      last_ch_time[chan] = time;
      last_bd_time[board] = time;

    /* --------------- CV830DecoderForEvent --------------- */
    } else if (board_type == dataIdC830e) {
      /*
        word 0-1 = dataID, evlen, cate, slot  (as usual)      evlen = 7
        word 2 = 32-bit-timestamp roll-over count
        word 3 = enable mask? = 3
        word 4 = veto event count? (mask off top-most byte?)
        word 5 = 32-bit-timestamp
        word 6 = 0
       */
      veto_evts++;
      if (0) {
        printf("## %d  V830 event type %2d  %8.8x %8.8x  evlen = %2d",
               totevts, board_type, head[0], head[1], evlen);
        printf("  %8.8x %8.8x %8.8x %8.8x %8.8x\n",
               evtdat[0], evtdat[1], evtdat[2], evtdat[3], evtdat[4]);
      }
      if (!DS0 && veto_evts > 3 && sub_veto_count != 3) // counts how many boards are in a veto event;
        printf("ERROR: Bad number of boards (%d) in a veto event\n\n", sub_veto_count);
      sub_veto_count = 1;
      if (!DS0 && (evtdat[1] != 3 || evlen != 7)) {
        printf("ERROR: Unexpected enabled mask (0x%x)\n"
               "        or evlen (%d) in V830 veto data!\n\n", evtdat[1], evlen);
        continue;
      }
      time3 = time2;  // time of *last* veto event
      time2 = evtdat[0];
      time2 = (time2 << 32) | evtdat[3];  // time of this veto event
      dt = (time-time2)/1000000;  // Ge-Veto time in 0.01 s (10 ms)
      if (veto_evts > 2) {
        if (dt > -4000 && dt < 4000) his[0][4000+dt]++;
        dtsum += dt;
        dt /= 10;                   // Ge-Veto time in 0.1 s
        if (dt > -4000 && dt < 4000) his[1][4000+dt]++;
        dt /= 10;                   // Ge-Veto time in 1 s
        if (dt > -4000 && dt < 4000) his[2][4000+dt]++;
      }
      if (0)
        printf("## %d V830 time = %lld  Gtime = %lld    Diff = %lld  %lld\n",
               totevts, time2, time, (time2-time), (time2-time3));
      k = (evtdat[2] & 0xffffff);         // veto event count
      if (veto_count > 0 && k != veto_count+1 && veto_error_count++ < 40) {
        printf("ERROR: Bad V830 veto event count %d (0x%x), should be %d (0x%x)\n"
               " record %d  veto_evt %d  time %lld (0x%llx), last time %lld (0x%llx)  (diff = %lld ms)\n\n",
               k, k, veto_count+1, veto_count+1, totevts, veto_evts,
               time2, time2, time3, time3, (time2-time3)/100000);
        if (veto_error_count == 40)
          printf("Future veto error messages will be suppressed!\n");
      }
      veto_count = k;
      if (start_time_veto < 0) start_time_veto = time2 / 100000000;
      if (!BUILD_VETO) continue;

      /* ------------- CAEN792[N]DecoderForQdc -------------- */
    } else if (board_type == dataIdC792N) {
      /*
        word 0-1 = dataID, evlen, cate, slot  (as usual)     evlen = 21
        word 2 = ?
        word 3 = ?  (header should include number of converted chs)
        word 4-19 = data (bits 0-11), overflow(12), underthresh(13), channel(17-20)
        word 20 = EndOfBlock, includes event counter (bits 0-23)
        see v792.pdf manual, page 44     GEO = f8?
        head[0]  head[1]  evlen  evtdat[0] evtdat[1] evtdat[2] evtdat[3]
        002c0015 006d0001  21    5829c344  000612d4  f80060c2  f8106192
        002c0015 00720001  21    5829c344  00061385  f8006023  f810604d
       */
      if (0) {
        printf("## %d  V792N event type %d  %8.8x %8.8x  evlen = %d",
               totevts, board_type, head[0], head[1], evlen);
        printf("  %8.8x %8.8x %8.8x %8.8x %8.8x ... %8.8x %8.8x %8.8x\n",
               evtdat[0], evtdat[1], evtdat[2], evtdat[3], evtdat[4],
               evtdat[16], evtdat[17], evtdat[evlen-3]);
      } else if (0) {
        printf("## %d  V792N event type %d  %8.8x %8.8x  evlen = %d\n",
               totevts, board_type, head[0], head[1], evlen);
        for (i=0; i<evlen-2; i++) printf(" %8.8x", evtdat[i]);
        printf("\n");

      }
      if (evlen != 21 && evlen != 37 && !(DS0 && evlen == 19)) {
        printf("ERROR: Unexpected vlen (%d) in V792 veto data!\n\n", evlen);
        continue;
      }
      k = (evtdat[evlen-3] & 0xffffff);     // event count
      if (veto_count2 > 0 &&
          ((sub_veto_count == 1 && k != veto_count2+1) ||
           (sub_veto_count == 2 && k != veto_count2))) {
        if (veto_error_count++ < 40) {
          printf("ERROR: Bad V792 veto event count %d (0x%x), should be %d (0x%x)\n"
                 " record %d  veto_evt %d  sub_veto_count = %d\n\n",
                 k, k, veto_count2+2-sub_veto_count, veto_count2+2-sub_veto_count,
                 totevts, veto_evts, sub_veto_count);
          if (veto_error_count == 40)
            printf("Future veto error messages will be suppressed!\n");
        }
      }
      veto_count2 = k;
      sub_veto_count++;
      
      if (VERBOSE && veto_evts < 7) {
        printf("\n ch   UN  OV     ADC   veto_evt %d\n", veto_evts);
        for (i=2; i<18; i++)
          printf("%3d %4d %3d  %6d\n", (evtdat[i]>>17 & 0xf),
                 (evtdat[i]>>13 & 1), (evtdat[i]>>12 & 1), (evtdat[i] & 0xff));
      }
      if (!BUILD_VETO) continue;

      /* -------- these two seem to be unused data types --------- */
    } else if (board_type == dataIdC830p) { // unexpected!
      printf("****  V830 polled read; type %d  %8.8x %8.8x  evlen = %d\n",
             board_type, head[0], head[1], evlen);
      continue;
    } else if (board_type == dataIdC792) { // unexpected!
      printf("****  V792 event type %d  %8.8x %8.8x  evlen = %d\n",
             board_type, head[0], head[1], evlen);
      continue;

      /* --------------- end of veto info ---------------- */

    } else if (board_type != dataIdRun) { 
      continue;   // FIXME; what else needs to be implemented??
    }
    
    if (VERBOSE && totevts < 50) printf("    time: %lld\n", time);

    /* ------------------------------------------------------------------ */
    /* store the event in the event-building buffer for this data channel */
    j = nevts[board] + iptr[board];
    if (j >= NEVTS) j -= NEVTS;

    /* first check that the timestamp is in order with earlier ones */
    if (0) {    // this doesn't seem to actually make any significant difference
      for (n = 0; n <= nevts[board]; n++) {
        k = j;
        if (--k < 0) k += NEVTS;
        if (modBuf[board][k].time < time) break;  // if order is ok, break
        // else shift already-stored data to later in the module buffer
        memcpy(&modBuf[board][j], &modBuf[board][k], sizeof(BdEvent));
        if (modBuf[board][k].siglen <= 2018) {  // signal was not presummed
          modBuf[board][j].sig = (short *) modBuf[board][j].evbuf + 32;
        } else {                                // signal was presummed
          memcpy(&sigBuf[board][j], &sigBuf[board][k],
                 (modBuf[board][j].siglen+10)*sizeof(short));
          modBuf[board][j].sig = sigBuf[board][j] + 10;
        }
      }
    }

    modBuf[board][j].evtID = totevts;
    modBuf[board][j].evlen = evlen;
    modBuf[board][j].crate = crate;
    modBuf[board][j].slot  = slot;
    modBuf[board][j].orca_type = board_type;
    modBuf[board][j].time = time;
    if (board_type == dataIdC830e ||
        board_type == dataIdC792N) modBuf[board][j].time = time2;
    modBuf[board][j].evbuf[0] = head[0];
    modBuf[board][j].evbuf[1] = head[1];
    memcpy(modBuf[board][j].evbuf+2, evtdat, (evlen-2)*4);

    if (crate > 0 && crate < 3 && oldCrate > 0 && oldCrate < 3) {
      if (oldTime == 0) oldTime = start_time;
      dt = (time - oldTime)/100000000; // time difference in seconds
      if (dt<-3999) dt = -3999;
      if (dt>3999)  dt = 3999;
      if (crate != oldCrate) {
        his[4+crate][dt+4000]++;  // crate-to-crate time differences (y=5,6)
        if (dt<-49) dt = -49;
        if (dt>49) dt = 49;
        his[9][dt + 4000*(crate-1) + 100*slot]++;    // slot-to-slot time differences (y=9)
      } else {
        his[6+crate][dt+4000]++;  // time differences within one crate (y=7,8)
        if (dt<-49) dt = -49;
        if (dt>49) dt = 49;
        his[10][dt + 4000*(crate-1) + 100*slot]++;   // slot-to-slot time differences (y=10)
      }
    }

    if (++nevts[board] == 1) nBdsAvail++;
    if (nevts[board] == NEVTS) {
      printf("Ackkk! Board %d has %d events stored! Flushing buffers...\n",
             board, nevts[board]);
      if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
#ifdef PRESORT
                        &nBdsAvail, &out_evts, &built_evts, ps_f_out,
#else
                        &nBdsAvail, &out_evts, &built_evts,
#endif
                        &badevts, his, lastBuiltTime) < 0)
        return -1;
      min_time = time;
      min_time_board = board;
    }
    if (crate < 3) his[3][205+crate]++;
    if (board_type == dataIdRun) continue;

    /* histogram buffer time difference and nevts */
    dt = (time - modBuf[board][iptr[board]].time)/20000000; // 0.2 s
    if (dt < 0 || dt > 99) dt = 99;
    his[4][dt + 2000 + 100*board]++;
    his[4][nevts[board]*100/NEVTS + 4000 + 100*board]++;
    if (nevts[board] > NEVTS-2) {  // the board buffer is full; histogram max time diff
      his[4][dt + 100*board]++;
    }

    /* nBdsAvail is the number of data channels for which
       we have at least one stored event */
    while (//ntOutOfOrder == 3 &&   // essentially the same as "DoNotBuild"?
           nBdsAvail > 0 &&
           totevts >= NEVTS &&
           (nevts[board] > NEVTS-2 ||
            modBuf[board][j].time - min_time > BUFFER_TIME)) {
      /* time to output some events; first sort the available events in time */
      n = 0;
      for (i=0; i<nModBuf; i++) {
	if (nevts[i] <= 0) continue;
	if (n == 0 ||
	    modBuf[order[n-1]][iptr[order[n-1]]].time <= modBuf[i][iptr[i]].time) {
	  order[n++] = i;
	} else {
	  for (j=0; j<n; j++) {
	    if (modBuf[order[j]][iptr[order[j]]].time > modBuf[i][iptr[i]].time) {
	      for (k=n; k>j; k--) order[k] = order[k-1];
	      order[j] = i;
	      n++;
	      break;
	    }
	  }
	}
      }
      if (n != nBdsAvail || n == 0) {
        printf("\nAckk! In eventbuild, n = %d != nBdsAvail = %d !\n\n",
               n, nBdsAvail);
        return -1;
      }
      min_time = modBuf[order[0]][iptr[order[0]]].time;  // earliest time in all module buffers
      min_time_board = order[0];

      if (VERBOSE && out_evts < 500)
        printf("%4d totevts; order[0] = %d\n", totevts, order[0]);

      /* now collect and assemble the built-event */
      nChData = 0;
      for (i=0; i<n; i++) {   // assemble and process board data
	j = order[i];
	while (//modBuf[j][iptr[j]].time > min_time - 20 &&   // FIXME
	       modBuf[j][iptr[j]].time - min_time < RESOLVING_TIME) {
          // extend event time by later times? if so, uncomment next line...
          // if (min_time < modBuf[j][iptr[j]].time) min_time = modBuf[j][iptr[j]].time;
          if (modBuf[j][iptr[j]].time - min_time < -1000000000) { // 10 s; clocks must have been reset?
            if (modBuf[j][iptr[j]].orca_type != dataIdRun)
              printf("Error: Large jump back in timestamps. Clocks reset?"
                     " chan, min_time, time: %d %lld %lld ms\n",
                     modBuf[j][iptr[j]].chan, min_time/100000, modBuf[j][iptr[j]].time/100000);
            break;
          }
          his[3][210+modBuf[j][iptr[j]].crate]++;
          if (nChData >= 200) {
            printf("ERROR: More than 200 sub-events in built event!\n\n");
            printf("   min_time = %12lld\n\n Ordering:\n", min_time);
            for (k=0; k<n; k++)
              printf("%12lld\n", modBuf[order[k]][iptr[order[k]]].time);
            printf("\n");
            for (k=0; k<200; k+=4)
              printf("%12lld %12lld %12lld %12lld\n", ChData[k]->time,
                     ChData[k+1]->time, ChData[k+2]->time, ChData[k+3]->time);
            //return -1;
            break;
          }
          out_evts++;
          if (VERBOSE && out_evts < 500)
            printf("%8d out; i,j,time,mintime: %2d %2d %lld %lld\n",
                   out_evts, i, j, modBuf[j][iptr[j]].time, min_time);
          ChData[nChData++] = &modBuf[j][iptr[j]];  // insert channel-event into built-event
          if (++iptr[j] >= NEVTS) iptr[j] -= NEVTS; // and remove from module buffer
          if (--nevts[j] == 0) {
            nBdsAvail--;
            break;
          }
        }
      }
      /* now process the built-event, and write to skim file for PRESORT */
      /* first histogram time since last built event */
      dt = (ChData[0]->time - lastBuiltTime) / 100;  // us
      lastBuiltTime = ChData[0]->time;
      if (dt > -1000 && dt < 7000) his[11][1000+dt]++;
      dt /= 1000; // ms
      if (dt > -1000 && dt < 7000) his[12][1000+dt]++;

      /* eventprocess return values:
            -1 on error,
             1 on pulser event,
             2 on bad event (bad = cannot be processed; not dirty)
             0 otherwise (good event, can be clean or dirty)
      */
      built_evts++;
      if ((k = eventprocess(Dets, runInfo, nChData, ChData)) < 0) return k;
#ifdef PRESORT
      if (k < 1) {  // good event, not a pulser, but before data cleaning; i.e. can be dirty
        for (i = 0; i < nChData; i++) {
          if (!DO_NOT_COMPRESS &&
              (ChData[i]->orca_type == dataIdGM ||
               ChData[i]->orca_type == dataIdGA) &&   // GRETINA digitizer
              ChData[i]->evlen == 1026) {  // compress the signal
            signal = (short *) ChData[i]->evbuf + 32;
            int m = compress_signal(signal, (unsigned short *) ucsig, 2020);
            memcpy(signal, ucsig, m*2);
            ChData[i]->evlen = 16 + m/2;
            ChData[i]->evbuf[0] = (ChData[i]->evbuf[0] & 0xfffc0000) + ChData[i]->evlen;
          }
          fwrite(ChData[i]->evbuf, sizeof(int), ChData[i]->evlen, ps_f_out);
        }
      }
#else
      if (k > 1) badevts++;
#endif
    }

    // if (totevts > 200) break;
    // FIXME - need to check for timestamp rollover?

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  end_time_veto = time2/100000000;
  end_time_Ge = time/100000000;

  if (VERBOSE)
    printf("FINISHED reading file\n");

  /* flush remaining data in all buffers (essentially a copy of code from above) */
#ifndef PRESORT
  if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                    &nBdsAvail, &out_evts, &built_evts,
                    &badevts, his, lastBuiltTime) < 0)
    return -1;
#endif

  printf("\n %d events in, %d events out, %d events built, %d * 3 veto events\n",
         totevts, out_evts, built_evts, veto_evts);
#ifndef QUIET
  printf("    Ge delta-time = %d s   Veto delta-time = %d s\n"
         "    Mean Ge-Veto time = %.0lf ms\n",
         end_time_Ge - start_time_Ge, end_time_veto - start_time_veto,
         10.0*((double) dtsum)/((double) (veto_evts-2)));
  if (end_time_Ge - start_time_Ge > end_time - start_time + 2)
    printf(" ****  Ge delta-time is %d s greater than run duration!\n",
           end_time_Ge - start_time_Ge - end_time + start_time);
  if (end_time_veto - start_time_veto > end_time - start_time + 2)
    printf(" ****  Veto delta-time is %d s greater than run duration!\n",
           end_time_veto - start_time_veto - end_time + start_time);
#endif
  printf(" %d bad events, %d sub-threshold events\n", badevts, subthreshevts);

  /* if needed, write out diagnostic histograms */
  if (1) {
    f_out = fopen("eb.rms", "w");
    for (i=0; i<16; i++) {
      sprintf(spname, "%d EB run %d", i, runInfo->runNumber);
      if (i==0)  strcat(spname, "; Ge-Veto dt [10 ms]");
      if (i==1)  strcat(spname, "; Ge-Veto dt [0.1 s]");
      if (i==2)  strcat(spname, "; Ge-Veto dt [s]");
      if (i==3)  strcat(spname, "; evt/crate/bd/delta-t stats");
      if (i==4)  strcat(spname, "; buff dt/nevts stats");
      if (i==5)  strcat(spname, "; buff crate-dt [s], crate 1-2");
      if (i==6)  strcat(spname, "; buff crate-dt [s], crate 2-1");
      if (i==7)  strcat(spname, "; buff crate-dt [s], crate 1-1");
      if (i==8)  strcat(spname, "; buff crate-dt [s], crate 2-2");
      if (i==9)  strcat(spname, "; buff slot-dt [s], diff crates");
      if (i==10) strcat(spname, "; buff slot-dt [s], same crate");
      if (i==11) strcat(spname, "; time between built events, us");
      if (i==12) strcat(spname, "; time between built events, ms");
      // if (i==13) strcat(spname, "; time between built events, s");
      write_his(his[i], 8192, i, spname, f_out);
    }
    fclose(f_out);
  }

#ifndef PRESORT
  /* special call to eventprocess() with nChData = -99 as flag
     to finish up processing, write out files, etc. */
  if ((k = eventprocess(Dets, runInfo, -99, ChData)) < 0) return k;
#endif

#ifdef PRESORT
  for (i=0; i<NBDS; i++)
    if (t_offset[i]) printf("prebuild t_offset[%d] = %lld\n", i, t_offset[i]);
#endif

  if (ntOutOfOrder > 5)
    printf("\n@@@@@ Time stamps were out of order. This run %d requires prebuilding!\n",
           runInfo->runNumber);

  return  out_evts;
} /* eventbuild() */

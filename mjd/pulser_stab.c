#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"
#include "runBits.h"

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

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define VERBOSE 0
#define DEBUG   0
#define NEVTS  50000               /* no. of events to use in each board's event buffer
                                      (should only need ~ 50 in low-rate situations) */
#define RESOLVING_TIME        400  /* up to 4 us delay allowed within one built event */ 
#define BUFFER_TIME (20*100000000) /* up to 20s delay allowed before forcing event building */
/*  ------------------------------------------------------------ */

int eventprescan(FILE *f_in, MJDetInfo *detInfo, MJRunInfo *runInfo);
int read_his(int *his, int idimsp, int spid, char *namesp, FILE *file);
int run = 0;


int main(int argc, char **argv) {

  FILE       *f_in, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  MJRunInfo  runInfo;
  BdEvent    **ChData = NULL;
  int        nDets, i, argn=1, nEnab=0;
  char       *c, fname[256], line[256];
  int        start = 0, nn = 0;

  
  if (argc < 2) {
    fprintf(stderr,
            "\nusage: %s DS_fname_in"
            "\n   or: %s <list_fname>.lis\n\n", argv[0], argv[0]);
    return -1;
  }

  /* get data file name from input list (either .lis file or command arguments) */
  while (argn < argc && argv[argn][0] == '-') argn += 2;
  if (strstr(argv[argn], ".lis")) {    // input argument is a list file
    if ((f_lis = fopen(argv[argn],"r")) == NULL) {
      fprintf(stderr, "\n Failed to open list input file %s\n", argv[argn]);
      return 0;
    }
    if (!fgets(fname, sizeof(fname), f_lis)) {
      fprintf(stderr, "\n Failed to read list input file %s\n", argv[argn]);
      return 0;
    }
    if (argc > 2) {
      start = atoi(argv[2]);  // continue from a previous set of runs; skip over <start> files
      if (start < 0) start = 0;
      if (start > 0) {
        run = start;
        for (i=0; i<start; i++) {
          if (!fgets(fname, sizeof(fname), f_lis)) {
            fprintf(stderr, "\n Failed to read list input file %s\n", argv[argn]);
            return 0;
          }
        }
      }
      if (argc > 3) {
        nn = atoi(argv[3]);
        if (nn < 0) nn = 0;
      }
    }
    printf(" >>>>>>>> note <<<<<<<<< run = %d  start = %d  nn = %d\n", run, start, nn);
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk

  } else {   // using command argument list for input files
    strncpy(fname, argv[argn], sizeof(fname));
  }
  /* open raw data file as input */
  if ((f_in = fopen(fname,"r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", fname);
    return 0;
    }
    printf("\n >>> Reading %s\n\n", fname);
    strncpy(runInfo.filename, fname, sizeof(runInfo.filename));
    runInfo.argc = argc;
  runInfo.argv = argv;

  /* get all the required information for the file header */
  nDets = decode_runfile_header(f_in, detInfo, &runInfo);
  if (nDets < 1) return 1;

  if (!runInfo.flashcam) {
    printf("\n# PT_ID  Dig\n");
    for (i=0; i<runInfo.nPT; i++) {
      printf(" %3d  %d,%2.2d,%d\n", i,
             runInfo.PTcrate[i], runInfo.PTslot[i], runInfo.PTchan[i]);
    }

    if (runInfo.dataIdGM == 0 && runInfo.dataIdGA == 0) {
      printf("\n No data ID found for Gretina4M or 4A data!\n");
      return 1;
    }

    if (runInfo.dataIdGM)
      printf("\n Data ID %d found for Gretina4M data\n", runInfo.dataIdGM);
    if (runInfo.dataIdGA)
      printf("\n Data ID %d found for Gretina4A data\n", runInfo.dataIdGA);
    printf(" Run number: %d in file %s\n"
           " Start time: %s  (%d)\n"
           " Run bits  : %d = 0x%8.8x\n",
           runInfo.runNumber, runInfo.filename, runInfo.date, runInfo.startTime,
           runInfo.runType, runInfo.runType);
    for (i=0; i<32; i++) {
      if (runInfo.runType & 1<<i) printf("  0x%8.8x - %s\n", 1<<i, runBitDesc[i]);
    }
  }      // if (!runInfo.flashcam) {

  for (i=0; i<nDets; i++)
    if (detInfo[i].HGChEnabled || detInfo[i].LGChEnabled) nEnab++;
  printf("\n  %d of %d detectors are enabled.\n", nEnab, nDets);

  runInfo.analysisPass = 0;
  fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go to start of data in input file

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file, to build and process them */
    i = eventprescan(f_in, detInfo, &runInfo);
    fclose(f_in);
    if (i == -1) {
      fprintf(stderr, "\n ERROR: Scan quit with error\n\n");
      return i;
    }
    if (nn > 0) {
      printf(" >>>>>>>> note <<<<<<<<< run = %d  start = %d  nn = %d\n", run, start, nn);
      if (run >= start+nn) break;
    }

    /* get next input file name */
    if (f_lis) {  // using list file for input names
      if (!fgets(fname, sizeof(fname), f_lis))  // no more lines in the input list file
        strncpy(fname, "0", strlen(fname));     // special string
      for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
    } else {  // using command argument list for input names
      if (argn == argc-1) {                     // just read last input file
        strncpy(fname, "0", strlen(fname));     // special string
      } else {
        strncpy(fname, argv[++argn], sizeof(fname));
      }
    }
    if (strlen(fname) < 2) {  // no more files to process
      runInfo.analysisPass = -1;
      //       special flag  ^^ to tell event builder to run end-of-run cleanup
      i = eventprescan(f_in, detInfo, &runInfo);
      if (i == -1) {  // error
        fprintf(stderr, "\n ERROR: Scan quit with error\n\n");
        return i;
      }
      break;
    }

    /* open next file and read header size and run number */
    if ((f_in = fopen(fname,"r")) == NULL) {
      fprintf(stderr, "\n Failed to open input file %s\n", fname);
      return 0;
    }
    strncpy(runInfo.filename, fname, sizeof(runInfo.filename));
    printf("\n >>> Reading %s\n\n", fname);
    fread(&runInfo.fileHeaderLen, 1, sizeof(int), f_in);
    /* loop through the lines of the XML data until we find the run number */
    while (fgets(line, sizeof(line), f_in) && strncmp(line, "</plist>", 8)) {
      if (strstr(line, "<key>RunNumber</key>")) {
        fgets(line, sizeof(line), f_in);
        if (!(c=strstr(line, "<integer>")) ||
            (1 != sscanf(c+9, "%d", &runInfo.runNumber))) {
          fprintf(stderr, "\n ERROR decoding run number:\n %s\n", line);
          return -1;
        }
        break;
      }
    }
    /* position at start of data */
    fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);

    /* see if we need to read new pulser-tag info */
    runInfo.analysisPass = -17;
    eventprocess(detInfo, &runInfo, 0, ChData);
    runInfo.analysisPass = 0;   
  }

  printf("\n All Done.\n\n");
  return i;
}

int flush_buffers(MJDetInfo *Dets, MJRunInfo *runInfo, BdEvent *modBuf[NBDS], int nModBuf,
                  int nevts[NBDS], int iptr[NBDS], int order[NBDS], int *nBdsAvail,
                  int *out_evts, int *built_evts, int *badevts, int his[400][16384]) {

  int     i, j, k, n, nChData;
  BdEvent *ChData[200];        // pointers for constructing a built-event
  long long int min_time;


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
      while (modBuf[j][iptr[j]].time - min_time < RESOLVING_TIME) {
        // extend event time by later times? if so, uncomment next line...
        // if (min_time < modBuf[j][iptr[j]].time) min_time = modBuf[j][iptr[j]].time;
        if (modBuf[j][iptr[j]].time - min_time < -1000000000) { // 10 s; clocks must have been reset?
          printf("Error: Large jump back in timestamps. Clocks reset?\n");
          break;
        }
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
    /* now analyze/process the built-event */
    /* eventprocess return values:
       -1 on error,
       1 on pulser event,
       2 on bad event (bad = cannot be processed; not dirty)
       0 otherwise (good event, can be clean or dirty)
    */
    (*built_evts)++;
    if ((k = eventprocess(Dets, runInfo, nChData, ChData)) < 0) return k;

    if (k == 1) { //pulser event
      his[98][nChData]++;
      for (i=0; i<nChData; i++) {
        his[98][100+ChData[i]->det]++;
        his[98][200+ChData[i]->chan]++;
        if ((j = ChData[i]->chan) < 60) {
          short *signal = ChData[i]->sig;
          int e = ChData[i]->e;
          if (e > 100 && e < 3000) {
            his[j][4000+run]++;
            his[100+j][e]++;
            his[100+j][4000+run] += e;
            int s = 50, tmax;
            for (int m=0; m<100; m++) s += signal[m+1150] - signal[m+1850];                // related to tau
            s /= 100;
            his[200+j][2000 + s]++;
            his[200+j][4000+run] += s;
            s = (int) (0.5 + 275.0 * trap_max_range(signal, &tmax, 8, 0, 800, 1200) / e);  // A/E
            if (s > 0 && s < 2000) {
              his[300+j][1000 + s]++;
              his[300+j][4000+run] += s;
            }
          }
        }
      }
    }

    if (k > 1) (*badevts)++;  // CHECKME
  }
  printf("Flushing finished... %d\n", run+1); // fflush(stdout);
  for (j=0; j<60; j++) {
    if (his[j][4000+run] > 5) {
      double a1, a2, a3;
      a1 = 10.0    * his[100+j][4000+run] / his[j][4000+run];      // height
      a2 = 10000.0 * his[200+j][4000+run] / his[100+j][4000+run];  // tail slope
      a3 = 10.0    * his[300+j][4000+run] / his[j][4000+run];      // A/E
      if (his[100+j][0] == 0) {
        his[100+j][0] = lrint(a1);                                 // first non-zero value, for normalization
        his[200+j][0] = lrint(a2);
        his[300+j][0] = lrint(a3);
      }
      his[100+j][4000+run] = lrint(1000.0 * a1 / his[100+j][0]);   // normalize starting height to 1000
      his[200+j][4000+run] = lrint(1000.0 * a2 / his[200+j][0]);
      his[300+j][4000+run] = lrint(1000.0 * a3 / his[300+j][0]);
    } else {
      his[100+j][4000+run] = his[200+j][4000+run] = his[300+j][4000+run] = 0;
    }
  }
  his[99][4000+run] = run+1;

  return 0;
}

/*  ------------------------------------------------------------ */


int eventprescan(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int    i, j, k, n, crate=0, slot=0, board, evlen, ch, chan=0, nmod, ee;
  unsigned int  head[2], evtdat[20000];
  static int    dataId[32], goodDataId[32], idNum, board_type;
  static int    dataIdGM=0, dataIdGA=0, dataIdRun=0;
  static int    totevts, badevts, subthreshevts, out_evts, built_evts, recordID;

  long long int time=0, min_time = -1, prev_ch_time=0;
  int           min_time_board=-1, istep1, step, step2, siglen = 2018;
  static int    ntOutOfOrder = 0, oldCrate=0;
  static int    start_time;
  static long long int dt, dtsum, oldTime=0;

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
  static int his[400][16384] = {{0}};  // histograms for mean energy, tau, and A/E for each run

  static int first = 1;
  static int first_runNumber = 0, current_runNumber = 0, run_evts_count = 0, run_gretina_evts_count = 0;
  char     *c, spname[64], line[256];
  short    *signal, mkrs[8192];
  FILE     *f_in2, *f_out;
  static int           last_ch_energy[200] = {0};
  static long long int last_ch_time[200] = {0}, last_bd_time[20] = {0};
  static long long int t_offset[200] = {0}, last_gretina_time = -1;
  int  cts=0, erro=0, chan_k;


  if (runInfo->analysisPass < 0) {
    /* end-of-run cleanup */
    if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                      &nBdsAvail, &out_evts, &built_evts, &badevts, his) < 0)
      return -1;

    /* change filename stored in runInfo for the benefit of ep_finalize and plot titles, etc */
    printf(" Flushing eventprocess...\n");
    sprintf(runInfo->filename, "DS%d", first_runNumber);
    if ((k = eventprocess(Dets, runInfo, -99, ChData)) < 0) return k;
    printf("\n %d events in, %d events out, %d events builts\n",
           totevts, out_evts, built_evts);
    printf(" %d bad events, %d sub-threshold events\n", badevts, subthreshevts);
    return  out_evts;
  }
  
  /* initialize */
  current_runNumber = runInfo->runNumber;

  /* initialize some values */
  if ((nmod = ep_init(Dets, runInfo, module_lu, det_lu, chan_lu)) < 0) return -1;
  if (nModBuf > 0 && nmod != nModBuf) {
    fprintf(stderr, "ERROR: Number of VME boards has changed from %d to %d!\n\n", nModBuf, nmod);
    // return -1;
  }

  if (first) {
    first = 0;
    i = runInfo->firstRunNumber = first_runNumber = current_runNumber = runInfo->runNumber;
    totevts = badevts = subthreshevts = out_evts = built_evts = dtsum = recordID = 0;
    min_time = -1;
    start_time = 0;

    if (run > 0) {  // continuing from a previous set of runs; read pulstab.rms
      f_out = fopen("pulstab.rms", "r");
      if (!(f_out = fopen("pulstab.rms", "r"))) {
        fprintf(stderr, "\n Failed to open pulstab.rms as input file\n");
        exit(-1);
      }
      for (i=0; i<400; i++) {
        read_his(his[i], 16384, i, spname, f_out);
      }
      fclose(f_out);
    }

    for (i=0; i<NBDS; i++) {
      nevts[i] = iptr[i] = 0;
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
                (Dets[chan-100].LGPostrecnt + Dets[chan-100].LGPrerecnt) < 2008) {
              presum[chan] = Dets[chan-100].LGPostrecnt + Dets[chan-100].LGPrerecnt;
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
    /* malloc space to the rundecoder event buffer */
    if (!(modBuf[nModBuf++] = malloc(NEVTS * sizeof(BdEvent)))) {
      printf("Malloc failed for Run modBuf[%d][%d]\n\n", nModBuf-1, NEVTS);
      return -1;
    }
    if (nModBuf > NBDS) {
      printf("ERROR: Too many VME boards (%d)!\n\n", nModBuf);
      return -1;
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
        runInfo->runNumber = runInfo->firstRunNumber = current_runNumber = evtdat[0];
        ntOutOfOrder = 0;
      } else if (head[1] & 0x8) {
        if (VERBOSE)
          printf(" -- %8.8x Heartbeat step %d at %d\n", head[1], evtdat[0], evtdat[1]);
      } else if ((head[1] & 0x29) == 0) {
        printf("------- END run %d at %d;  duration = %d s;  %d events -------\n",
               evtdat[0], evtdat[1], evtdat[1]-start_time, totevts);

        // flush buffers at end of each run
        his[98][4000 + run] = evtdat[0];
        printf("Flushing buffers...\n");
        if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                          &nBdsAvail, &out_evts, &built_evts, &badevts, his) < 0) return -1;
        run++;
        min_time = -1;
        min_time_board = -1;
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

      /* ========== read in the rest of the eventdata ========== */
      if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
        printf("  No more data...\n");
        break;
      }
      totevts++;
      if (totevts % 100000 == 0)
        printf(" %8d evts in, %d out, %d built\n", totevts, out_evts, built_evts); fflush(stdout);
    }

    /* --------------- Gretina4M or Gretina4A digitizer ---------------- */
    if (board_type == dataIdGM || board_type == dataIdGA) {   // GRETINA digitizer
      // here we could also add a check for 0xaaaaaaaa, other integrity checks
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
      /* if time is less than 0, then this is probably a bad event */
      if (time < 0) continue;

      ch = (evtdat[1] & 0xf);
      if (board >= 0 && ch < 10) {
        chan = chan_lu[board][ch];
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
          siglen = 2*evlen - 32;

          if ((k = presum[chan]) > 1) {  // presumming with factor k = presum[chan] // TEMPORARY 0
            istep1 = step = step2 = 0;
            for (i=0; i<siglen; i++) {
              mkrs[i] = (unsigned short) signal[i] >> 14;
              if (step == 0 && mkrs[i] == 1) step = i;
              if (step2 == 0 && istep1 > 0 && mkrs[i] == 1) step2 = i;
              if (mkrs[i] > 1) {
                if (istep1 == 0) istep1 = i;
              }
              signal[i] = (signal[i] & 0x3fff) - 8192; // + (mkrs[i]/2) * 200; // TEMPORARY
            }
          } else {                 // no presumming
            for (i=0; i<siglen; i++) signal[i] = (signal[i] & 0x3fff) - 8192;
          }
        }

        ee = trap_max(signal, &j, 401, 200);
        if (e_thresh[chan] > 0 && ee < e_thresh[chan] * 401) {
          subthreshevts++;
          continue;                // discard as a bad event
        }
        if (board_type == dataIdGM && erro) {
          printf("Hmmm... step not found in %d - %d in chan %3d, counts %5d, ee = %.1lf\n",
                 step-1, step, chan, cts, (double) ee/401.0);
        }

        /* if necessary, save uncompressed signal in sigBuf */
        j = nevts[board] + iptr[board];
        if (j >= NEVTS) j -= NEVTS;
        if (siglen > 2018) {
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
                          &nBdsAvail, &out_evts, &built_evts, &badevts, his) < 0)
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
          }
          if (ntOutOfOrder == 10)
            printf(" Supressing further out-of-order messages. This run requires prebuilding!\n");
        }
        min_time = time;
        min_time_board = board;
      }
      last_gretina_time = time;

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

    } else if (board_type != dataIdRun) { 
      continue;   // FIXME; what else needs to be implemented??
    }
    
    if (VERBOSE && totevts < 50) printf("    time: %lld\n", time);

    /* ------------------------------------------------------------------ */
    /* store the event in the event-building buffer for this data channel */
    j = nevts[board] + iptr[board];
    if (j >= NEVTS) j -= NEVTS;

    modBuf[board][j].evtID = totevts;
    modBuf[board][j].evlen = evlen;
    modBuf[board][j].crate = crate;
    modBuf[board][j].slot  = slot;
    modBuf[board][j].orca_type = board_type;
    modBuf[board][j].time = time;
    modBuf[board][j].evbuf[0] = head[0];
    modBuf[board][j].evbuf[1] = head[1];
    memcpy(modBuf[board][j].evbuf+2, evtdat, (evlen-2)*4);

    if (crate > 0 && crate < 3 && oldCrate > 0 && oldCrate < 3) {
      if (oldTime == 0) oldTime = start_time;
    }

    if (++nevts[board] == 1) nBdsAvail++;
    if (nevts[board] == NEVTS) {
      printf("Ackkk! Board %d has %d events stored! Flushing buffers...\n",
             board, nevts[board]);
      if (flush_buffers(Dets, runInfo, modBuf, nModBuf, nevts, iptr, order,
                        &nBdsAvail, &out_evts, &built_evts, &badevts, his) < 0)
        return -1;
      min_time = time;
      min_time_board = board;
    }
    if (board_type == dataIdRun) continue;

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
	while ( modBuf[j][iptr[j]].time - min_time < RESOLVING_TIME) {
          // extend event time by later times? if so, uncomment next line...
          // if (min_time < modBuf[j][iptr[j]].time) min_time = modBuf[j][iptr[j]].time;
          if (modBuf[j][iptr[j]].time - min_time < -1000000000) { // 10 s; clocks must have been reset?
            if (modBuf[j][iptr[j]].orca_type != dataIdRun)
              printf("Error: Large jump back in timestamps. Clocks reset?"
                     " chan, min_time, time: %d %lld %lld ms\n",
                     modBuf[j][iptr[j]].chan, min_time/100000, modBuf[j][iptr[j]].time/100000);
            break;
          }
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
      /* now process the built-event */
      /* eventprocess return values:
            -1 on error,
             1 on pulser event,
             2 on bad event (bad = cannot be processed; not dirty)
             0 otherwise (good event, can be clean or dirty)
      */
      built_evts++;
      if ((k = eventprocess(Dets, runInfo, nChData, ChData)) < 0) return k;

      if (k == 1) { //pulser event
        his[98][nChData]++;
        for (i=0; i<nChData; i++) {
          his[98][100+ChData[i]->det]++;
          his[98][200+ChData[i]->chan]++;
          if ((j = ChData[i]->chan) < 60) {
            short *signal = ChData[i]->sig;
            int e = ChData[i]->e;
            if (e > 100 && e < 3000) {
              his[j][4000+run]++;
              his[100+j][e]++;
              his[100+j][4000+run] += e;
              int s = 50, tmax;
              for (int m=0; m<100; m++) s += signal[m+1150] - signal[m+1850];                // related to tau
              s /= 100;
              his[200+j][2000 + s]++;
              his[200+j][4000+run] += s;
              s = (int) (0.5 + 275.0 * trap_max_range(signal, &tmax, 8, 0, 800, 1200) / e);  // A/E
              if (s > 0 && s < 2000) {
                his[300+j][1000 + s]++;
                his[300+j][4000+run] += s;
              }
            }
          }
        }
      }

      if (k > 1) badevts++;
    }

    // if (totevts > 200) break;
    // FIXME - need to check for timestamp rollover?

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  if (VERBOSE)
    printf("FINISHED reading file\n");

  /* flush remaining data in all buffers (essentially a copy of code from above) */

  printf("\n %d events in, %d events out, %d events built\n",
         totevts, out_evts, built_evts);
  printf(" %d bad events, %d sub-threshold events\n", badevts, subthreshevts);

  /* if needed, write out stability histograms */
  f_out = fopen("pulstab.rms", "w");
  for (i=0; i<400; i++) {
    sprintf(spname, "%d run %d", i, runInfo->runNumber);
    //if (i==0)  strcat(spname, "; Ge dt [10 ms]");
    //if (i==1)  strcat(spname, "; Ge dt [0.1 s]");
    write_his(his[i], 16384, i, spname, f_out);
  }
  fclose(f_out);

  /* special call to eventprocess() with nChData = -99 as flag
     to finish up processing, write out files, etc. */
  if ((k = eventprocess(Dets, runInfo, -99, ChData)) < 0) return k;

  return  out_evts;
} /* eventbuild() */



/* ------------------------------------------------------- */

#define ERR {printf("ERROR reading rms file at %ld!\n", ftell(file)); return 1;}
/* ======================================================================= */
int read_his(int *his, int idimsp, int spid, char *namesp, FILE *file)
{

  int  i, iy, nsp, dptr, spdir[2];
  int  id, spmode, xlen, expand;
  int  x0, nx, numch;
  char txt[64];
  char *modes[5] = {"shorts", "ints", "floats", "doubles", "bytes"};


  numch = 0;
  iy = spid;

  /* read spectrum ID directory */
  /* directory(ies):   2*(n+1) ints
     int n = number_of_IDs_in_header; int pointer_to_next_directory
     int sp_ID_1    int pointer_to_sp_data_1
     int sp_ID_2    int pointer_to_sp_data_2
     int sp_ID_3    int pointer_to_sp_data_3
     . . .
     int sp_ID_n    int pointer_to_sp_data_n
     [-1            0]  if no more spectra
  */
  //printf(" >> Looking for ID %d\n", iy); fflush(stdout);
  rewind(file);
  while (1) {
    if (1 != fread(&nsp,  sizeof(int), 1, file) ||
        1 != fread(&dptr, sizeof(int), 1, file)) ERR;
    for (i=0; i<nsp; i++) {  // read through directory entries
      if (2 != fread(spdir, sizeof(int), 2, file)) ERR;
      // printf("... found ID %d at %d\n", spdir[0], spdir[1]); fflush(stdout);
      if (spdir[0] == iy) break;  // found it!
    }
    if (spdir[0] == iy || dptr <= 0) break;
    fseek(file, dptr, SEEK_SET);
  }
  if (spdir[0] != iy) {
    printf("Spectrum ID %d not found!\n", iy);
    return 1;
  }

  /* if we are here, then we have found the spectrum we want */
  /* initialize output spectrum to zero */
  memset(his, 0, sizeof(int) * idimsp);
  if (spdir[1] < 1) return 1; // spectrum has zero length / does not exist

  fseek(file, spdir[1], SEEK_SET);
  /* now read the spectrum header */
  if (1 != fread(&id,      sizeof(int), 1, file) || // sp ID again
      1 != fread(&spmode,  sizeof(int), 1, file) || // storage mode (short/int/float...)
      1 != fread(&xlen,    sizeof(int), 1, file) || // x dimension
      1 != fread(&expand,  sizeof(int), 1, file) || // extra space for expansion
      1 != fread(txt,      sizeof(txt), 1, file)) ERR;  // description/name text
  if (id != iy) {
    printf("Spectrum ID mismatch (%d, %d)!\n", iy, id);
    return 1;
  }
  if (spmode != 2) {
    printf("spmode = %d\n", spmode);
    ERR;  // spectrum is not ints
  }

  if (0) printf("Sp ID %d; %d %s; %s\n", iy, xlen, modes[spmode-1], txt);
  if ((numch = xlen) > idimsp) {
    printf("First %d (of %d) chs only taken.", idimsp, numch);
    numch = idimsp;
  }

  i = 0;
  /* now read the actual spectrum data */
  while (i < numch) {
    if (1 != fread(&x0, sizeof(int), 1, file) ||    // starting ch for this chunk of spectrum
        1 != fread(&nx, sizeof(int), 1, file)) ERR; // number of chs in this chunk of spectrum
    if (x0 < 0 || nx < 1) break;   // no more non-zero bins in the spectrum
    if (nx > numch - x0) nx = numch - x0;

    if (nx != fread(his+x0, sizeof(int), nx, file)) ERR;
    i = x0 + nx;
  }
  strncpy(namesp, txt, 8);
  return 0;
#undef ERR
} /* rmsread */

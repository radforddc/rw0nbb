#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MJDSort.h"

/* ---------- dataId types we want to process ---------- */
char  *GoodDataTypes[] =
  {
    "ORGretina4MWaveformDecoder",
    "ORGretina4AWaveformDecoder",
    "-"
  };

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define VERBOSE   0

/*  ------------------------------------------------------------ */

int eventprescan(FILE *f_in, FILE *ps_f_out, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int   i, j, crate=0, slot=0, board, evlen;
  int   dataId[32], goodDataId[32], idNum, board_type;
  int   dataIdGM=0;
  unsigned int  head[2], evtdat[20000];
  short  *signal;

  int totevts = 0, out_evts = 0;
  int module_lu[NCRATES+1][21]; // lookup table to map VME crate&slot into module IDs
  int det_lu[NBDS][16];         // lookup table to map module&chan into detector IDs
  int chan_lu[NBDS][16];        // lookup table to map module&chan into parameter IDs

  int    iss;
  unsigned short ss[3000];
  
  if (runInfo->analysisPass < 0) return 0;

  /* initialize */
  if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
  
  if (0) {
    unsigned short db[2];
    unsigned int   *dd = (unsigned int *) db;
    dd[0] = 1;
    printf("%d %d %d\n", dd[0], db[0], db[1]);
    dd[0] = dd[0]<<8;
    printf("%d %d %d\n", dd[0], db[0], db[1]);
    dd[0] = dd[0]<<8;
    printf("%d %d %d\n", dd[0], db[0], db[1]);
    dd[0] = dd[0]<<8;
    printf("%d %d %d\n", dd[0], db[0], db[1]);
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
        if (strstr(runInfo->decoder[i], "ORGretina4MWaveformDecoder"))  dataIdGM     = dataId[i];
      }
  }
  goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below

  /* start loop over reading events from input file
     ============================================== */

  while (fread(head, sizeof(head), 1, f_in) == 1) {
    board_type = head[0] >> 18;
    evlen = (head[0] & 0x3ffff);
    if (board_type != dataIdGM)                  // if not waveform data, then
      fwrite(head, sizeof(head), 1, ps_f_out);  // copy input header data to output file

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

    /* ========== read in the rest of the event data ========== */
    if (fread(evtdat, 4, evlen-2, f_in) != evlen-2) {
      printf("  No more data...\n");
      break;
    }

    totevts++;
    if (board_type != dataIdGM) {
      fwrite(evtdat, 4, evlen-2, ps_f_out);  // copy input data to output file
      continue;   // (not) GRETINA digitizer
    }

    /* --------------- Gretina4M digitizer ---------------- */
    if (DS0 && evlen > 2000) {
      printf(" ----- Discarding event with length (%d) too big\n", evlen);
      continue;
    }
    if (evlen != 1026) {
      printf("\n -------- ERROR: record length (%d) is not 1026!\n"
             " -------- Is this compressed data??\n", evlen);
      return -1;
    }

    /* ------------ do compression of signal ------------ */
    signal = (short *) evtdat + 28;
    iss = compress_signal(signal, ss, 2*(evlen - 2) - 28);

    // do final write/copy to output
    //evlen = (head[0] & 0x3ffff);
    evlen = 16 + iss/2;               // new length of event, in 4-byte words
    head[0] = (head[0] & 0xfffc0000) + evlen;
    fwrite(head, sizeof(head), 1, ps_f_out);  // write header data to output file
    fwrite(evtdat, 2, 28, ps_f_out);  // write record header to output file
    fwrite(ss, 2, iss, ps_f_out);     // write compressed signal to output file

  } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  return  out_evts;
} /* eventprescan() */

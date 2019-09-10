#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0

void listruns(FILE *f_in, MJDetInfo *detInfo, MJRunInfo *runInfo);

int main(int argc, char **argv) {

  FILE       *f_in, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  int        nDets;
  MJRunInfo  runInfo;
  int        argn=1;
  char       *c, fname[256], line[256];

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname_in\n\n", argv[0]);
    return -1;
  }
  /* open raw data file as input */
  if (strstr(argv[argn], ".lis")) {    // input argument is a list file
    if ((f_lis = fopen(argv[argn],"r")) == NULL) {
      fprintf(stderr, "\n Failed to open list input file %s\n", argv[argn]);
      return 0;
    }
    if (!fgets(fname, sizeof(fname), f_lis)) {
      fprintf(stderr, "\n Failed to read list input file %s\n", argv[argn]);
      return 0;
    }
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
  } else {                                      // use command argument as input file
    strncpy(fname, argv[argn], sizeof(fname));
  }
  while (argn < argc && argv[argn][0] == '-') argn ++;
  if ((f_in = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", fname);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", fname);
  strncpy(runInfo.filename, fname, 256);
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* read file header */
  nDets = decode_runfile_header(f_in, detInfo, &runInfo);
  if (nDets < 1) return 1;

  printf(" Run number: %d in file %s\n", runInfo.runNumber, runInfo.filename);

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file */
    runInfo.analysisPass = 0;
    listruns(f_in, detInfo, &runInfo);
    fclose(f_in);
    if (!f_lis) break;

    /* get next input file name */
    if (!fgets(fname, sizeof(fname), f_lis))  // no more lines in the input list file
      strncpy(fname, "0", strlen(fname));     // special string
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
    if (strlen(fname) < 2) break;  // no more files to process

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
  }

  printf("\n All Done.\n\n");
  return 0;
}

/* ========================================================== */

void listruns(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int           i, evlen, board_type, dataIdRun=0;
  unsigned int  head[2], evtdat[20000];
  static int    totevts=0, current_runNumber = 0;

  
  /* initialize */
  for (i=0; i<runInfo->idNum; i++) {
    if (strstr(runInfo->decoder[i], "ORRunDecoderForRun")) {
      dataIdRun = runInfo->dataId[i];
      printf("dataIdRun = %d %s %d\n", dataIdRun, runInfo->decoder[i], i);
    }
  }
  printf("dataIdRun = %d\n", dataIdRun);

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

    if (board_type == dataIdRun) {
      fread(evtdat, 8, 1, f_in);
      if (head[1] & 0x21) {
        printf("------- START Run %d at %d", evtdat[0], evtdat[1]);
        current_runNumber = evtdat[0];
      }
      continue;
    }

    if (evlen > 10000) {
      printf("\n >>>> ERROR: Event length too long??\n"
             " >>>> This file is probably corruped, ending scan!\n");
      break;
    }
    fseek(f_in, 4*(evlen-2), SEEK_CUR);
    if (++totevts % 200000 == 0) {
      printf(" %8d evts\n", totevts); fflush(stdout);
    }
    continue;
  }

  //fclose(f_out);
  printf(" %8d evt\n", totevts);
  return;
} /* listruns() */

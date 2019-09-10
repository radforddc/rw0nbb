#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>

#include "MJDSort.h"
#include "runBits.h"

#define VERBOSE 0


int main(int argc, char **argv) {

  FILE       *f_in, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  MJRunInfo  runInfo;
  int        nDets, i, argn=1, file, file_is_a_pipe=0;
  char       *c, fname[256], line[256];

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [-s[p0/p1][c0/c1] channum] fname_in\n"
            "    p0/p1: select signals from channel channum flagged as not-pulser/pulser\n"
            "    c0/c1: select signals from channel channum flagged as dirty/clean\n\n", argv[0]);
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
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
  } else {   // using command argument list for input files
    strncpy(fname, argv[argn], sizeof(fname));
  }
  /* open raw data file as input */
  file_is_a_pipe=0;
  if (strstr(fname, "stdin")) {               // use stdin for input data
    f_in = stdin;
    file_is_a_pipe=1;
  } else if (strstr(fname, "Fifo") || strstr(fname, "Pipe") ||
             strstr(fname, "FIFO") || strstr(fname, "PIPE") ||
             strstr(fname, "fifo") || strstr(fname, "pipe")) {
    file_is_a_pipe=1;
    printf("\nAttempting to open FIFO or pipe: %s\n", fname);
    if ((file = open(fname, O_RDONLY)) <= 0 ||
        (f_in = fdopen(file,"r")) == NULL) {
      fprintf(stderr, "\n Failed to open input FIFO or pipe: %s\n", fname);
      return 0;
    }
  } else if ((f_in = fopen(fname,"r")) == NULL) {
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

  /* write out results */
#ifndef QUIET
  for (i=0; i<nDets; i++) {
    if (i%20 == 0)
      printf("#    DetID      pos      name     HiGain GAT Enab Thresh     HVch  MaxV Target"
             "     Pulser times enab  ampl atten  CC\n");
    printf(" %3d  was %2d %8s %9s   %d,%2.2d,%d %4d %3d %6d %3d,%2.2d,%d"
           " %5d %6d %9d %7d %3d %5d %3d %d %3d\n", i,
           detInfo[i].OrcaDetID,       detInfo[i].StrName,
           detInfo[i].DetName,         detInfo[i].crate,
           detInfo[i].slot,            detInfo[i].chanHi,
           detInfo[i].crate*512 + detInfo[i].slot*16 + detInfo[i].chanHi, 
           detInfo[i].HGChEnabled,     detInfo[i].HGTrapThreshold,
           detInfo[i].HVCrate,         detInfo[i].HVCard,
           detInfo[i].HVChan,          detInfo[i].HVMax,
           detInfo[i].HVtarget,
           detInfo[i].pulseHighTime,   detInfo[i].pulseLowTime,
           detInfo[i].pulserEnabled,   detInfo[i].amplitude,
           detInfo[i].attenuated,      detInfo[i].finalAttenuated,
           detInfo[i].CCnum);
  }
#endif

  if (runInfo.dataIdGM == 0 && runInfo.dataIdGA == 0) {
    printf("\n No data ID found for Gretina4M or 4A data!\n");
    return 1;
  }
#ifndef QUIET
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

  if (VERBOSE) {
    printf(" Quickstart: %d\n"
           " Ref   time: %d\n"
           " ORCA setup: %s\n"
           " %d dataIds:\n",
           runInfo.quickStart, runInfo.refTime, runInfo.orcaname, runInfo.idNum);
    for (i=0; i<runInfo.idNum; i++)
      printf("  %2d  %s\n", runInfo.dataId[i], runInfo.decoder[i]);
  }

  int nEnab=0;
  for (i=0; i<nDets; i++)
    if (detInfo[i].HGChEnabled || detInfo[i].LGChEnabled) nEnab++;
  printf("\n  %d of %d detectors are enabled.\n", nEnab, nDets);
#endif

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file, to build and process them */
    /* loop to allow for multiple-pass analysis of the same file */
    runInfo.analysisPass = 0;
    while ((i = eventbuild(f_in, detInfo, &runInfo)) < 0) {  // i >= 0 : normal exit
      if (i == -1) {                                         // i = -1: error
        fclose(f_in);
        fprintf(stderr, "\n ERROR: Scan quit with error\n\n");
        return i;
      }                                                      // i < -1: do another pass
      runInfo.analysisPass++;
      if (!file_is_a_pipe)
        fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go back to start of data in file
    }
    fclose(f_in);
    if (file_is_a_pipe) break;  // assume no more files to read if input was a pipe 

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
      i = eventbuild(f_in, detInfo, &runInfo);
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
  }
#ifndef QUIET
  printf("\n All Done.\n\n");
#endif
  return i;
}

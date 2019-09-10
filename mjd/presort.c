#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"
#include "runBits.h"

#define VERBOSE 0


int eventprescan(FILE *f_in, FILE *f_out, MJDetInfo *detInfo, MJRunInfo *runInfo);

int main(int argc, char **argv) {

  FILE       *f_in, *f_out, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  MJRunInfo  runInfo;
  int        nDets, i, argn=1, nEnab=0;
  char       *c, fname[256], line[256], buf[10000];

  if (argc < 2) {
    fprintf(stderr,
            "\nusage: %s fname_in ... fname_in"
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

  printf("\n# PT_ID  Dig\n");
  for (i=0; i<runInfo.nPT; i++) {
    printf(" %3d  %d,%2.2d,%d\n", i,
           runInfo.PTcrate[i], runInfo.PTslot[i], runInfo.PTchan[i]);
  }

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

  for (i=0; i<nDets; i++)
    if (detInfo[i].HGChEnabled || detInfo[i].LGChEnabled) nEnab++;
  printf("\n  %d of %d detectors are enabled.\n", nEnab, nDets);
#endif

  /* open data output file */
  runInfo.analysisPass = 0;
  sprintf(fname, "DS%d", runInfo.runNumber);
  if ((f_out = fopen(fname, "r"))) {  // output file exists!
    fclose(f_out);
    sprintf(fname+strlen(fname), "_%d", (int) getpid());
    printf("\n  Using output file %s\n\n", fname);
  }
  if (!(f_out = fopen(fname, "w"))) {
    fprintf(stderr, "\n  ERROR: Cannot open output file %s\n\n", fname);
    return 1;
  }
  /* copy input file header to output file */
  rewind(f_in);
  for (i=0; i<=4*runInfo.fileHeaderLen/sizeof(buf); i++) {
    fread(buf, sizeof(buf), 1, f_in);
    fwrite(buf, sizeof(buf), 1, f_out);
  }
  fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go back to start of data in input file
  fseek(f_out, 4*runInfo.fileHeaderLen, SEEK_SET);  //  and in output file

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file, to build and process them */
    i = eventprescan(f_in, f_out, detInfo, &runInfo);
    fclose(f_in);
    if (i == -1) {
      fprintf(stderr, "\n ERROR: Scan quit with error\n\n");
      return i;
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
      i = eventprescan(f_in, f_out, detInfo, &runInfo);
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
  fclose(f_out);
  return i;
}

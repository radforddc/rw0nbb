#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"
#include "runBits.h"

/* ---------------------------------------------
   makes get_run_nums
   = code to list run numbers in a presorted file (in case we need to re-presort from scratch)

   OR makes get_run_times
   = code to list start and end times of runs (e.g. for use in interpolating PSA parameters)
   --------------------------------------------- */


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

#define VERBOSE 0
#define DEBUG   0
/*  ------------------------------------------------------------ */

int eventprescan(FILE *f_in, FILE *f_out, MJDetInfo *detInfo, MJRunInfo *runInfo);
int run = 0;
unsigned int start_time=0, end_time=0, duration=0, last_duration=0;


int main(int argc, char **argv) {

  FILE       *f_in, *f_out, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  MJRunInfo  runInfo;
  int        nDets, i, argn=1;
  char       *c, fname[256], line[256];

  if (argc < 2) {
    fprintf(stderr,
            "\nusage: %s DS_fname_in\n\n", argv[0]);
    return -1;
  }

  /* open output file */
#ifdef GET_NUMS
  if ((f_out = fopen("run_nums.sh", "w")) == NULL) {
    fprintf(stderr, "\n Failed to open output file run_nums.sh\n");
    return 0;
  }
  fprintf(f_out, "# Runs in %s:\n", argv[1]);
#endif
#ifdef GET_TIMES
  if ((f_out = fopen("run_times.txt", "w")) == NULL) {
    fprintf(stderr, "\n Failed to open output file run_times.txt\n");
    return 0;
  }
#endif

#ifdef GET_EXPOSURE
  float  a, b, d, e, g, h, mass[100]={0}, enr_mass[100]={0};
  float  enr_exp=0, tot_exp=0, venr_exp=0, vtot_exp=0, csenr_exp=0, cstot_exp=0;
  int    enab[100]={0}, veto[100]={0}, csel[100]={0}, csel_id[1000], csel_start[1000], csel_end[1000];
  char   name[100][8]={""}, pos[8];
  int    j, k, l, m, n_csel = 0;

  /* read detector mass info */
  if ((f_in = fopen("Detector_masses.csv", "r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file Detector_masses.csv\n");
    return 0;
  }
  i = 0;
  while (i < 100 && fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#') continue;
    if (sscanf(line, "%d,%f,%d,%d,%d,%6s,%7s", &j, &mass[i], &k, &l, &m, &pos, &name[i]) != 7 || i != j) {
      fprintf(stderr, "\n Failed to read input file Detector_masses.csv\n");
      return 0;
    }
    //enab[i] = l;
    enab[i] = 1;
    veto[i] = m;
    mass[i] /= 1000.0;
    enr_mass[i] = k * mass[i];
    i++;
  }
  fclose(f_in);
  printf("%d detector masses read from file Detector_masses.csv\n", i);
  printf("Det  0 : %8s has mass %.4f kg\n", name[0], mass[0]);
  printf("Det %d : %8s has mass %.4f kg\n", i-1, name[i-1], mass[i-1]);

  /* read channel selection info */
  if ((f_in = fopen("veto_only_runs.txt", "r"))) {
    printf("Reading channel selection info from file veto_only_runs.txt\n");
    i = 0;
    while (i < 1000 && fgets(line, sizeof(line), f_in)) {
      if (line[0] == '#') continue;
      c = strstr(line, "ds") + 2;
      sscanf(c, "%5d", &csel_start[i]);
      c = strstr(c, "ds") + 2;
      sscanf(c, "%5d", &csel_end[i]);
      c = strstr(c, " ") + 2;
      sscanf(c, "%d", &csel_id[i]);
      if (csel_id[i] >= 0 && csel_id[i] < 100) i++;
    }
    n_csel = i;
    fclose(f_in);
    printf("Read %d channel selection lines\n", i);
  }

  if ((f_out = fopen("exposure.txt", "w")) == NULL) {
    fprintf(stderr, "\n Failed to open output file exposure.txt\n");
    return 0;
  }
  fprintf(f_out,
          "#             tot_enr   veto_enr   csel_enr    tot_nat   veto_nat   csel_nat    total\n"
          "# duration   exposure   exposure   exposure   exposure   exposure   exposure   exposure     file_name\n"
          "#  [days]     [kg-d]     [kg-d]     [kg-d]     [kg-d]     [kg-d]     [kg-d]     [kg-d]\n");
#endif

  /* get data file name from command arguments) */
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

#ifdef GET_EXPOSURE
  FILE *fp2;
  float dcratio[100];
  for (i=0; i<100; i++) dcratio[i] = 0;

  /* see if there's a file a3.txt, created by count_dc, matching the unput data file */
  strncpy(line, fname, sizeof(line));
  if ((c = strstr(line, "Data")) || (c = strstr(line, "Data"))) {
    sprintf(c, "a3.txt");
    if ((fp2 = fopen(line, "r"))) {
      printf(" ... Reading DC ratio data from %s\n", line);
      while (fgets(line, sizeof(line), fp2)) {
        if (line[0] != '#' && sscanf(line, "%d %f", &i, &a) == 2) dcratio[i] = a;
      }
      fclose(fp2);
    }
  }
#endif

  printf("\n >>> Reading %s\n\n", fname);
  strncpy(runInfo.filename, fname, sizeof(runInfo.filename));
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* get all the required information for the file header */
  nDets = decode_runfile_header(f_in, detInfo, &runInfo);
  if (nDets < 1) return 1;

  if (!runInfo.flashcam)
    printf(" Run number: %d in file %s\n", runInfo.runNumber, runInfo.filename);

  runInfo.analysisPass = 0;
  fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go to start of data in input file

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file, to build and process them */
    i = eventprescan(f_in, f_out, detInfo, &runInfo);
    fclose(f_in);
    if (i == -1) {
      fprintf(stderr, "\n ERROR: Scan quit with error\n\n");
      return i;
    }
#ifdef GET_EXPOSURE
    sscanf(fname+2, "%5d", &j);  // j = run number for current data subset or file
    for (i=0; i<nDets; i++) csel[i] = 0;
    for (i=0; i<n_csel; i++) {
      if (j >= csel_start[i] && j <= csel_end[i] && !veto[csel_id[i]]) {
        csel[csel_id[i]] = 1;
        printf("ChannelSelection rejecting  %s  for DetID %2d\n", fname, csel_id[i]);
      }
    }

    a = b = d = e = g = h = 0;
    for (i=0; i<nDets; i++) {
      float f = enab[i] * (duration - last_duration) / (3600.0*24.0) * detInfo[i].HGChEnabled;
      a += mass[i]     * f;
      b += enr_mass[i] * f;
      d += mass[i]     * f * veto[i];
      e += enr_mass[i] * f * veto[i];
      g += mass[i]     * f * fmax(csel[i], dcratio[i]);
      h += enr_mass[i] * f * fmax(csel[i], dcratio[i]);
    }
    tot_exp += a;
    enr_exp += b;
    vtot_exp += d;
    venr_exp += e;
    cstot_exp += g;
    csenr_exp += h;
    printf( "%8.2f  %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f      %s\n",
            (duration - last_duration)/(3600.0*24.0), b, e, g, a-b, d-e, g-h, a, fname);
    fprintf(f_out,
            "%8.2f  %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f      %s\n",
            (duration - last_duration)/(3600.0*24.0), b, e, g, a-b, d-e, g-h, a, fname);
#endif
  last_duration = duration;

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
#ifdef GET_NUMS
        fprintf(f_out, "# Runs in %s:\n", fname);
#endif
      }
    }
    if (strlen(fname) < 2) break;  // no more files to process

    /* open next file and read header size and run number */
    if ((f_in = fopen(fname,"r")) == NULL) {
      fprintf(stderr, "\n Failed to open input file %s\n", fname);
      return 0;
    }

#ifdef GET_EXPOSURE
    for (i=0; i<100; i++) dcratio[i] = 0;
    /* see if there's a file a3.txt, created by count_dc, matching the unput data file */
    strncpy(line, fname, sizeof(line));
    if ((c = strstr(line, "Data")) || (c = strstr(line, "Data"))) {
      sprintf(c, "a3.txt");
      if ((fp2 = fopen(line, "r"))) {
        printf(" ... Reading DC ratio data from %s\n", line);
        while (fgets(line, sizeof(line), fp2)) {
          if (line[0] != '#' && sscanf(line, "%d %f", &i, &a) == 2) dcratio[i] = a;
        }
        fclose(fp2);
      }
    }
#endif

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

#ifdef GET_TIMES
  printf("\n>>>  %u %u %u %u : start end difference duration (%.2f %.2f hrs)\n",
         start_time, end_time, end_time-start_time, duration,
         (float) (end_time-start_time)/3600.0, duration/3600.0);
  fprintf(f_out, "%u %u %u %u : start end difference duration (%.2f %.2f hrs)\n",
          start_time, end_time, end_time-start_time, duration,
          (float) (end_time-start_time)/3600.0, duration/3600.0);
#endif
#ifdef GET_EXPOSURE
    printf( "#-------------------------------------------------------------------------------------\n"
            "%8.2f  %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f  %s\n",
            duration/(3600.0*24.0), enr_exp, venr_exp, csenr_exp,
            tot_exp-enr_exp, vtot_exp-venr_exp, cstot_exp-csenr_exp, tot_exp, "TOTAL");
    fprintf(f_out,
            "#-------------------------------------------------------------------------------------\n"
            "%8.2f  %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f  %s\n",
            duration/(3600.0*24.0), enr_exp, venr_exp, csenr_exp,
            tot_exp-enr_exp, vtot_exp-venr_exp, cstot_exp-csenr_exp, tot_exp, "TOTAL");
    printf( "#-------------------------------------------------------------------------------------\n"
            "#\n# total clean exposure: %10.2f enriched, %10.2f natural  [kg-days]\n",
            enr_exp - venr_exp - csenr_exp,
            (tot_exp-enr_exp) - (vtot_exp-venr_exp) - (cstot_exp-csenr_exp));
    fprintf(f_out,
            "#-------------------------------------------------------------------------------------\n"
            "#\n# total clean exposure: %10.2f enriched, %10.2f natural  [kg-days]\n",
            enr_exp - venr_exp - csenr_exp,
            (tot_exp-enr_exp) - (vtot_exp-venr_exp) - (cstot_exp-csenr_exp));
#endif
  fclose(f_out);

  printf("\n All Done.\n\n");
  return i;
}

/*  ------------------------------------------------------------ */

int eventprescan(FILE *f_in, FILE *f_out, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int           i, j, n, evlen;
  unsigned int  head[2], evtdat[20000];
  static int    dataId[32], goodDataId[32], idNum, board_type;
  static int    dataIdRun=0;

  static int first = 1;
  static int first_runNumber = 0, current_runNumber = 0;
  static int last_start_time = 0, last_end_time = 0;


  if (runInfo->analysisPass < 0) return -1;
  
  /* initialize */
  current_runNumber = runInfo->runNumber;

  if (first || runInfo->analysisPass > 0) {
    first = 0;
    i = runInfo->firstRunNumber = first_runNumber = current_runNumber = runInfo->runNumber;

    /* identify the data types that we want to decode */
    idNum = runInfo->idNum;
    for (i=0; i<idNum; i++) {
      dataId[i] = runInfo->dataId[i];
      goodDataId[i] = 0;
      for (j=0; GoodDataTypes[j][0]=='O'; j++)
        if (strstr(runInfo->decoder[i], GoodDataTypes[j])) {
          goodDataId[i] = 1;
          if (strstr(runInfo->decoder[i], "ORRunDecoderForRun")) dataIdRun   = dataId[i];
        }
    }
    goodDataId[i] = 0;  // in case we don't find a valid dataID when scanning events below
  }

  /* start loop over reading events from input file
     ============================================== */

  while (1) {

    if (fread(head, sizeof(head), 1, f_in) != 1) break;
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
    if (j >= idNum) {  // the dataID type is not found in the file header
      n = evlen-2;
      if (n > 32) n = 32;
      fread(evtdat, 4, n, f_in);
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
#ifdef GET_NUMS
        fprintf(f_out, "grep %d dir.lis\n", evtdat[0]);
#endif
        current_runNumber = evtdat[0];
        last_start_time   = evtdat[1];
        if (abs(last_start_time - last_end_time) < 5) last_start_time = last_end_time;
      } else if (head[1] & 0x8) {
        if (VERBOSE)
          printf(" -- %8.8x Heartbeat step %d at %d\n", head[1], evtdat[0], evtdat[1]);
      } else if ((head[1] & 0x29) == 0) {
        printf("   ------- END run %d at %d ------- ( %.1f minutes)\n",
               evtdat[0], evtdat[1], (float) (evtdat[1] - last_start_time)/60.0);
        last_end_time   = evtdat[1];
        duration += evtdat[1] - last_start_time;

      } else {
        printf("***** ORRunDecoderForRun %8.8x %d %d\n",
               head[1], evtdat[0], evtdat[1]);
      }
      if (start_time == 0 || start_time > evtdat[1]) start_time = evtdat[1];
      if (end_time == 0 || end_time < evtdat[1]) end_time = evtdat[1];
      /* ------------- end ORRunDecoderForRun end -------------- */

    } else {   // all other event types
      if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) break;
      continue;
    }

   } /* +--+--+--+--+--+--+--+--+ END of reading events +--+--+--+--+--+--+--+--+ */

  if (VERBOSE) printf("FINISHED reading file\n");
  return 1;
} /* eventprescan */

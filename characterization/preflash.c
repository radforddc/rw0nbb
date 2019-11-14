#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 1
#define E_THRESH  1500  // threshold in keV, for speed of processing;
                        // should be less than 1550 to get DE peak for A/E
                        // can be over-ridden by specifying elo != 3000

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1;
  char       *c, data_file_name[256], list_file_name[256], out_file_name[256];
  FILE       *f_in = NULL, *f_file_list = NULL, *f_out;


  if (argc < 3) {
    fprintf(stderr, "\nusage: %s fname_in fname_out [e_lo]\n\n", argv[0]);
    return -1;
  }
  /* open raw data file as input */
  while (argn < argc && argv[argn][0] == '-') argn += 2;
  if (strstr(argv[argn], ".lis")) {
    /* input file name is a .lis file, contaning more than one data file name */
    printf("Reading list of input data files from list file %s\n", list_file_name);
    strncpy(list_file_name, argv[argn], sizeof(list_file_name));
    if (!(f_file_list = fopen(list_file_name, "r"))) {
      fprintf(stderr, "\n Failed to open input file %s\n", list_file_name);
      return 1;
    }
    if (!fgets(data_file_name, sizeof(data_file_name), f_file_list)) {
      fprintf(stderr, "\n Failed to read input file %s\n", list_file_name);
      return 1;
    }
    for (c = data_file_name + strlen(data_file_name); *c < '0'; c--) *c = 0;  // removing trailing junk
  } else {
    strncpy(data_file_name, argv[argn], sizeof(data_file_name));
  }

  f_in = fopen(data_file_name, "r");
  if (f_in == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", data_file_name);
    return 1;
  }
  printf("\n >>> Reading %s\n\n", data_file_name);
  strncpy(runInfo.filename, data_file_name, 256);
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* read file header */
  int nDets = decode_runfile_header(f_in, Dets, &runInfo);
  if (nDets < 1) return 1;
  if (runInfo.flashcam != 1) {
    fprintf(stderr, "\n This program only handles raw flashcam data.\n");
    return 1;
  }

  if (!runInfo.flashcam) {
    printf(" This program is for FlashCan data ONLY. Sorry.\n");
    return 1;
  }

  /* open output file for presorted data */
  strncpy(out_file_name, argv[2], sizeof(out_file_name) - 8);
  strncat(out_file_name, ".fciops", 8);
  f_out = fopen(out_file_name, "w");
  if (f_out == NULL) {
    fprintf(stderr, "\n Failed to open output file %s\n", out_file_name);
    return 1;
  }

/* ---------------------------------------------------- */

  int totevts=0, out_evts=0;

  unsigned int   evtdat[20000];
  short          *signal;
  unsigned short sigc[2048];
  int    i, j, k, t70, t100, sig_len, his[8192] = {0};

  int elo = 3000;
  // see if low-energy limit is defined in the command line
  if (runInfo.argc > 3) elo = atoi(runInfo.argv[3]);
  printf("\n e_trapmax limit is %d\n\n", elo);

  // start loop over reading events from input file
  while (1) {
    int evlen = 0, last_read = 0;

    while ((last_read = fread(&k, sizeof(int), 1, f_in)) == 1 && k < 1200) {
      if (k > 1) {
        if ((last_read = fread(evtdat, k, 1, f_in)) != 1) break;
      }
    }
    evlen = k/4 + 2;

    /*  read in the rest of the event data  */
    if (last_read != 1 || fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
      /* if no more data, *break* to end
         or read a new file name from the list file to open a new data file */
      printf("\n  No more data in file %s\n\n", data_file_name);
      fclose(f_in);
      if (!f_file_list) break;
      if (!fgets(data_file_name, sizeof(data_file_name), f_file_list)) {
        printf("   End of list file %s\n\n", list_file_name);
        break;
      }
      for (c = data_file_name + strlen(data_file_name); *c < '0'; c--) *c = 0;  // removing trailing junk
      f_in = fopen(data_file_name, "r");
      if (f_in == NULL) {
        fprintf(stderr, "     Failed to open new input file %s\n\n", data_file_name);
        break;
      }
      fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go to start of data in input file
      printf(" >>> Reading %s\n\n", data_file_name);
      continue;
    }
    if (++totevts % 50000 == 0) {   
      printf(" %8d evts in, %d out\n", totevts, out_evts); fflush(stdout);
    }

    sig_len = 2*(evlen-2);
    signal = (short *) evtdat;
    for (i=0; i < sig_len; i++) signal[i] = ((unsigned short) signal[i]) / 4 - 8192;

    /* sticky-bit fix */
    int d = 128;  // basically a sensitivity threshold; max change found will be d/2
    for (i=20; i<2000; i++) {
      // calculate second derivatives
      int dd0 = abs((int) signal[i+1] - 2*((int) signal[i]) + (int) signal[i-1]);
      int dd1 = abs((int) signal[i+2] - 2*((int) signal[i+1]) + (int) signal[i]);
      if (dd0 > d && dd0 > dd1 && dd0 > abs(signal[i+1] - signal[i-1])) {
        // possible occurrence; make sure it's not just high-frequency noise
        for (k=i-8; k<i+8; k++) {
          if (k==i-1 || k==i || k == i+1) continue;
          dd1 = abs((int) signal[k+1] - 2*((int) signal[k]) + (int) signal[k-1]);
          if (dd0 < dd1*3) break;
        }
        if (k<i+8) continue;
        dd0 = (int) signal[i+1] - 2*((int) signal[i]) + (int) signal[i-1];
        j = lrintf((float) dd0 / (float) d);
        printf("Fixing sticky bit in signal %d, t=%d, change %d\n",
               out_evts-1, i, j*d/2);
        signal[i] += j*d/2;
        // break;
      }
    }

    /* find t100 and t70*/
    t100 = 700;                 // FIXME? arbitrary 700?
    for (i = t100+1; i < sig_len; i++)
      if (signal[t100] < signal[i]) t100 = i;
    /* get mean baseline value */
    int bl = 0;
    for (i=100; i<150; i++) bl += signal[i];
    bl /= 50;
    for (t70 = t100-1; t70 > 700; t70--)
      if ((signal[t70] - bl) <= (signal[t100] - bl)*7/10) break;

    int e = (signal[t100] - bl)/2;
    if (e < elo ) continue;
    if (e < 8192) his[e]++;
    out_evts++;

    int start = t70 - 1000;    // starting point in signal to write out
    int len = 2048;            // length signal to write out
    if (start < 0) start = 0;
    if (len > sig_len)   len = sig_len;
    if (start+len > sig_len) start = sig_len - len;

    sigc[0] = sigc[1] = len;
    if (1) {
      // compress signal before writing to output file
      sigc[0] = compress_signal(signal+start, sigc+2, len);
      fwrite(sigc, sizeof(short), sigc[0]+2, f_out);
    } else {
      // write uncompressed signal
      fwrite(sigc, sizeof(short), 2, f_out);
      fwrite(signal+start, sizeof(short), len, f_out);
    }
  }
  fclose(f_out);
  printf(" %8d evts in, %d saved\n", totevts, out_evts);

  f_out = fopen("preflash.sec", "w");
  fwrite(his, sizeof(int), 8192, f_out);
  fclose(f_out);

  return 0;
}

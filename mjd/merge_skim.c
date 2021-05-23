#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 1

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  CTCinfo    CTC;
  int        argn=1, nsd_start = 0;
  int        adjust_aoe_pos = 1;     // change 1 to 0 to NOT change skimmed A/E values to a common peak position
  float      aoe_pos1[200] = {1000}, aoe_pos[200]={1000}, a, b, pos;
  char       *c, psa_fname[256], line[256], ds_fname[256];
  FILE       *f_in, *f_in2, *f_out, *fp;

  SavedData2 **sd;
  int     nsd = 0;           // number of data to be read from one file
  int     sdchunk = 1000000; // number of SavedData events to malloc at one time
  int     numsd = 0, maxsd = sdchunk;  // number of saved data, and pointer to saved data id
  int     i, j, sd_version = 1;
  double  calib[200] = {0};


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname1_in fname2_in [fname3_in ...] \n\n", argv[0]);
    return -1;
  }
  /* malloc initial space for SavedData */
  if ((sd = malloc(sdchunk*sizeof(*sd))) == NULL ||
      (sd[0] = malloc(sdchunk*sizeof(**sd))) == NULL) {
    printf("ERROR in merge_skim.c; cannot malloc SavedData!\n");
    exit(-1);
  }
  for (i=1; i<sdchunk; i++) sd[i] = sd[i-1] + 1;
  nsd = sdchunk;

  /* open skim data file as input */
  while (argn < argc && (f_in = fopen(argv[argn], "r"))) {

    /* read charge-trapping correction values corresponding to new skim file */
    strncpy(ds_fname, argv[argn], sizeof(ds_fname));
    if ((c = strstr(ds_fname, "/skim.dat"))) {
      *c = 0;
    } else {
      printf("\n ERROR: No directory found in %s; c = %d\n", ds_fname, (int) c);
      exit(-1);
    }
    sprintf(CTC.ctc_fname, "%s/ctc.input", ds_fname);
    /* read energy correction factors from ctc.input */
    if (!CTC_info_read(&runInfo, &CTC)) {
      printf("\n ERROR: No charge-trapping correction data read. Does %s exist?\n", CTC.ctc_fname);
      exit(-1);
    }

    // read detector and run info from the saved skim data, from f_in
    fread(&nsd, sizeof(int), 1, f_in);
    if (nsd == -2) {
      sd_version = 2;
      fread(&nsd, sizeof(int), 1, f_in);
    }
    fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);  // NOTE that this overwrites energy calibs
    fread(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_in);
    if (runInfo.idNum == 0) {
      runInfo.flashcam = 1;
      fread(&(runInfo.flashcam), 8*sizeof(int), 1, f_in);
    }

    /* read energy calibration gains corresponding to new skim file */
    /* this part must be done after reading detector info from skim.dat,
       since that would override these calibs. */
    sprintf(line, "%s/gains.input", ds_fname);
    if ((fp = fopen(line, "r"))) {
      printf("\n Reading energy calibrations from %s, argn = %d\n", line, argn);
      double g1 = 0.5, g2=1.5;
      j = 0;
      while (fgets(line, sizeof(line), fp)) {
        if (*line != '#' &&
            sscanf(line, "%d %lf %lf", &i, &g1, &g2) == 3 &&
            i >= 0 && i < 100) {
          Dets[i].HGcalib[0] = g1;
          Dets[i].LGcalib[0] = g2;
          if (argn == 1 || calib[i] ==  0) {
            calib[i] = g1;
            calib[i+100] = g2;
            printf("i, calib: %d %.4f\n", i, calib[i]);
          }
          if (calib[i+100] == 0) calib[i+100] = g2;
          j++;
        }
      }
      printf("\n %d energy calibrations read, argn = %d\n", j, argn);
      fclose(fp);
    } else {
      printf("ERROR: Cannot open %s !\n", line);
      exit(-1);
    }

    if (adjust_aoe_pos) {
      strncpy(psa_fname, argv[argn], sizeof(psa_fname)-5);
      if ((c = strstr(psa_fname, "skim.dat"))) {
        strncpy(c, "psa.input", 12);              // replace skim.dat in psa_fname with psa.input
        if (!(f_in2 = fopen(psa_fname, "r"))) {
          printf("Error: Cannot open file %s. Has script1 been run?\n", psa_fname);
          return 0;
        }
        /* read a/e positions from psa.input */
        while (fgets(line, 256, f_in2)) {
          if (line[0] == '#') continue;  // # header line
          if (sscanf(line, "%d %f %f %f", &i, &a, &b, &pos) < 4) continue;
          if (i >= 0 && i < 200) aoe_pos[i] = pos;
        }
        fclose(f_in2);
        printf("A/E positions read from %s\n", psa_fname);
        if (argn == 1) {
          for (i=0; i<200; i++) aoe_pos1[i] = aoe_pos[i];
          printf("   ... and stored for position matching.\n");
        }
      }
    }
    printf("Reading skim file %s\n", argv[argn]);

    // read saved skim data from f_in
    printf("Skim data mode = %d;  %d detectors, %d skimmed events\n",
           sd_version, runInfo.nGe, nsd);
    if (sd_version == 1) {
      printf("Error: This version of merge_skim works only for newer type skim files.\n");
      exit(-1);
    } else {   // sd_version == 2
      while (numsd + nsd >= maxsd) {
        //printf("0  nsd: %d  numsd: %d    maxsd: %d\n", nsd, numsd, maxsd); fflush(stdout);
        fread(sd[numsd], sizeof(**sd), maxsd - numsd, f_in);
        nsd -= maxsd - numsd;
        numsd = maxsd;
        if ((sd = realloc(sd, (maxsd + sdchunk)*sizeof(*sd))) == NULL ||
            (sd[maxsd] = malloc(sdchunk*sizeof(**sd))) == NULL) {
          printf("ERROR in skim.c; cannot realloc SavedData! nsd = %d\n", nsd);
          exit(-1);
        }
        maxsd += sdchunk;
        //printf("1  nsd: %d  numsd: %d    maxsd: %d\n", nsd, numsd, maxsd); fflush(stdout);
        for (i=numsd+1; i<maxsd; i++) sd[i] = sd[i-1]+1;
      }
      //printf("2  nsd: %d  numsd: %d    maxsd: %d\n", nsd, numsd, maxsd); fflush(stdout);
      fread(sd[numsd], sizeof(**sd), nsd, f_in);
      numsd += nsd;
      //printf("3  nsd: %d  numsd: %d    maxsd: %d\n", nsd, numsd, maxsd); fflush(stdout);
    }
    printf(" Skim is data from runs starting at number %d from file %s\n"
           "  Total events now %d\n",
           runInfo.runNumber, runInfo.filename, numsd); 

    if (argn > 0) {
      printf("Adjusting E values by differences in gains (e.g. ID 1 %f -> %f\n", calib[1], Dets[1].HGcalib[0]);
      for (i = nsd_start; i < numsd; i++) {
        if ((j = sd[i]->chan) >= 0 && j < 100 && calib[j] > 0.1) sd[i]->e *= Dets[j].HGcalib[0] / calib[j];
        if (              j >= 100 && j < 200 && calib[j] > 0.1) sd[i]->e *= Dets[j-100].LGcalib[0] / calib[j];
        //if (j >= 0 && j < 200 && calib[j] < 0.1) printf("i, j, calib: %d %d %.4f\n", i, j, calib[j]);
      }
    }
    if (argn > 1 && adjust_aoe_pos) {
      printf("Adjusting A/E values by differences in positions\n");
      for (i = nsd_start; i < numsd; i++) {
        if ((j = sd[i]->chan) >= 0 && j < 200) sd[i]->a_over_e += aoe_pos1[j] - aoe_pos[j];
      }
    }

    argn++;
    nsd_start = numsd;
    fclose(f_in);
  }
  // save skim data (SavedData) to disk
  if (adjust_aoe_pos) {
    printf("\n >>> Writing skim file merged_adj_skim.dat\n"); fflush(stdout);
    f_out = fopen("merged_adj_skim.dat", "w");
  } else {
    printf("\n >>> Writing skim file merged_skim.dat\n"); fflush(stdout);
    f_out = fopen("merged_skim.dat", "w");
  }

  i = -2;
  fwrite(&i, sizeof(int), 1, f_out);    // flag to tell reading programs to use SavedData2 instead of SavedData
  fwrite(&numsd, sizeof(int), 1, f_out);
  for (i = 0; i < 100; i++) { // replace energy calibrations with those from the first merged skim
    if (calib[i] != 0)     Dets[i].HGcalib[0] = calib[i];
    if (calib[i+100] != 0) Dets[i].LGcalib[0] = calib[i+100];
  }
  fwrite(&Dets[0], sizeof(Dets[0]), NMJDETS, f_out);
  if (runInfo.flashcam) {
    runInfo.idNum = 0;
    fwrite(&runInfo, sizeof(runInfo), 1, f_out);  // for backwards compatibility
  } else {
    fwrite(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_out);
  }

  for (i = 0; i < numsd; i += sdchunk) {
    if (numsd - i < sdchunk) {
      fwrite(sd[i], sizeof(**sd), numsd - i, f_out);
    } else {
      fwrite(sd[i], sizeof(**sd), sdchunk, f_out);
    }
  }
  fclose(f_out);
  printf(" Wrote %d skimmed events of size %lu\n\n", numsd, sizeof(**sd));

  return 0;
}

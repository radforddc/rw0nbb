#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>

/*
 * code to average values in a set of psa.input files
 * to use: cd <path_to_data_set>/cal; psa_getmean
 */

#include "MJDSort.h"
int PSA_read(char *fn, PSAinfo *PSA);
int PSA_write(char *fn, PSAinfo *PSA);

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  int      i, n_psa_chan[200] = {0};
  PSAinfo  PSA1, PSA2;
  char     filename[256];
  char     calib_subset_dir[256], cwd[FILENAME_MAX];
  struct dirent *de;  // Pointer for directory entry
  DIR           *dp;


  /* initialize PSA1 to zero */
  for (i=0; i<200; i++) {
    PSA1.ae_dt_slope[i]    = 0;
    PSA1.ae_e_slope[i]     = 0;
    PSA1.ae_pos[i]         = 0;
    PSA1.ae_cut[i]         = 0;
    PSA1.dcr_dt_slope[i]   = 0;
    PSA1.dcr_lim[i]        = 0;
    PSA1.lamda_dt_slope[i] = 0;
    PSA1.lamda_lim[i]      = 0;
    PSA1.lq_dt_slope[i]    = 0;
    PSA1.lq_lim[i]         = 0;
  }

  /* get current working directory */
  getcwd(cwd, FILENAME_MAX);

  /* get lists of calibration directories */
  if ((dp = opendir(".")) == NULL) {
    printf("Could not open current directory" );
    return 0;
  }
  while ((de = readdir(dp)) != NULL) {
    if (strstr(de->d_name, "ds"))  {
      strncpy(calib_subset_dir, de->d_name, sizeof(calib_subset_dir));
    
      /* read A/E, DCR, and lamda values from psa.input */
      sprintf(filename, "%s/psa.input", calib_subset_dir);
      if (PSA_read(filename, &PSA2)) {
        printf("\n ERROR: No PSA data read. Does %s exist?\n", filename);
        return 0;
      }
      /* add PSA2 to PSA1 */
      for (i=0; i<200; i++) {
        if (PSA2.ae_pos[i]!= 500 || PSA2.ae_cut[i] != 500) {
          PSA1.ae_dt_slope[i]    += PSA2.ae_dt_slope[i]   ;
          PSA1.ae_e_slope[i]     += PSA2.ae_e_slope[i]    ;
          PSA1.ae_pos[i]         += PSA2.ae_pos[i]        ;
          PSA1.ae_cut[i]         += PSA2.ae_cut[i]        ;
          PSA1.dcr_dt_slope[i]   += PSA2.dcr_dt_slope[i]  ;
          PSA1.dcr_lim[i]        += PSA2.dcr_lim[i]       ;
          PSA1.lamda_dt_slope[i] += PSA2.lamda_dt_slope[i];
          PSA1.lamda_lim[i]      += PSA2.lamda_lim[i]     ;
          PSA1.lq_dt_slope[i]    += PSA2.lq_dt_slope[i]   ;
          PSA1.lq_lim[i]         += PSA2.lq_lim[i]        ;
          n_psa_chan[i]++;
        }
      }
    }
  }
  closedir(dp);

  /* compute mean values */
  for (i=0; i<200; i++) {
    if (n_psa_chan[i] > 1) {
      PSA1.ae_dt_slope[i]    /= (float) n_psa_chan[i];
      PSA1.ae_e_slope[i]     /= (float) n_psa_chan[i];
      PSA1.ae_pos[i]         /= (float) n_psa_chan[i];
      PSA1.ae_cut[i]         /= (float) n_psa_chan[i];
      PSA1.dcr_dt_slope[i]   /= (float) n_psa_chan[i];
      PSA1.dcr_lim[i]        /= (float) n_psa_chan[i];
      PSA1.lamda_dt_slope[i] /= (float) n_psa_chan[i];
      PSA1.lamda_lim[i]      /= (float) n_psa_chan[i];
      PSA1.lq_dt_slope[i]    /= (float) n_psa_chan[i];
      PSA1.lq_lim[i]         /= (float) n_psa_chan[i];
    } else if (n_psa_chan[i] < 1) {
      PSA1.ae_pos[i]         = 500;
      PSA1.ae_cut[i]         = 500;
      PSA1.lq_lim[i]         = -1;
    }
  }

  // write PSA data to psa.output
  sprintf(filename, "psa_mean.input");
  PSA_write(filename, &PSA1);

  return 0;
}


/*  ------------------------------------------------------------ */

int PSA_read(char *fn, PSAinfo *PSA) {

  int    i, n;
  char   line[256];
  float  a, b, c, d, e, f, g, h, s, q;
  FILE   *f_in;

  /*
    read A/E and DCR data from file psa.input
  */
  for (i=0; i<200; i++) {
    PSA->ae_dt_slope[i]    = 0;
    PSA->ae_e_slope[i]     = 0;
    PSA->ae_pos[i]         = 500;
    PSA->ae_cut[i]         = 500;
    PSA->dcr_dt_slope[i]   = 0;
    PSA->dcr_lim[i]        = 0;
    PSA->lamda_dt_slope[i] = 0;
    PSA->lamda_lim[i]      = 0;
    PSA->lq_dt_slope[i]    = 0;
    PSA->lq_lim[i]         = 0;
  }

  if (!(f_in = fopen(fn, "r"))) return 1;
  printf("\n Reading PSA values from file %s\n", fn);
  i = -1;

  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#' || strlen(line) < 4) continue;
    e = 0; f = 1; s = 0; q = -1;
    if ((n=sscanf(line,
                  "%d %f %f %f %f %f %f %f %f %f %f",
                  &i,&a,&b,&c,&d,&e,&f,&g,&h,&s,&q)) < 9 ||
        i < 0 || i > 199) {
      printf(" ERROR; reading of psa.input file failed! line: %s\n", line);
      fclose(f_in);
      return 1;
    }
    PSA->ae_dt_slope[i]    = a;
    PSA->ae_e_slope[i]     = b;
    PSA->ae_pos[i]         = c;
    PSA->ae_cut[i]         = d;
    PSA->dcr_dt_slope[i]   = e;
    PSA->dcr_lim[i]        = f;
    PSA->lamda_dt_slope[i] = g;
    PSA->lamda_lim[i]      = h;
    PSA->lq_dt_slope[i]    = s;
    PSA->lq_lim[i]         = q;
  }
  printf("Last channel read: %d\n", i);
  fclose(f_in);
  return 0;
} /* PSA_read */

/* ---------------------------------------- */

int PSA_write(char *fn, PSAinfo *PSA) {

  int   i;
  FILE  *f_out;
  int   nGe = 57; 

  /*
    write A/E and DCR data to file psa.output
  */
  if (!(f_out = fopen(fn, "w"))) {
    printf("\n ERROR: Cannot open file %st for writing!\n", fn);
    return 1;
  }
  printf("\n Writing PSA values to file %s\n", fn);

  for (i=0; i<200; i++) {
    if (i%100 == 0)
      fprintf(f_out,
              "#Chan  ae_dt_slope ae_e_slope  ae_pos  ae_cut  dcr_dt_slope dcr_lim  lamda_dt_slope lamda_lim  lq_dt_slope lq_lim\n");
    if (i%100 < nGe)
      fprintf(f_out, "%5d  %11.2f %10.2f %7.2f %7.2f  %12.2f %7.2f  %14.2f %9.2f %12.2f %6.1f\n",
              i, PSA->ae_dt_slope[i], PSA->ae_e_slope[i], PSA->ae_pos[i], PSA->ae_cut[i],
              PSA->dcr_dt_slope[i], PSA->dcr_lim[i], PSA->lamda_dt_slope[i], PSA->lamda_lim[i], PSA->lq_dt_slope[i], PSA->lq_lim[i]);
  }

  fclose(f_out);
  return 0;
} /* PSA_write */

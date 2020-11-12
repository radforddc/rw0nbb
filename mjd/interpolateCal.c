#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>

#include "MJDSort.h"

#define VERBOSE 0

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS], Dets1[NMJDETS], Dets2[NMJDETS];
  MJRunInfo  runInfo;
  CTCinfo    CTC, CTC1, CTC2;
  PSAinfo    PSA, PSA1, PSA2;
  PZinfo     PZI, PZI1, PZI2;
  int        i, nsd, n_ds=0;
  char       ds1[256] = "", ds2[256] = "", cwd[FILENAME_MAX];
  char       fn1[256], fn2[256], line[1024], line2[1024];
  double     g1 = 0.5, g2=1.5;
  FILE       *f_in, *f_out;
  struct dirent *de;  // Pointer for directory entry
  DIR           *dp;


  /* get current working directory */
  getcwd(cwd, FILENAME_MAX);

  /* get lists of calibration directories */
  if ((dp = opendir(".")) == NULL) {
    printf("Could not open current directory" );
    return 0;
  }
  while ((de = readdir(dp)) != NULL) {
    if (strstr(de->d_name, "ds"))  {
      strncpy(ds1, ds2, sizeof(ds1));
      n_ds++;
      strncpy(ds2, de->d_name, sizeof(ds2));
      if (n_ds < 2) continue;

      printf("\nReading calib input files from %s and %s\n", ds1, ds2);

      /* start by reading just the run info from one of the skim.dat files */
      sprintf(fn1, "%s/skim.dat", ds1);
      if (!(f_in = fopen(fn1, "r"))) {
        fprintf(stderr, "\n Failed to open file %s; does directory exist?\n", fn1);
        return 0;
      }
      fread(&nsd, sizeof(int), 1, f_in);
      if (nsd == -2) fread(&nsd, sizeof(int), 1, f_in);
      fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
      fread(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_in);
      if (runInfo.idNum == 0) {
        runInfo.flashcam = 1;
        fread(&(runInfo.flashcam), 8*sizeof(int), 1, f_in);
      }
      fclose(f_in);
  
      /* read gains from gains.input */
      memcpy(&Dets1[0], &Dets[0], sizeof(Dets[0])*NMJDETS);
      sprintf(fn1, "%s/gains.input", ds1);
      if (!(f_in = fopen(fn1, "r"))) {
        fprintf(stderr, "\n Failed to open file %s; does directory exist?\n", fn1);
        return 0;
      }
      printf("Reading energy calibrations from %s\n",fn1);
      while (fgets(line, sizeof(line), f_in)) {
        if (*line != '#' &&
            sscanf(line, "%d %lf %lf", &i, &g1, &g2) == 3 &&
            i >=0 && i < runInfo.nGe) {
          Dets1[i].HGcalib[0] = g1;
          Dets1[i].LGcalib[0] = g2;
        }
      }
      fclose(f_in);

      memcpy(&Dets2[0], &Dets[0], sizeof(Dets[0])*NMJDETS);
      sprintf(fn2, "%s/gains.input", ds2);
      if (!(f_in = fopen(fn2, "r"))) {
        fprintf(stderr, "\n Failed to open file %s; does directory exist?\n", fn2);
        return 0;
      }
      printf("Reading energy calibrations from %s\n",fn2);
      while (fgets(line, sizeof(line), f_in)) {
        if (*line != '#' &&
            sscanf(line, "%d %lf %lf", &i, &g1, &g2) == 3 &&
            i >=0 && i < runInfo.nGe) {
          Dets2[i].HGcalib[0] = g1;
          Dets2[i].LGcalib[0] = g2;
        }
      }
      fclose(f_in);

      /* copy data from the first of the two directories */
      sprintf(line, "mv PZ.input PZsave.input; mv ctc.input ctcsave.input; mv psa.input psasave.input");
      sprintf(line2, "mv PZsave.input PZ.input; mv ctcsave.input ctc.input; mv psasave.input psa.input");
      printf(">>> %s\n", line);
      system(line);
      sprintf(line, "cp %s/PZ.input .; cp %s/ctc.input .; cp %s/psa.input . ", ds1, ds1, ds1);
      printf(">>> %s\n", line);
      system(line);

      /* read PZ correction info from PZ.input */
      if (!PZ_info_read(&runInfo, &PZI1)) {
        printf("\nERROR: No PZ.input data read. Does %s/PZ.input exist?\n", ds1);
        system(line2);
        return 0;
      }
      /* read energy correction factors from ctc.inpt */
      if (!CTC_info_read(&runInfo, &CTC1)) {
        printf("\nERROR: No ctc.input data read. Does %s/ctc.input exist?\n", ds1);
        system(line2);
        return 0;
      }
      /* read A/E, DCR, and lamda values from psa.input */
      if (!PSA_info_read(&runInfo, &PSA1)) {
        printf("\nERROR: No PSA data read. Does %s/psa.input exist?\n", ds1);
        system(line2);
        return 0;
      }
      system(line);

      /* copy data from the second of the two directories */
      sprintf(line, "cp %s/PZ.input .; cp %s/ctc.input .; cp %s/psa.input . ", ds2, ds2, ds2);
      printf(">>> %s\n", line);
      system(line);

      /* read PZ correction info from PZ.input */
      if (!PZ_info_read(&runInfo, &PZI2)) {
        printf("\nERROR: No PZ.input data read. Does %s/PZ.input exist?\n", ds2);
        system(line2);
        return 0;
      }
      /* read energy correction factors from ctc.inpt */
      if (!CTC_info_read(&runInfo, &CTC2)) {
        printf("\nERROR: No ctc.input data read. Does %s/ctc.input exist?\n", ds2);
        system(line2);
        return 0;
      }
      /* read A/E, DCR, and lamda values from psa.input */
      if (!PSA_info_read(&runInfo, &PSA2)) {
        printf("\nERROR: No PSA data read. Does %s/psa.input exist?\n", ds2);
        system(line2);
        return 0;
      }
      printf(">>> %s\n", line2);
      system(line2);

      /* ---------------------- do interpolation ------------------*/
      memcpy(&PZI, &PZI1, sizeof(PZI));
      memcpy(&CTC, &CTC1, sizeof(CTC));
      memcpy(&PSA, &PSA1, sizeof(PSA));

      double dt = (PSA.e_ctc_rise[1] + PSA.e_ctc_flat[1]) / 200.0;    // μs
      double s1=0, s2=0;

      for (i=0; i<runInfo.nGe; i++) {
        PZI.tau[i]            = (PZI1.tau[i]            + PZI2.tau[i]           ) / 2.0;
        PZI.baseline[i]       = (PZI1.baseline[i]       + PZI2.baseline[i]      ) / 2.0;
        PZI.frac2[i]          = (PZI1.frac2[i]          + PZI2.frac2[i]         ) / 2.0;
        PZI.bl_rms[i]         = (PZI1.bl_rms[i]         + PZI2.bl_rms[i]        ) / 2.0;

        CTC.e_dt_slope[i]     = (CTC1.e_dt_slope[i]     + CTC2.e_dt_slope[i]    ) / 2.0;
        CTC.e_lamda_slope[i]  = (CTC1.e_lamda_slope[i]  + CTC2.e_lamda_slope[i] ) / 2.0;
        CTC.e_lamda_gain[i]   = (CTC1.e_lamda_gain[i]   + CTC2.e_lamda_gain[i]  ) / 2.0;
        CTC.best_dt_lamda[i]  = (CTC1.best_dt_lamda[i]  + CTC2.best_dt_lamda[i] ) / 2;

        PSA.ae_dt_slope[i]    = (PSA1.ae_dt_slope[i]    + PSA2.ae_dt_slope[i]   ) / 2.0;
        PSA.ae_e_slope[i]     = (PSA1.ae_e_slope[i]     + PSA2.ae_e_slope[i]    ) / 2.0;
        PSA.ae_pos[i]         = (PSA1.ae_pos[i]         + PSA2.ae_pos[i]        ) / 2.0;
        PSA.ae_cut[i]         = (PSA1.ae_cut[i]         + PSA2.ae_cut[i]        ) / 2.0;
        PSA.dcr_dt_slope[i]   = (PSA1.dcr_dt_slope[i]   + PSA2.dcr_dt_slope[i]  ) / 2.0;
        PSA.dcr_lim[i]        = (PSA1.dcr_lim[i]        + PSA2.dcr_lim[i]       ) / 2.0;
        PSA.lamda_dt_slope[i] = (PSA1.lamda_dt_slope[i] + PSA2.lamda_dt_slope[i]) / 2.0;
        PSA.lamda_lim[i]      = (PSA1.lamda_lim[i]      + PSA2.lamda_lim[i]     ) / 2.0;
        PSA.lq_dt_slope[i]    = (PSA1.lq_dt_slope[i]    + PSA2.lq_dt_slope[i]   ) / 2.0;
        PSA.lq_lim[i]         = (PSA1.lq_lim[i]         + PSA2.lq_lim[i]        ) / 2.0;

        Dets[i].HGcalib[0]    = (Dets1[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI1.tau[i])) +
                                 Dets2[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI2.tau[i]))) / 2.0;  // CHECKME
        Dets[i].LGcalib[0]    = (Dets1[i].LGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI1.tau[i])) +
                                 Dets2[i].LGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI2.tau[i]))) / 2.0;
        if (0) {
          printf("%3d %13.10f %13.10f   %13.10f\n",
                 i, Dets[i].HGcalib[0] - (Dets1[i].HGcalib[0] + Dets2[i].HGcalib[0])/2.0,
                 (Dets1[i].HGcalib[0] * exp((dt/PZI1.tau[i] - dt/PZI.tau[i])) +
                  Dets2[i].HGcalib[0] * exp((dt/PZI2.tau[i] - dt/PZI.tau[i]))) / 2.0 -
                 (Dets1[i].HGcalib[0] + Dets2[i].HGcalib[0])/2.0,
                 exp((dt/PZI1.tau[i] - dt/PZI.tau[i])));
          if (i!= 32) {
            s1 += fabs(Dets[i].HGcalib[0] - (Dets1[i].HGcalib[0] + Dets2[i].HGcalib[0])/2.0);
            s2 += fabs((Dets1[i].HGcalib[0] * exp((dt/PZI1.tau[i] - dt/PZI.tau[i])) +
                        Dets2[i].HGcalib[0] * exp((dt/PZI2.tau[i] - dt/PZI.tau[i]))) / 2.0 -
                       (Dets1[i].HGcalib[0] + Dets2[i].HGcalib[0])/2.0);
          }
        }
        /* check for big changes that imply veto-only status for ββ runs */
        // gain
        if (Dets[i].DetName[0] == 'P') {
          if (0.0002 < fabs(Dets1[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI1.tau[i])) -
                            Dets2[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI2.tau[i]))))
            printf(" ---------- Detector  %2d %s %s veto-only due to gain change!  %6.2f ppt\n",
                   i, Dets[i].StrName, Dets[i].DetName,
                   1000.0 * (Dets1[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI1.tau[i])) -
                             Dets2[i].HGcalib[0] * exp((dt/PZI.tau[i] - dt/PZI2.tau[i])))/Dets[i].HGcalib[0]);
          // A/E
          if (fabs(PSA1.ae_pos[i] - PSA2.ae_pos[i]) >
              (PSA1.ae_pos[i] - PSA1.ae_cut[i] + PSA1.ae_pos[i] - PSA1.ae_cut[i]) / 3.0)
            printf(" ---------- Detector  %2d %s %s veto-only due to A/E  change!  %6.2f, pos-cut = %5.2f\n",
                   i, Dets[i].StrName, Dets[i].DetName,
                   PSA1.ae_pos[i] - PSA2.ae_pos[i],  
                   (PSA1.ae_pos[i] - PSA1.ae_cut[i] + PSA1.ae_pos[i] - PSA1.ae_cut[i])/2.0);
          // DCR
          if (fabs(PSA1.dcr_lim[i] - PSA2.dcr_lim[i]) > 2.0*8.0/2.355)   // ~ 2σ
            printf(" ---------- Detector  %2d %s %s veto-only due to DCR  change!  %6.2f  > %4.2f\n",
                   i, Dets[i].StrName, Dets[i].DetName, PSA1.dcr_lim[i] - PSA2.dcr_lim[i], 2.0*8.0/2.355);
          // LQ
          if (fabs(PSA1.lq_lim[i] - PSA2.lq_lim[i]) > 2.0*5.0/2.355)     // ~ 2σ
            printf(" ---------- Detector  %2d %s %s veto-only due to LQ   change!  %6.2f  > %4.2f\n",
                   i, Dets[i].StrName, Dets[i].DetName, PSA1.lq_lim[i] - PSA2.lq_lim[i], 2.0*5.0/2.355);
          // PZ tau
          if (fabs(PZI1.tau[i] - PZI2.tau[i]) > 1.0*0.4/2.355)           // ~ 1σ
            printf(" ---------- Detector  %2d %s %s veto-only due to tau  change!  %6.2f  > %4.2f\n",
                   i, Dets[i].StrName, Dets[i].DetName, PZI1.tau[i] - PZI2.tau[i], 1.0*0.4/2.355);
        }
      }
      printf("\ndt = %lf, tau = %lf μs, exp(dt/tau) = %lf\n", dt, PZI.tau[1], exp(dt/PZI.tau[1]));
      if (0) printf("s1, s2 = %13.10lf, %13.10lf\n", s1, s2);

      /* write data to PZ.output, ctc.output, and psa.output */
      PZ_info_write(&runInfo,  &PZI);
      CTC_info_write(&runInfo, &CTC);
      PSA_info_write(&runInfo, &PSA);

      /* write data to gains.output */
      if ((f_out = fopen("gains.output", "w"))) {
        printf("\n Writing energy calibrations to gains.output\n");
        for (i=0; i<runInfo.nGe; i++)
          fprintf(f_out,
                  "%3d %10.8lf %10.8lf %s\n",
                  i, Dets[i].HGcalib[0], Dets[i].LGcalib[0], Dets[i].StrName);
        fclose(f_out);
      }
      sprintf(line,
              "mv PZ.output %s/PZ.interp; mv ctc.output %s/ctc.interp; "
              "mv psa.output %s/psa.interp; mv gains.output %s/gains.interp",
              ds1, ds1, ds1, ds1);
      printf(">>> %s\n", line);
      system(line);
      printf("\n --------------------------------------- \n");
    }
  }

  printf("\n>>  All done... processed %d pairs of calibrations.\n\n", n_ds-1);
  return 0;
}

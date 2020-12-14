#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>

#include "MJDSort.h"

#define VERBOSE 0
unsigned int get_skim_time(char * fname);

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS], Dets1[NMJDETS], Dets2[NMJDETS];
  MJRunInfo  runInfo;
  CTCinfo    CTC, CTC1, CTC2;
  PSAinfo    PSA, PSA1, PSA2;
  PZinfo     PZI, PZI1, PZI2;
  int        i, j, nsd, veto, n_veto[200][10] = {{0}};
  int        ic, nc=0, order[100];
  char       calib_subset_dir[100][256]={{0}}, ds1[256] = "", ds2[256] = "", cwd[FILENAME_MAX];
  char       fn1[256], fn2[256], line[1024], line2[1024], veto_str[10];
  double     g1 = 0.5, g2=1.5, da, dd, dl, dq, ea, ed, el, eq;
  FILE       *f_in, *f_out, *f_out_veto;
  struct dirent *de;  // Pointer for directory entry
  DIR           *dp;
  float      limit_factor = 1.5;
  unsigned int startTime1, startTime2;

  int    take_mean = 0;
  double opt_bkgnd = -1.0;  // adjust mean to prioritize acceptance over background rejection
                            // (gives less uncertainty in exposure)


  if (argc > 1) {
    if (strstr(argv[1], "-m")) {
      take_mean = 1;
      printf("\n ---------------- Taking mean value of cuts ----------\n\n");
    } else if (strstr(argv[1], "-a")) {
      printf("\n ---------------- Prioritizing acceptance of cuts ----------\n\n");
    } else if (strstr(argv[1], "-b")) {
      opt_bkgnd = 1.0;
      printf("\n ---------------- Prioritizing background rejection of cuts ----------\n\n");
    } else {
      printf("\n Usage: %s [-m -a -b]\n"
             "        -m : Take mean value of the two A/E cuts\n"
             "        -b : Prioritize background rejection (default)\n", argv[0]);

      printf("\n ---------------- Using default: Prioritizing acceptance of cuts ----------\n\n");
    }
  }

  /* get current working directory */
  getcwd(cwd, FILENAME_MAX);

  /* get list of calibration directories */
  if ((dp = opendir(".")) == NULL) {
    printf("Could not open current directory" );
    return 0;
  }
  while ((de = readdir(dp)) != NULL) {
    if (strstr(de->d_name, "ds")) strncpy(calib_subset_dir[nc++], de->d_name, sizeof(calib_subset_dir[0]));
  }
  closedir(dp);

  /* sort list in order of run number */
  for (i=0; i<nc; i++) {
    order[i] = i;
    for (j=i-1; j>=0; j--) {
      if (strncmp(calib_subset_dir[i], calib_subset_dir[order[j]], 7) < 0) {
        order[j+1] = order[j];
      } else {
        order[j+1]= i;
        break;
      }
    }
    if (j < 0) order[0] = i;
  }
  printf("%d calibrations:\n", nc);
  for (i=0; i<nc; i++)
    printf("%3d %s\n", i, calib_subset_dir[order[i]]);

  /* open file for veto-only * detector list */
  f_out_veto = fopen("veto_only_runs.txt", "w");
  fprintf(f_out_veto, "###   cal_ds to cal_ds   det_ID   detector\n");

  /* --------- loop over consecutive pairs of calibration data subsets -------- */
  for (ic=0; ic<nc-1; ic++) {
    strncpy(ds1, calib_subset_dir[order[ic]], sizeof(ds1));
    strncpy(ds2, calib_subset_dir[order[ic+1]], sizeof(ds2));
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
    system(line2);  // restore any saved input files within this working directory

    /* get the start run toimes from the two skim files */
    sprintf(fn1, "%s/skim.dat", ds1);
    sprintf(fn2, "%s/skim.dat", ds2);
    startTime1 = get_skim_time(fn1);
    startTime2 = get_skim_time(fn2);
    if (!startTime1 || !startTime2) {
      printf("\nERROR: No start time for one of the calibrations!\n");
      return 0;
    }

    /* ---------------------- do interpolation ------------------*/
    memcpy(&PZI, &PZI1, sizeof(PZI));
    memcpy(&CTC, &CTC1, sizeof(CTC));
    memcpy(&PSA, &PSA1, sizeof(PSA));

    double dt = (PSA.e_ctc_rise[1] + PSA.e_ctc_flat[1]) / 200.0;    // μs
    double hg1, hg2, lg1, lg2;

    for (i=0; i<runInfo.nGe; i++) {
      PSA.ae_t0[i]      = 0;
      /* decide whether these data are valid (detector working okay) */
      if (!Dets[i].HGChEnabled) continue;
      if (CTC2.e_dt_slope[i] == 1.0 &&
          CTC2.e_lamda_slope[i] == 1.0 &&
          CTC2.e_lamda_gain[i] == 1.0) {          // Second calibration is invalid, so just take the first one
        Dets[i].HGcalib[0]    = Dets1[i].HGcalib[0];
        Dets[i].LGcalib[0]    = Dets1[i].LGcalib[0];
        printf("!!! Channel %d: Second calibration is invalid, so just taking the first one\n", i);
        continue;
      } else if (CTC1.e_dt_slope[i] == 1.0 &&
                 CTC1.e_lamda_slope[i] == 1.0 &&
                 CTC1.e_lamda_gain[i] == 1.0) {   // first calibration is invalid, so take the second one
        PZI.tau[i]            = PZI2.tau[i];
        PZI.baseline[i]       = PZI2.baseline[i];
        PZI.frac2[i]          = PZI2.frac2[i];
        PZI.bl_rms[i]         = PZI2.bl_rms[i];

        CTC.e_dt_slope[i]     = CTC2.e_dt_slope[i];
        CTC.e_lamda_slope[i]  = CTC2.e_lamda_slope[i];
        CTC.e_lamda_gain[i]   = CTC2.e_lamda_gain[i];
        CTC.best_dt_lamda[i]  = CTC2.best_dt_lamda[i];

        PSA.ae_dt_slope[i]    = PSA2.ae_dt_slope[i];
        PSA.ae_e_slope[i]     = PSA2.ae_e_slope[i];
        PSA.ae_pos[i]         = PSA2.ae_pos[i];
        PSA.ae_cut[i]         = PSA2.ae_cut[i];
        PSA.dcr_dt_slope[i]   = PSA2.dcr_dt_slope[i];
        PSA.dcr_lim[i]        = PSA2.dcr_lim[i];
        PSA.lamda_dt_slope[i] = PSA2.lamda_dt_slope[i];
        PSA.lamda_lim[i]      = PSA2.lamda_lim[i];
        PSA.lq_dt_slope[i]    = PSA2.lq_dt_slope[i];
        PSA.lq_lim[i]         = PSA2.lq_lim[i];
        Dets[i].HGcalib[0]    = Dets2[i].HGcalib[0];
        Dets[i].LGcalib[0]    = Dets2[i].LGcalib[0];
        printf("!!! Channel %d: First calibration is invalid, so just taking the second one\n", i);
        continue;
      }
          
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
      PSA.dcr_dt_slope[i]   = (PSA1.dcr_dt_slope[i]   + PSA2.dcr_dt_slope[i]  ) / 2.0;
      PSA.lamda_dt_slope[i] = (PSA1.lamda_dt_slope[i] + PSA2.lamda_dt_slope[i]) / 2.0;
      PSA.lq_dt_slope[i]    = (PSA1.lq_dt_slope[i]    + PSA2.lq_dt_slope[i]   ) / 2.0;

      da = PSA2.ae_cut[i]    - PSA1.ae_cut[i];
      dd = PSA2.dcr_lim[i]   - PSA1.dcr_lim[i];
      dl = PSA2.lamda_lim[i] - PSA1.lamda_lim[i];
      dq = PSA2.lq_lim[i]    - PSA1.lq_lim[i];

      ea = (PSA1.ae_pos[i] - PSA1.ae_cut[i] + PSA2.ae_pos[i] - PSA2.ae_cut[i]) / 6.0;  //  1/3 of mean A/E (pos-cut)
      ed = 8.0/2.355;                // typical σ for DCR peak
      el = 4.0/2.355;                // typical σ for lamda peak
      eq = 5.0/2.355;                // typical σ for LQ peak

      // calculate mean values of cuts
      PSA.ae_cut[i]    += da / 2.0;
      PSA.dcr_lim[i]   += dd / 2.0;
      PSA.lamda_lim[i] += dl / 2.0;
      PSA.lq_lim[i]    += dq / 2.0; // take the mean LQ cut value
      if (!take_mean) {
        // adjust mean to prioritize either acceptance or background rejection
        double factor = 0.5;  // fudge factor for how quickly the mean is adjusted; range should be 0.5 - 1.0
        PSA.ae_cut[i]    += opt_bkgnd * da * erf(factor * da / ea) / 2.0;
        PSA.dcr_lim[i]   -= opt_bkgnd * dd * erf(factor * dd / ed) / 2.0;
        PSA.lamda_lim[i] -= opt_bkgnd * dl * erf(factor * dl / el) / 2.0;
        PSA.lq_lim[i]    -= opt_bkgnd * dq * erf(factor * dq / eq) / 2.0;
      }
      //instead of mean or acceptance/bgnd, do linear time-interpolation of A/E cut value
      PSA.ae_cut[i]     = PSA1.ae_cut[i];
      PSA.ae_t0[i]      = startTime1;
      PSA.ae_t_slope[i] = 3600.0 * (PSA2.ae_cut[i] - PSA1.ae_cut[i]) / (startTime2 - startTime1); // change in cut value per hour

      hg1 = Dets1[i].HGcalib[0] * exp(dt/PZI.tau[i] - dt/PZI1.tau[i]);
      hg2 = Dets2[i].HGcalib[0] * exp(dt/PZI.tau[i] - dt/PZI2.tau[i]);
      lg1 = Dets1[i].LGcalib[0] * exp(dt/PZI.tau[i] - dt/PZI1.tau[i]);
      lg2 = Dets2[i].LGcalib[0] * exp(dt/PZI.tau[i] - dt/PZI2.tau[i]);
      Dets[i].HGcalib[0] = (hg1 + hg2) / 2.0;
      Dets[i].LGcalib[0] = (lg1 + lg2) / 2.0;
      Dets[i].HGcalib_unc[0] = fabs(hg1 - hg2) / 2.0;  // additional uncertainty in gain from gain or tau changes
      Dets[i].LGcalib_unc[0] = fabs(lg1 - lg2) / 2.0;

      /* ===================================================================================== */
      /* =========== check for big changes that imply veto-only status for ββ runs =========== */
      // gain
      sprintf(line, " ----- Detector  %2d %s %s veto-only at %s - %s (%d) due to",
              i, Dets[i].StrName, Dets[i].DetName, ds1, ds2, ic+1);
      veto = 0;
      if (Dets[i].DetName[0] == 'P' &&
          (1 || (i != 30 && i != 31 && i != 32))) {          // list of veto-only detectors for all runs; 1 to ignore
        // gain  - note multiplication by limit_factor
        // if (fabs(hg1 - hg2) > limit_factor * 0.0002) {                   // ~ 0.5 ppt ~ 1 keV at Qbb, assuming gain = 0.4
        if (fabs(hg1 - hg2) > limit_factor * 0.002 * Dets[i].HGcalib[0]) {  // 4 keV * limit_factor at Qbb
          printf("%s gain change!  %6.2f ppt\n", line, 1000.0 * (hg1 - hg2)/Dets[i].HGcalib[0]);
          n_veto[i][1]++;
          if (veto==1) n_veto[i][6]++;  // multiple reasons
          veto_str[veto++] = 'g';
        }
        // A/E  - note multiplication by limit_factor
        if (fabs(PSA1.ae_pos[i] - PSA2.ae_pos[i]) > limit_factor * 3.0*ea) {   //  = mean A/E (pos-cut) * limit_factor
          printf("%s A/E  change!  %6.2f, pos-cut = %5.2f\n", line,  PSA1.ae_pos[i] - PSA2.ae_pos[i], 3.0*ea);
          n_veto[i][2]++;
          if (veto==1) n_veto[i][6]++;  // multiple reasons
          veto_str[veto++] = 'A';
        }
        // DCR
        if (fabs(dd) > 3.0*ed) {              // ~ 3σ for DCR
          printf("%s DCR  change!  %6.2f  > %4.2f\n", line, PSA1.dcr_lim[i] - PSA2.dcr_lim[i], 2.0*8.0/2.355);
          n_veto[i][3]++;
          if (veto==1) n_veto[i][6]++;  // multiple reasons
          veto_str[veto++] = 'D';
        }
        // LQ
        if (fabs(dq) > 4.0*eq) {             // ~ 4σ for LQ
          printf("%s LQ   change!  %6.2f  > %4.2f\n", line, PSA1.lq_lim[i] - PSA2.lq_lim[i], 2.0*5.0/2.355);
          n_veto[i][4]++;
          if (veto==1) n_veto[i][6]++;  // multiple reasons
          veto_str[veto++] = 'Q';
        }
        // PZ tau
        if (fabs(PZI1.tau[i] - PZI2.tau[i]) > 2.0*0.4/2.355) {        // ~ 2σ for tau
          printf("%s tau  change!  %6.2f  > %4.2f\n", line, PZI1.tau[i] - PZI2.tau[i], 1.0*0.4/2.355);
          n_veto[i][5]++;
          if (veto==1) n_veto[i][6]++;  // multiple reasons
          veto_str[veto++] = 't';
        }
      }
      if (veto) {
        n_veto[i][0]++;
        veto_str[veto] = 0;
        fprintf(f_out_veto, "%3d %9s %9s   %3d  %s %s   %s\n",
                ic+1, ds1, ds2, i, Dets[i].StrName, Dets[i].DetName, veto_str);
      }
    }
    printf("\ndt = %lf, tau = %lf μs, exp(dt/tau) = %lf\n", dt, PZI.tau[1], exp(dt/PZI.tau[1]));

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

  printf("Veto-only datasets per detector:\n  ID  sets   gain A/E DCR  LQ tau   multiple\n");
  fprintf(f_out_veto, "#\n# ------------------------------\n"
          "# Veto-only datasets per detector:\n#   ID  sets   gain A/E DCR  LQ tau   multiple\n");
  for (i=0; i<100; i++) {
    if (n_veto[i][0]) {
      printf("%4d  %4d %6d %3d %3d %3d %3d %8d\n",
             i, n_veto[i][0], n_veto[i][1], n_veto[i][2],
             n_veto[i][3], n_veto[i][4], n_veto[i][5], n_veto[i][6]);
      fprintf(f_out_veto, "# %4d  %4d %6d %3d %3d %3d %3d %8d\n",
             i, n_veto[i][0], n_veto[i][1], n_veto[i][2],
             n_veto[i][3], n_veto[i][4], n_veto[i][5], n_veto[i][6]);
      for (j=0; j<7; j++) n_veto[100][j] += n_veto[i][j];
    }
  }
  printf("Total %4d %6d %3d %3d %3d %3d %8d\n",
         n_veto[100][0], n_veto[100][1], n_veto[100][2], n_veto[100][3],
         n_veto[100][4], n_veto[100][5], n_veto[100][6]);
  fprintf(f_out_veto, "# Total %4d %6d %3d %3d %3d %3d %8d\n",
         n_veto[100][0], n_veto[100][1], n_veto[100][2], n_veto[100][3],
         n_veto[100][4], n_veto[100][5], n_veto[100][6]);
    
  fclose(f_out_veto);  // close veto-only * detector list
  printf("\n>>  All done... processed %d pairs of calibrations.\n\n", nc-1);
  return 0;
}

/* ----------------------------------------------------------------------- */

unsigned int get_skim_time(char * fname) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        nsd;
  FILE      *f_in;


  if (!(f_in = fopen(fname,"r"))) {
    fprintf(stderr, "\n Failed to open input skim file %s\n", fname);
    return 0;
  }
  fread(&nsd, sizeof(int), 1, f_in);
  if (nsd == -2) fread(&nsd, sizeof(int), 1, f_in);
  fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
  fread(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_in);
  if (runInfo.idNum == 0) {
    printf("Flashcam data...\n");
    return 1;
  }
  printf(" Skim is data starting at %d seconds, run number %d from file %s\n",
         runInfo.startTime, runInfo.runNumber, runInfo.filename);

  return (unsigned int) runInfo.startTime;
}
  

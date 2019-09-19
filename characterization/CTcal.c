#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define MAKE_2D 1       // make a file with A/E, drift, and energy data for channel CHAN_2D, for 2d plots
#define CHAN_2D 0      // channel of interest for 2d plotting
#define HIS_COUNT 2200  // number of spectra in his[][] array
#define FWHM_RATIO 0.95 // ratio of FWHM to decide between E_dt and E_lamda charge-trapping correction

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1;


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname_in [chnum_lo] [chnum_hi] [e_lo] [e_hi] [-n]\n\n", argv[0]);
    return -1;
  }
  /* open skim data file as input */
  while (argn < argc && argv[argn][0] == '-') argn += 2;
  FILE *f_in = fopen(argv[argn],"r");
  if (f_in == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", argv[argn]);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", argv[argn]);
  runInfo.argc = argc;
  runInfo.argv = argv;

/* ---------------------------------------------------- */

  int clo=0, chi, elo=3000, ehi=7200;
  CTCinfo CTC;

  // skim data
  SavedData **sd;
  int      chan;
  float    drift, lamda;
  double   e_raw;
  int      nsd = 0, isd = 0;  // number of saved data, and pointer to saved data id

  double e_ctc, e_ctc_adc, e_lamda, e_lamda_adc, gain;
  float  pos, area, fwhm;
  int    i, j, k, n, roi_elo;
  int    *his[HIS_COUNT];
  FILE   *f_out, *f_out_2d = NULL;


  /* initialize */
  /* malloc and clear histogram space */
  if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
    printf("ERROR in CTcal.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1]+8192;

  // see if channel and energy limits are defined in the command line
  chi=100+runInfo.nGe-1;
  elo = 3000;
  ehi = 7200;
  if (runInfo.argc > 2) clo = atoi(runInfo.argv[2]);
  if (runInfo.argc > 3) chi = atoi(runInfo.argv[3]);
  if (runInfo.argc > 4) elo = atoi(runInfo.argv[4]);
  if (runInfo.argc > 5) ehi = atoi(runInfo.argv[5]);
  if (clo < 0) clo = 0;
  if (chi > 100+runInfo.nGe) chi = 100+runInfo.nGe;

  // read saved skim data from f_in
  fread(&nsd, sizeof(int), 1, f_in);
  fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
  fread(&runInfo, sizeof(runInfo), 1, f_in);
  /* malloc space for SavedData */
  if ((sd = malloc(nsd*sizeof(*sd))) == NULL ||
      (sd[0] = malloc(nsd*sizeof(SavedData))) == NULL) {
    printf("ERROR in CTcal.c; cannot malloc SavedData!\n");
    exit(-1);
  }
  for (i=1; i<nsd; i++) sd[i] = sd[i-1] + 1;
  fread(*sd, sizeof(SavedData), nsd, f_in);
  printf(" Skim is data from runs starting at number %d from file %s\n",
         runInfo.runNumber, runInfo.filename);

  /* read energy correction factors from ctc.input */
  if (!CTC_info_read(&runInfo, &CTC)) {
    printf("\n Warning: No inital charge-trapping correction data read. Does ctc.input exist?\n");
  }

  printf("\nChs %d to %d, e_trapmax %d to %d\n\n", clo, chi, elo, ehi);

  if (MAKE_2D &&
      ((CHAN_2D < 100 && Dets[CHAN_2D].HGChEnabled) ||
       (CHAN_2D >  99 && Dets[CHAN_2D%100].LGChEnabled))) {
    char fname[64];
    sprintf(fname, "CT_ch%3.3d_2d.dat", CHAN_2D);
    f_out_2d = fopen(fname, "w");
    fprintf(f_out_2d, "#chan    E_ctc    E_raw       DT   DT_corr    lamda lamda_corr  A/E\n");
  }

  // end of initialization
  // start loop over reading events from input file

  // ---------------------- steps 1-5 ----------------------
  for (int step = 1; step <= 5; step++) {
    printf(" ******************** Step %d of 5 ********************\n", step);

    /*
     * Step 1: Histogram: E_raw
     *            to get: E_raw gain
     * Step 2: Histogram: E_ctc vs. CTC (2614 peak), E_lamda vs. lamda
     *            to get: E_ctc slope/factor, E_lamda slope/factor
     * Step 3: Histogram: E_ctc_adc, E_lamda_adc
     *            to get: E_ctc gain, E_lamda gain
     * Step 4: Histogram: E_ctc, E_lamda
     */

    for (isd = 0; isd < nsd; isd++) {
      // if (isd%(nsd/10) == 0) printf(">>  event %7d (%d/10\n", isd, isd*10/nsd);
      chan   = sd[isd]->chan;
      e_raw  = sd[isd]->e;
      drift  = sd[isd]->drift;
      lamda  = sd[isd]->lamda;

      // histogram raw energy in [ADC] units
      if (step == 1) {
        if (e_raw < 8192) his[chan][(int)(e_raw + 0.5)]++;
        continue;
        // end of step-1 processing for this event
      }

      e_ctc_adc   = e_raw + drift * CTC.e_dt_slope[chan];

      /* ----------------------------------------------------------------------------------------
       * I tried subtracting the mean value of lamda, as a fuction of E, to correct residual INL's
       *   before using lamda to do a charge-trapping orrection to the energy.
       *   But we can't do this!  It artificially modifies the 2614 peak shape and FWHM in an
       *   unphysical way. It can shift the energy of an event by an amount proportional to the
       *   deriviative of the spectrum.
       *
       * So we can only do this mean-value subtraction to improve the resolution of lamda itself,
       *   for alpha rejection. And it must be done AFTER calculating e_lamda.
       *
       * if (subtract_dcr_mean && e_ctc_adc < 8000)   // correct for mean lamda vs. energy
       *           lamda -= (float) dcr_mean[chan][4000 + (int) e_ctc_adc/2] / 80.0;
       * ----------------------------------------------------------------------------------------  */

      lamda *= e_raw/6000.0;  // take into account energy dependence for charge-trapping correction
      e_lamda_adc = e_raw + lamda * CTC.e_lamda_slope[chan];

      if (chan < 100) {
        gain = Dets[chan].HGcalib[0];
      } else {
        gain = Dets[chan-100].LGcalib[0];
      }
      e_ctc   = e_ctc_adc   * gain;
      e_lamda = e_lamda_adc * gain * CTC.e_lamda_gain[chan];
      if (step == 2 && (k = (int) (e_raw * gain * 2.0 + 0.5)) > 0 && k < 8192)
          his[1800+chan][k]++; // raw energy, [0.5 keV]

      float dtc = gain * drift*2.0*2614.0/e_ctc;   // in (x)us, for DT- correction to A/E
      if (chan%100 == 50) dtc *= 3.0;              // special hack for detector 50; has especially strong variation

      // histogram energies in [ADC] and [0.5 keV] units
      if (step == 3 && e_ctc_adc < 8192 && e_lamda_adc < 8192) {
        his[600+chan][(int) (e_ctc_adc + 0.5)]++;
        his[800+chan][(int) (e_lamda_adc + 0.5)]++;
      }
      if (step == 4 && e_ctc < 4090 && e_lamda < 4090) {
        his[1000+chan][(int) (2.0*e_ctc + 0.5)]++;
        his[1200+chan][(int) (2.0*e_lamda + 0.5)]++;
        if (CTC.best_dt_lamda[chan]) {
          his[1400+chan][(int) (2.0*e_lamda + 0.5)]++;
        } else {
          his[1400+chan][(int) (2.0*e_ctc + 0.5)]++;
        }
        // test for effect of 0.5-keV binning on FWHM:
        if (e_ctc > 2500 && e_ctc < 2750) his[2000+chan][(int) (4.0*(e_ctc-1307.25) + 0.5)]++;
      }
      // make file for 2D plots  of E vs CTC
      roi_elo =  2574;
      if (f_out_2d && step == 4 && chan == CHAN_2D && e_ctc >= roi_elo && e_ctc <= roi_elo+700)
        fprintf(f_out_2d, "%4d %9.3f %8.3f  %8.3f %6.3f %10.3f %6.3f %10.2f\n",
                chan, e_ctc, e_raw*gain, drift, drift*CTC.e_dt_slope[chan]*gain,
                lamda, drift*CTC.e_lamda_slope[chan]*gain, sd[isd]->a_over_e);

      /* find optimum drift-time correction for energy */
      /* energy drift-time / trapping correction... */
      roi_elo = 2574;
      if (e_ctc >= roi_elo && e_ctc <= roi_elo+80) { // wide gate on 2614-keV peak
        if (step == 2) {
          /* try 40 options to find optimum charge-trapping correction for photopeak energy */
          for (j=0; j<40; j++) {
            float e = gain * (e_raw + drift * (j-5)/2.0);
            int  e2 = (0.5 + 2.0 * e) - 2*roi_elo;
            if (e2 > 0 && e2 < 200) his[200+chan][e2 + j*200]++;
            e = gain * (e_raw + lamda * (j-5)/2.0);
            e2 = (0.5 + 2.0 * e) - 2*roi_elo;
            if (e2 > 0 && e2 < 200) his[400+chan][e2 + j*200]++;
          }
        } else if (step == 4) {
          // histogram drift and lamda distributions for 2614 peak
          if (drift > -20 && drift < 20) {
            his[1600+chan][(int) (1000.5 + 10.0*drift*CTC.e_dt_slope[chan])]++;
            his[1600+chan][(int) (2000.5 + 10.0*drift)]++;
          }
          if (lamda > -20 && lamda < 20) {
            his[1600+chan][(int) (3000.5 + 10.0*lamda*CTC.e_lamda_slope[chan])]++;
            his[1600+chan][(int) (4000.5 + 10.0*lamda)]++;
          }
        }
      }
    }
    // --------------------------------------------------------------------------
    // end of event processing for this step
    // now process the current histograms to find calibrations, optimal correction factors, etc

    if (step == 1) {
      // extract E_raw gain
      /* find position of 2614.5-keV peak and use that to compute the calibration gains
         this code essentially taken from ep_finalize.c */
      for (i=0; i<=200; i++) {
        if (i%100 >= runInfo.nGe) continue;
        if (i%100 == 0) printf("\n");
        fwhm = 20;
        j = 5000;
        if (i > 99) {
          fwhm = 8;
          j = 1700;
        }
        if ((pos = autopeak3(his[i], j, 7000, &area, &fwhm))) {
          if (i < 100) {
            Dets[i].HGcalib[0] = 2614.5/pos;
            printf(" E_raw 2615: %3d P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   i, pos, area, fwhm*Dets[i].HGcalib[0]);
          } else {
            Dets[i-100].LGcalib[0] = 2614.5/pos;
            printf(" E_raw 2615: %3d P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   i, pos, area, fwhm*Dets[i-100].LGcalib[0]);
          }
        }
      }
      continue;
    }

    if (step == 2) {
      /* find narrowest charge-trapping-corrected peak in test spectra (2614.5-keV peak) */
      for (int spnum = 200; spnum <= 400; spnum += 200) {  // use his[200+chan} then his[400+chan]
        if (spnum == 200) printf("\n chan  old -> new optimum_E_DT_factor     cts    fwhm\n");
        if (spnum == 400) printf("\n chan  old -> new optimum_E_lamda_factor  cts    fwhm\n");

        for (chan=0; chan<200; chan++) {
          if (chan%100 >= runInfo.nGe) continue;
          j = k = n = 0;
          float f1 = 999;
          for (j=0; j<40; j++) {
            fwhm = 5;
            if ((pos = autopeak3(his[spnum+chan], 200*j, 200+200*j, &area, &fwhm)) &&
                area > 100 && fwhm < f1) {
              f1 = fwhm;
              k = area;
              n = j;
            }
          }
          if (k > 99) {
            if (spnum == 200) {
              printf("%3d %5.1f -> %5.1f  (%2d) %18d %9.2f\n",
                     chan, CTC.e_dt_slope[chan], (n-5)/2.0, n+1, k, f1);
              CTC.e_dt_slope[chan] = (n-5)/2.0;
            } else {
              printf("%3d %5.1f -> %5.1f  (%2d) %18d %9.2f\n",
                     chan, CTC.e_lamda_slope[chan], (n-5)/2.0, n+1, k, f1);
              CTC.e_lamda_slope[chan] = (n-5)/2.0;
            }
          }
        }
      }

    } else if (step == 3) {
      /* find position of 2614.5-keV peak and use that to compute the calibration gains
         this code essentially taken from ep_finalize.c */
      FILE *fp = fopen("fwhm_ctc.txt", "w");
      fprintf(fp, "#chan   DTC    lamdaC\n");
      float fwhm0 = 0;
      for (chan=0; chan<200; chan++) {
        CTC.best_dt_lamda[chan] = 0;
        if (chan%100 >= runInfo.nGe) continue;
        if (chan%100 == 0) printf("\n");
        // first for E_ctc
        fwhm = 5;
        j = 5000;
        if (chan > 99) {
          fwhm = 3;
          j = 1700;
        }
        if ((pos = autopeak3(his[600+chan], j, 7000, &area, &fwhm))) {
          fwhm0 = fwhm;
          if (chan < 100) {
            Dets[chan].HGcalib[0] = 2614.5/pos;
            printf("Ch %3d  E_ctc_adc   2615:  P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   chan, pos, area, fwhm*Dets[chan].HGcalib[0]);
            fprintf(fp, "%4d %7.3f", chan, fwhm*Dets[chan].HGcalib[0]);
          } else {
            Dets[chan-100].LGcalib[0] = 2614.5/pos;
            printf("Ch %3d  E_ctc_adc   2615:  P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   chan, pos, area, fwhm*Dets[chan-100].LGcalib[0]);
            fprintf(fp, "%4d %7.3f", chan, fwhm*Dets[chan-100].LGcalib[0]);
          }
        }

        // now for E_lamda
        fwhm = 5;
        j = 5000;
        if (chan > 99) {
          fwhm = 3;
          j = 1700;
        }
        if ((pos = autopeak3(his[800+chan], j, 7000, &area, &fwhm))) {
          if (fwhm < fwhm0 * FWHM_RATIO) CTC.best_dt_lamda[chan] = 1;
          if (chan < 100) {
            CTC.e_lamda_gain[chan] = 2614.5/pos/Dets[chan].HGcalib[0];
            printf("Ch %3d  E_lamda_adc 2615:  P = %7.1f A = %7.0f; FWHM = %7.2f keV ; rel gain = %9.7lf\n",
                   chan, pos, area, fwhm*Dets[chan].HGcalib[0]*CTC.e_lamda_gain[chan], CTC.e_lamda_gain[chan]);
            fprintf(fp, " %7.3f\n", fwhm*Dets[chan].HGcalib[0]*CTC.e_lamda_gain[chan]);
          } else {
            CTC.e_lamda_gain[chan] = 2614.5/pos/Dets[chan-100].LGcalib[0];
            printf("Ch %3d  E_lamda_adc 2615:  P = %7.1f A = %7.0f; FWHM = %7.2f keV ; rel gain = %9.7lf\n",
                   chan, pos, area, fwhm*Dets[chan-100].LGcalib[0]*CTC.e_lamda_gain[chan], CTC.e_lamda_gain[chan]);
            fprintf(fp, " %7.3f\n", fwhm*Dets[chan-100].LGcalib[0]*CTC.e_lamda_gain[chan]);
          }
        }
      }
      fclose(fp);

    }
  }

  printf(">>  All done...\n"); fflush(stdout);

  // write energy calibrations to gains.output
  if ((f_out = fopen("gains.output","w"))) {
    printf("\n Writing energy calibrations to gains.output\n");
    for (i=0; i<runInfo.nGe; i++)
      fprintf(f_out,
              "%3d %10.8lf %10.8lf %s\n",
              i, Dets[i].HGcalib[0], Dets[i].LGcalib[0], Dets[i].StrName);
    fclose(f_out);
  }

  // write charge-trapping data to ctc.output
  CTC_info_write(&runInfo, &CTC);

  // write out histograms
  f_out = fopen("ctc.rms", "w");
  for (i=0; i<HIS_COUNT; i++) {
   char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d raw uncorrected energy [ADC]; run %d", i, i, runInfo.runNumber);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d Energy DT-corrected with test factors", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d Energy lamda-corrected with test factors", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d DT-corrected energy [ADC]", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d lamda-corrected energy [ADC]", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d DT-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d lamda-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d optimally corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; ch %d drift time, misc", i, i%200);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d raw (non-DT-corrected) energy [0.5 keV]", i, i%200);
    } else {
      sprintf(spname, "%d;", i);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  return 0;
}

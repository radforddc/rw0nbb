#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define MAKE_2D 0       // make a file with A/E, drift, and energy data for channel CHAN_2D, for 2d plots
#define CHAN_2D 0       // channel of interest for 2d plotting
#define HIS_COUNT 3000  // number of spectra in his[][] array
#define FWHM_RATIO 0.95 // ratio of FWHM to decide between E_dt and E_lamda charge-trapping correction

int CTC2_info_write(MJRunInfo *runInfo, CTC2info *CTC2);

/* -------------------------------------------------------------------

ctc.rms:
Sp ID        n:   detID n,  raw uncorrected energy [ADC]
Sp ID  200 + n:   detID n,  Energy DT-corrected with test factors
Sp ID  400 + n:   detID n,  Energy lamda-corrected with test factors
Sp ID  600 + n:   detID n,  DT-corrected energy [ADC]
Sp ID  800 + n:   detID n,  lamda-corrected energy [ADC]
Sp ID 1000 + n:   detID n,  DT-corrected energy [0.5 keV]
Sp ID 1200 + n:   detID n,  lamda-corrected energy [0.5 keV]
Sp ID 1400 + n:   detID n,  optimally corrected energy [0.5 keV]
Sp ID 1600 + n:   detID n,  drift time, misc
Sp ID 1800 + n:   detID n,  raw (non-DT-corrected) energy [0.5 keV]
Sp ID 2000 + n:   detID n,  optimally corrected energy [0.25 keV]

ctc2.rms:
Sp ID        n:   detID n,  Energy DT-corrected with test factors    (by CTcal) [0.5 keV]
Sp ID  200 + n:   detID n,  Energy lamda-corrected with test factors (by CTcal) [0.5 keV]
Sp ID  400 + n:   detID n,  Mean energy (DT- and lamda-corrected)    (by CTcal) [0.5 keV]
Sp ID  600 + n:   detID n,  New Quadratic-DT-corrected energy [ADC]
Sp ID  800 + n:   detID n,  New DT-corrected with linear fit          [0.5 keV]
Sp ID 1000 + n:   detID n,  New lamda-corrected with linear fit       [0.5 keV]
Sp ID 1200 + n:   detID n,  New Mean energy (DT- and lamda-corrected) [0.5 keV]
Sp ID 1400 + n:   detID n,  New Quadratic-DT-corrected energy         [0.5 keV]
Sp ID 1600 + n:   detID n,  New Quadratic Mean energy (QDT and lamda) [0.5 keV]
Sp ID 1800 + n:   detID n,  New DT-corrected with linear fit          [0.25 keV]
Sp ID 2000 + n:   detID n,  New lamda-corrected with linear fit       [0.25 keV]
Sp ID 2200 + n:   detID n,  New Mean energy (DT- and lamda-corrected) [0.25 keV]
Sp ID 2400 + n:   detID n,  New Quadratic-DT-corrected energy         [0.25 keV]
Sp ID 2600 + n:   detID n,  New Quadratic Mean energy (QDT and lamda) [0.25 keV]
Sp ID 2800 + n:   detID n,  Best option for corrected energy          [0.25 keV]

The detID is  0-57 for high gain channels, 100-157 for low gain.

   ------------------------------------------------------------------- */


int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1, make_2d = MAKE_2D, chan_2d = CHAN_2D;


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [-l 2d_chan_num] fname_in [chnum_lo] [chnum_hi] [e_lo] [e_hi]\n\n", argv[0]);
    return -1;
  }
  /* open skim data file as input */
  while (argn < argc && argv[argn][0] == '-') {
    if (argv[argn][1] == 'l') {
      make_2d = 1;
      chan_2d = atoi(argv[argn+1]);
    }
    argn += 2;
  }
  FILE *f_in = fopen(argv[argn],"r");
  if (f_in == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", argv[argn]);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", argv[argn]);
  runInfo.argc = argc;
  runInfo.argv = argv;

/* ---------------------------------------------------- */

  int clo=0, chi=-1, elo=1, ehi=16000;
  CTCinfo CTC;
  CTC2info CTC2;

  // skim data
  SavedData  **sd1;
  SavedData2 **sd2;
  int      sd_version = 1;
  int      chan;
  float    drift, lamda, aovere, dcr;
  double   e_raw;
  int      nsd = 0, isd = 0;  // number of saved data, and pointer to saved data id

  double e_ctc, e_ctc_adc, e_lamda, e_lamda_adc, gain;
  float  pos, area, fwhm;
  int    i, j, roi_elo;
  int    *his[HIS_COUNT];
  FILE   *fp, *f_out, *f_out_2d = NULL;

  double s00[2][200]={{0}}, s01[2][200]={{0}}, s10[2][200]={{0}}, s11[2][200]={{0}};
  double s20[2][200]={{0}}, s21[2][200]={{0}}, s30[2][200]={{0}}, s40[2][200]={{0}};
  float  e_qdt_slope[200];      // linear factor for quadratic drift-time correction of energy
  float  e_qdt_quad[200]={0};   // quad factor for quadratic drift-time correction of energy
  double e_qdt_gain[200];       // gain for quadratic drift-time correction of energy
  double e_qdtc, e_qdtc_adc;


  /* initialize */
  /* malloc and clear histogram space */
  if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
    printf("ERROR in CTcal.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1]+8192;
  strncpy(CTC2.ctc2_fname, "ctc2.input", sizeof(CTC2.ctc2_fname));
  for (i=1; i<200; i++) {
    CTC2.e_qdt_slope[i] = 1;     // linear factor for quadratic drift-time correction of energy
    CTC2.e_qdt_quad[i] = 0;      // quadratic factor for quadratic drift-time correction of energy
    CTC2.e_qdt_gain[i] = 1;      // gain for quadratic drift-time-corrected energy
    CTC2.e_dt_slope[i] = 1;      // new factor for mean (qdt + lamda)/2 correction of energy
    CTC2.e_dt_gain[i] = 1;       // new gain for linear dt corrected energy
    CTC2.e_lamda_slope[i] = 1;   // new factor for lamda correction of energy
    CTC2.e_lamda_gain[i] = 1;    // new gain for lamda-corrected energy
    CTC2.best_ctc2_res[i] = 0;   // 0 = dt, 1 = lamda, 2 = linear mean, 3 = qdt, 4 = quadratic mean
  }

  // see if channel and energy limits are defined in the command line
  // chi=100+runInfo.nGe-1;   // runInfo.nGe not yet set! See later.
  argn++;
  if (runInfo.argc > argn) clo = atoi(runInfo.argv[argn++]);
  if (runInfo.argc > argn) chi = atoi(runInfo.argv[argn++]);
  if (runInfo.argc > argn) elo = atoi(runInfo.argv[argn++]);
  if (runInfo.argc > argn) ehi = atoi(runInfo.argv[argn++]);
  if (clo < 0) clo = 0;

  // read saved skim data from f_in
  fread(&nsd, sizeof(int), 1, f_in);
  if (nsd == -2) {
    sd_version = 2;
    fread(&nsd, sizeof(int), 1, f_in);
  }
  fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
  fread(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_in);
  if (runInfo.idNum == 0) {
    runInfo.flashcam = 1;
    fread(&(runInfo.flashcam), 8*sizeof(int), 1, f_in);
  }
  /* malloc space for SavedData */
  if ((sd_version == 1 &&
       ((sd1 = malloc(nsd*sizeof(*sd1))) == NULL ||
        (sd1[0] = malloc(nsd*sizeof(*sd1[0]))) == NULL)) ||
      (sd_version == 2 &&
       ((sd2 = malloc(nsd*sizeof(*sd2))) == NULL ||
        (sd2[0] = malloc(nsd*sizeof(*sd2[0]))) == NULL))) {
    printf("ERROR in CTcal.c; cannot malloc SavedData!\n");
    exit(-1);
  }
  printf("Skim data mode = %d;  %d detectors, %d skimmed events\n",
         sd_version, runInfo.nGe, nsd);
  if (sd_version == 1) {
    for (i=1; i<nsd; i++) sd1[i] = sd1[i-1] + 1;
    fread(*sd1, sizeof(**sd1), nsd, f_in);
  } else {   // sd_version == 2
    for (i=1; i<nsd; i++) sd2[i] = sd2[i-1] + 1;
    fread(*sd2, sizeof(**sd2), nsd, f_in);
  }
  printf(" Skim is data from runs starting at number %d from file %s\n",
         runInfo.runNumber, runInfo.filename);

  /* read gains from gains.input (will have changed since skim was created) */
  if ((fp = fopen(ECAL_FILENAME,"r"))) {
    printf("\n Reading energy calibrations from %s\n", ECAL_FILENAME);
    char line[256];
    double g1 = 0.5, g2=1.5;
    while (fgets(line, sizeof(line), fp)) {
      if (*line != '#' &&
          sscanf(line, "%d %lf %lf", &i, &g1, &g2) == 3 &&
          i >=0 && i < runInfo.nGe) {
        Dets[i].HGcalib[0] = g1;
        Dets[i].LGcalib[0] = g2;
      }
    }
    fclose(fp);
  } else {
    printf("ERROR: Cannot open %s !\n", ECAL_FILENAME);
    exit(-1);
  }

  /* read energy correction factors from ctc.input */
  if (!CTC_info_read(&runInfo, &CTC)) {
    printf("\n Warning: No initial charge-trapping correction data read. Does ctc.input exist?\n");
  }

  if (chi == -1) chi=100+runInfo.nGe-1;   // runInfo.nGe is now set!
  printf("\nChs %d to %d, e_trapmax %d to %d\n\n", clo, chi, elo, ehi);

  if (make_2d &&
      ((chan_2d < 100 && Dets[chan_2d].HGChEnabled) ||
       (chan_2d >  99 && Dets[chan_2d%100].LGChEnabled))) {
    char fname[64];
    sprintf(fname, "CT_ch%3.3d_2d.dat", chan_2d);
    f_out_2d = fopen(fname, "w");
    fprintf(f_out_2d, "#chan    E_raw    E_ctc   E_qdtc  E_lamda        DT    lamda    DCR       A/E\n");
  }

  // end of initialization
  // start loop over reading events from input file

  // ---------------------- steps 5-8 ---------------------- continues from CTcal
  for (int step = 5; step <= 8; step++) {
    printf(" ******************** Step %d of 8 ********************\n", step);

    /*
     * Step 5: Histogram: E_ctc, E_lamda
     *         Quadratic  fit of  E_raw vs. DT     to get e_qdt_slope and _quad;
     *         New linear fit of  E_raw vs. DT     to get e_dt_slope;
     *         New linear fit of  E_raw vs. lamda  to get e_lamda_slope;
     * Step 6: Histogram: E_qdtc_adc
     *            to get: E_qdtc gain
     * Step 7: Histogram: new E_ctc, E_lamda, E_mean, E_qdtc, E_qmean
     *            to get: FWHM
     * Step 7: Histogram: new best option for corrected energy, at 0.25 keV/bin
     *            to get: FWHM, summed FWHM, etc
     */

    for (isd = 0; isd < nsd; isd++) {
      // if (isd%(nsd/10) == 0) printf(">>  event %7d (%d/10\n", isd, isd*10/nsd);
      if (sd_version == 1) {
        chan   = sd1[isd]->chan;
        e_raw  = sd1[isd]->e;
        drift  = sd1[isd]->drift;
        aovere = sd1[isd]->a_over_e;
        lamda  = sd1[isd]->lamda;
        dcr    = sd1[isd]->dcr;
      } else {                       // sd_version == 2
        chan   = sd2[isd]->chan;
        e_raw  = sd2[isd]->e;
        drift  = sd2[isd]->drift;
        aovere = sd2[isd]->a_over_e;
        lamda  = sd2[isd]->lamda;
        dcr    = sd2[isd]->dcr;
      }
      if (chan < clo || chan > chi) continue;
      if (chan < 100 && (e_raw < elo || e_raw > ehi)) continue;
      if (chan > 99 && (e_raw < elo/3.4 || e_raw > ehi/3.2)) continue;

      lamda *= e_raw/6000.0;  // take into account energy dependence for charge-trapping correction
      e_ctc_adc   = e_raw + drift * CTC.e_dt_slope[chan];
      e_lamda_adc = e_raw + lamda * CTC.e_lamda_slope[chan];
      if (chan < 100) {
        gain = Dets[chan].HGcalib[0];
      } else {
        gain = Dets[chan-100].LGcalib[0];
      }
      e_ctc   = e_ctc_adc   * gain;
      e_lamda = e_lamda_adc * gain * CTC.e_lamda_gain[chan];
      if (e_ctc > 4000 || e_lamda > 4000) continue;

      if (step == 5) {
        // histogram energies corrected in earlier CTcal
        his[    chan][(int) (2.0*e_ctc + 0.5)]++;
        his[200+chan][(int) (2.0*e_lamda + 0.5)]++;
        his[400+chan][(int) (e_ctc + e_lamda + 0.5)]++;
        if (e_ctc > 2611.5 && e_ctc < 2617.5) {       //   +/- 3 keV energy window
          // evaluate sums for linear and quadratic fits of e_raw vs. drift time
          s00[0][chan]++;
          s10[0][chan] += drift;
          s20[0][chan] += drift * drift;
          s30[0][chan] += drift * drift * drift;
          s40[0][chan] += drift * drift * drift * drift;
          s01[0][chan] += e_raw;
          s11[0][chan] += e_raw * drift;
          s21[0][chan] += e_raw * drift * drift;
          // evaluate sums for linear fit of e_raw vs. lamda
          // if (lamda < -2 || lamda > 4) continue;
          s00[1][chan]++;
          s10[1][chan] += lamda;
          s20[1][chan] += lamda * lamda;
          s01[1][chan] += e_raw ;
          s11[1][chan] += e_raw * lamda;
        }
        continue;       // end of step-5 processing for this event
      }

      e_qdtc_adc = e_raw + (e_qdt_slope[chan] + e_qdt_quad[chan]*drift) * drift;

      // histogram energies in [ADC] and [0.5 keV] units
      if (step == 6) {
        // histogram energies corrected with new quadratic fit
        if (e_qdtc_adc < 8192) his[600+chan][(int) (e_qdtc_adc + 0.5)]++;
      } else if (step == 7 && e_ctc < 4000) {
        e_qdtc = e_qdtc_adc * e_qdt_gain[chan];
        // histogram energies corrected with new linear fit
        his[800+chan][(int) (2.0*e_ctc + 0.5)]++;
        his[1000+chan][(int) (2.0*e_lamda + 0.5)]++;
        his[1200+chan][(int) (e_ctc + e_lamda + 0.5)]++;
        // histogram energies corrected with new quadratic fit and mean
        his[1400+chan][(int) (2.0*e_qdtc + 0.5)]++;
        his[1600+chan][(int) (e_qdtc + e_lamda + 0.5)]++;

        // now at 0.25 keV/bin
        if (e_ctc > 1310 && e_ctc < 2750) {
          his[1800+chan][(int) (4.0*(e_ctc - 1307.25) + 0.5)]++;
          his[2000+chan][(int) (4.0*(e_lamda - 1307.25) + 0.5)]++;
          his[2200+chan][(int) (2.0*(e_ctc + e_lamda - 2.0*1307.25) + 0.5)]++;
          his[2400+chan][(int) (4.0*(e_qdtc - 1307.25) + 0.5)]++;
          his[2600+chan][(int) (2.0*(e_qdtc + e_lamda - 2.0*1307.25) + 0.5)]++;
        }

      } else if (step == 8 && e_ctc < 2750 && e_ctc > 1320) {
        double *sgain = Dets[chan].HGcalib;
        if  (chan >= 100) sgain = Dets[chan-100].LGcalib;
        float e0, e1, e3, ebest;
        e0 = e_raw + drift * CTC2.e_dt_slope[chan];
        e1 = e_raw + lamda * CTC2.e_lamda_slope[chan];
        e3 = e_raw + (CTC2.e_qdt_slope[chan] + CTC2.e_qdt_quad[chan]*drift) * drift;
        int best = CTC2.best_ctc2_res[chan];
        if      (best == 0) ebest = e0 * sgain[0];
        else if (best == 1) ebest = e1 * sgain[1];
        else if (best == 2) ebest = (e0 + e1)/2.0 * sgain[2];
        else if (best == 3) ebest = e3 * sgain[3];
        else if (best == 4) ebest = (e3 + e1)/2.0 * sgain[4];
        else {
          printf("ERROR! best_res option for chan %d = %d\n\n", chan, best);
          exit(-1);
        }
        his[2800+chan][(int) (4.0*(ebest - 1307.25) + 0.5)]++;
      }

      if (step == 8) {
        // make file for 2D plots  of E vs CTC
        // roi_elo = DEP_E - 40.0;
        // if (f_out_2d && step == 8 && chan == chan_2d && e_ctc >= roi_elo && e_ctc <= roi_elo+80)
        roi_elo = CAL_E - 40.0;
        if (f_out_2d && step == 8 && chan == chan_2d && e_ctc >= roi_elo && e_ctc <= roi_elo+700) {
          fprintf(f_out_2d, "%4d %9.2f %8.2f %8.2f %8.2f %10.3f %7.3f %7.3f %9.2f\n",
                  chan, e_raw*gain, e_ctc, e_qdtc, e_lamda, drift, lamda, dcr, aovere);
        }
      }
    }

    // --------------------------------------------------------------------------
    // end of event processing for this step
    // now process the current histograms to find calibrations, optimal correction factors, etc

    if (step == 5) {
      /* analyse fit sums to extract quadratice fits of energies vs drift time and lamda */
      for (chan=0; chan<200; chan++) {
        e_qdt_slope[chan] = CTC.e_dt_slope[chan];
        if (s00[0][chan] > 500) {
          // quadratic fit against drift time:
          double d = (s00[0][chan]*s20[0][chan]*s40[0][chan] - s10[0][chan]*s10[0][chan]*s40[0][chan] -
                      s00[0][chan]*s30[0][chan]*s30[0][chan] + s10[0][chan]*s20[0][chan]*s30[0][chan]*2.0 -
                      s20[0][chan]*s20[0][chan]*s20[0][chan]);
          if (fabs(d) < 0.1) continue;
          e_qdt_slope[chan] = -(s11[0][chan]*s00[0][chan]*s40[0][chan] - s01[0][chan]*s10[0][chan]*s40[0][chan] +
                                s01[0][chan]*s20[0][chan]*s30[0][chan] - s21[0][chan]*s00[0][chan]*s30[0][chan] -
                                s11[0][chan]*s20[0][chan]*s20[0][chan] + s21[0][chan]*s10[0][chan]*s20[0][chan]) / d;
          e_qdt_quad[chan]  = -(s01[0][chan]*s10[0][chan]*s30[0][chan] - s11[0][chan]*s00[0][chan]*s30[0][chan] -
                                s01[0][chan]*s20[0][chan]*s20[0][chan] + s11[0][chan]*s10[0][chan]*s20[0][chan] +
                                s21[0][chan]*s00[0][chan]*s20[0][chan] - s21[0][chan]*s10[0][chan]*s10[0][chan]) / d;
          CTC2.e_qdt_slope[chan] = e_qdt_slope[chan];
          CTC2.e_qdt_quad[chan]  = e_qdt_quad[chan];
          // linear fit:
          double a = -(s10[0][chan]*s01[0][chan] - s11[0][chan]*s00[0][chan])/
                      (s10[0][chan]*s10[0][chan] - s20[0][chan]*s00[0][chan]);
          printf("%3d %5.1f -> %5.1f  -> %5.1f %6.2f;  ",
                 chan, CTC.e_dt_slope[chan], a, e_qdt_slope[chan], e_qdt_quad[chan]);
          CTC.e_dt_slope[chan] = CTC2.e_dt_slope[chan] = a;

          // linear fit against lamda:
          a = -(s10[1][chan]*s01[1][chan] - s11[1][chan]*s00[1][chan])/
               (s10[1][chan]*s10[1][chan] - s20[1][chan]*s00[1][chan]);
          printf("%5.1f -> %5.1f\n", CTC.e_lamda_slope[chan], a);
          CTC.e_lamda_slope[chan] = CTC2.e_lamda_slope[chan] = a;
          if (0 && (chan == 110 || chan == 10 || chan == 112))
            printf("s10 s01 s11 s00 = %.2e %.2e %.2e %.2e\n"
                   "s10 s10 s20 s00 = %.2e %.2e %.2e %.2e   %.2e\n",
                   s10[1][chan], s01[1][chan], s11[1][chan], s00[1][chan],
                   s10[1][chan], s10[1][chan], s20[1][chan], s00[1][chan], s00[0][chan]);
        }
      }

    } else if (step == 6) {
      /* find position of 2614.5-keV peak and use that to compute the qdt calibration gain */
      printf(" Finding gain for Quad DT correction\n");
      for (chan=0; chan<200; chan++) {
        e_qdt_gain[chan] = 1;
        if (s00[0][chan] < 501) continue;  // not enough events
        if (chan%100 >= runInfo.nGe) continue;
        if (chan%100 == 0) printf("\n");
        // for E_qdtc
        fwhm = 5;
        j = 4000;
        if (runInfo.flashcam == 3) j = 3000;  // HADES data has lower effective gain
        if (chan > 99) {
          fwhm = 3;
          j = 1700;
        }
        if ((pos = autopeak3(his[600+chan], j, 8000, &area, &fwhm))) {
          e_qdt_gain[chan] = CAL_E/pos;
          printf("Ch %3d %6.1f:   E_qdtc_adc  P = %7.1f A = %7.0f  FWHM = %7.2f keV\n",
                 chan, CAL_E, pos, area, fwhm*e_qdt_gain[chan]);
        }
      }

    } else if (step == 7) {
      /* find widths of 2614.5-keV peak for each correction option and use that to determine best option */
      printf(" Finding peak widths for each correction option\n");
      printf("    ...Results will be written to fwhm_ctc2.txt\n");
      fp = fopen("fwhm_ctc2.txt", "w");
      fprintf(fp, "#chan DT_lin  lamda mean_lin DT_quad Mean_quad Best BestOption");
      float penalty[5] = {1.0, 0.95, 0.975, 0.998, 0.975};   // penalty factor for FWHM from each option
                                                             //    (FWHM_RATIO in CTcal.c)
      for (chan=0; chan<200; chan++) {
        double gain[5] = {1};
        float fwhm0[5] = {99}, area0[5] = {100};
        int   best = 0;
        if (chan%100 == 99) fprintf(fp, "\n");
        if (s00[0][chan] < 501) continue;  // not enough events
        if (chan%100 >= runInfo.nGe) continue;
        fprintf(fp, "\n%4d ", chan);
        for (int iopt=0; iopt < 5; iopt++) {
          fwhm = 6;
          j = 2 * CAL_E - 200;
          // if ((pos = autopeak3(his[200*iopt + 800 + chan], j, j+400, &area, &fwhm))) {  // equiv to autopeak4(,,,2.0f,,)
          if ((pos = autopeak4(his[200*iopt + 800 + chan], j, j+400, 1.5f, &area0[iopt], &fwhm))) {
            // printf("Ch %3d option %d  P = %7.1f A = %7.0f  FWHM = %7.2f keV\n",
            //        chan, iopt, pos, area, fwhm/2.0);
            fprintf(fp, "%7.3f", fwhm/2.0);
            if (fwhm0[best] * penalty[iopt] > fwhm * penalty[best]) best = iopt;
            fwhm0[iopt] = fwhm;
            gain[iopt] = 2.0 * CAL_E/pos;
            // Dets[chan].HGcalib[2] = CAL_E/pos;
          } else {
            fprintf(fp, "  -----");            
            area0[iopt] = 1;
          }
        }
        fprintf(fp, "%12.3f   %d", fwhm0[best]/2.0, best);
        // save new gain values
        double *sgain = Dets[chan].HGcalib, *sgunc = Dets[chan].HGcalib_unc;
        if  (chan >= 100) {
          sgain = Dets[chan-100].LGcalib;
          sgunc = Dets[chan-100].LGcalib_unc;
        }
        gain[0] *= sgain[0];
        gain[1] *= sgain[0] * CTC.e_lamda_gain[chan];
        gain[2] *= sgain[0] * (1.0+CTC.e_lamda_gain[chan])/2.0;
        gain[3] *= e_qdt_gain[chan];
        gain[4] *= e_qdt_gain[chan] * (1.0+CTC.e_lamda_gain[chan])/2.0;
        for (i=0; i<5; i++) {
          sgain[i] = gain[i];
          sgunc[i] = gain[i] * (1.0 / 2.355) * fwhm0[i] / sqrt(area0[i]) / CAL_E;
        }
        CTC2.e_dt_gain[chan] = gain[0];
        CTC2.e_lamda_gain[chan] = gain[1];
        CTC2.e_qdt_gain[chan] = gain[3];
        CTC2.best_ctc2_res[chan] = best;
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
  }

  printf(">>  All done...\n"); fflush(stdout);

  // write energy calibrations to gains.output
  if ((f_out = fopen("gains2.output","w"))) {
    printf("\n Writing energy calibrations to gains.output\n");
    fprintf(f_out, "chan detector"
            " HG: lin_dt    lamda   lin_mean   quad_dt quad_mean  "
            //  9.7lfxxxx 9.7lfxxxx 9.7lfxxxx 9.7lfxxxx 9.7lfxxxx
            " LG: lin_dt    lamda   lin_mean   quad_dt quad_mean\n           "
            "    ............... HG uncertainties ................"
            "    ............... LG uncertainties ................\n");
            //        9.7lfxxxx 9.7lfxxxx 9.7lfxxxx 9.7lfxxxx 9.7lfxxxx
    for (i=0; i<runInfo.nGe; i++) 
      fprintf(f_out,
              "%3d  %6s    %9.7lf %9.7lf %9.7lf %9.7lf %9.7lf    %9.7lf %9.7lf %9.7lf %9.7lf %9.7lf\n"
              "               %9.7lf %9.7lf %9.7lf %9.7lf %9.7lf    %9.7lf %9.7lf %9.7lf %9.7lf %9.7lf\n",
              i, Dets[i].StrName,
              Dets[i].HGcalib[0], Dets[i].HGcalib[1], Dets[i].HGcalib[2], Dets[i].HGcalib[3], Dets[i].HGcalib[4],
              Dets[i].LGcalib[0], Dets[i].LGcalib[1], Dets[i].LGcalib[2], Dets[i].LGcalib[3], Dets[i].LGcalib[4],
              Dets[i].HGcalib_unc[0], Dets[i].HGcalib_unc[1], Dets[i].HGcalib_unc[2], Dets[i].HGcalib_unc[3],
              Dets[i].HGcalib_unc[4], Dets[i].LGcalib_unc[0], Dets[i].LGcalib_unc[1], Dets[i].LGcalib_unc[2],
              Dets[i].LGcalib_unc[3], Dets[i].LGcalib_unc[4]);
    fclose(f_out);
  }

  // write charge-trapping data to ctc.output
  CTC2_info_write(&runInfo, &CTC2);

  // write out histograms
  f_out = fopen("ctc2.rms", "w");
  for (i=0; i<HIS_COUNT; i++) {
    char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d old linear DT-corrected energy [0.5 keV]", i, i);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d old lamda-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d old mean linear corrected energy [0.5 keV]", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d New QDT-corrected energy [ADC]", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d New DT-corrected with linear fit          [0.5 keV]", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d New lamda-corrected with linear fit       [0.5 keV]", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d New Mean energy (DT- and lamda-corrected) [0.5 keV]", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d New Quadratic-DT-corrected energy         [0.5 keV]", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; ch %d New Quadratic Mean energy (QDT and lamda) [0.5 keV]", i, i%200);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d New DT-corrected with linear fit          [0.25 keV]", i, i%200);
    } else if (i < 2200) {
      sprintf(spname, "%d; ch %d New lamda-corrected with linear fit       [0.25 keV]", i, i%200);
    } else if (i < 2400) {
      sprintf(spname, "%d; ch %d New Mean energy (DT- and lamda-corrected) [0.25 keV]", i, i%200);
    } else if (i < 2600) {
      sprintf(spname, "%d; ch %d New Quadratic-DT-corrected energy         [0.25 keV]", i, i%200);
    } else if (i < 2800) {
      sprintf(spname, "%d; ch %d New Quadratic Mean energy (QDT and lamda) [0.25 keV]", i, i%200);
    } else {
      sprintf(spname, "%d; ch %d Best option for corrected energy          [0.25 keV]", i, i%200);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  return 0;
}

/* ---------------------------------------- */

int CTC2_info_write(MJRunInfo *runInfo, CTC2info *CTC2) {

  int   i;
  FILE  *f_out;

  /*
    write new charge-trapping correction data to file ctc2.output
  */
  if (!(f_out = fopen("ctc2.output", "w"))) {
    printf("\n ERROR: Cannot open file ctc2.output for writing!\n");
    return 1;
  }
  printf("\n Writing charge-trapping correction values to file ctc2.output\n");

  for (i=0; i<200; i++) {
    if (i%100 == 0)
      fprintf(f_out,
              "#Chan  lin_dt_slope  lamda_slope quad_dt_slope quad_dt_quad    lin_dt_gain   lamda_gain     qdt_gain  best_option\n");
    //         %5dxx  %12.4fxxxxxx %12.4fxxxxxx  %12.4fxxxxxx %12.4fxxxxxx   %12.8lfxxxxx %12.8lfxxxxx %12.8lfxxxxx %10dxxxxxx
    if (i%100 < runInfo->nGe)
      fprintf(f_out, "%5d  %12.4f %12.4f  %12.4f %12.4f   %12.8lf %12.8lf %12.8lf %10d\n",
              i, CTC2->e_dt_slope[i], CTC2->e_lamda_slope[i], CTC2->e_qdt_slope[i], CTC2->e_qdt_quad[i],
              CTC2->e_dt_gain[i], CTC2->e_lamda_gain[i], CTC2->e_qdt_gain[i], CTC2->best_ctc2_res[i]);
  }

  fclose(f_out);
  return 0;
} /* CTC_info_write */

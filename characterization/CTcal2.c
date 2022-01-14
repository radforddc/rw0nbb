#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define MAKE_2D 0       // make a file with A/E, drift, and energy data for channel CHAN_2D, for 2d plots
#define CHAN_2D 0      // channel of interest for 2d plotting
#define HIS_COUNT 3200  // number of spectra in his[][] array
#define FWHM_RATIO 0.95 // ratio of FWHM to decide between E_dt and E_lamda charge-trapping correction
#define FIT_CTC_QLC 0   // For qdratic lamda energy fir, switches between e_raw (0) and e_ctc_adc (1)


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
  float  his2[16384]={0};
  int    *his3[1000];
  FILE   *fp, *f_out, *f_out_2d = NULL;

  double s00[2][200]={{0}}, s01[2][200]={{0}}, s10[2][200]={{0}}, s11[2][200]={{0}};
  double s20[2][200]={{0}}, s21[2][200]={{0}}, s30[2][200]={{0}}, s40[2][200]={{0}};
  float  e_qdtc_slope[200];      // linear factor for quadratic drift-time correction of energy
  float  e_qdtc_quad[200]={0};   // quad factor for quadratic drift-time correction of energy
  float  e_qlc_slope[200];       // linear factor for quadratic lamda correction of e_ctc
  float  e_qlc_quad[200]={0};    // quad factor for quadratic lamda correction of e_ctc
  double e_qdtc_gain[200];       // gain for quadratic drift-time correction of energy
  double e_qlc_gain[200];        // gain for quadratic lamda correction of e_ctc
  int    best_qdt_ql[200]={0};   // indicates which option has better resolution
  double e_qdtc, e_qdtc_adc, e_qlc, e_qlc_adc;
  // double dcr2lamda[2][200] = {{0}};


  /* initialize */
  /* malloc and clear histogram space */
  if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
    printf("ERROR in CTcal.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1]+8192;
  if ((his3[0] = calloc(1000*16384, sizeof(int))) == NULL) {
    printf("ERROR in CTcal.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<1000; i++) his3[i] = his3[i-1]+16384;

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
     *            to get: FWHM
     * Step 6: Quadratic fit to: E_raw vs. DT;        E_ctc vs. lamda
     *                   to get: E_qdtc slope, quad;  E_qlc slope, quad
     * Step 7: Histogram: E_qdtc_adc, E_qlc_adc
     *            to get: E_qdtc gain, E_qlc gain
     * Step 8: Histogram: E_qdtc, E_qlc
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
        /*
        if (e_ctc > 2611.5 && e_ctc < 2617.5) {
          if ( dcr > -30 && dcr < 70 && lamda > -2 && lamda < 4) {
            // evaluate sums for linear fit of DCR vs lamda
            s00[2][chan]++;
            s10[2][chan] += lamda;
            s20[2][chan] += lamda * lamda;
            s01[2][chan] += dcr;
            s11[2][chan] += lamda * dcr;
            //} else {
            // printf("chan DCR lamda:  %3d %7.1f %7.2f\n", chan, dcr, lamda);
          }
        }
        */
        continue;       // end of step-1 processing for this event
      }

      /*
      if (dcr2lamda[0][chan] > 0.01) {
        // try using mean of lamda and dcr' for e_lamda
        lamda = (lamda + dcr2lamda[0][chan] * dcr + dcr2lamda[1][chan]) / 2.0;
        e_lamda_adc = e_raw + lamda * CTC.e_lamda_slope[chan];
        e_lamda = e_lamda_adc * gain * CTC.e_lamda_gain[chan];
        if (step == 6 && e_lamda < 4000) {
          his[2800+chan][(int) (2.0*e_lamda + 0.5)]++;
          his[3000+chan][(int) (e_ctc + e_lamda + 0.5)]++;
        }
      }
      */

      if (step == 6 && e_ctc > 2611.5 && e_ctc < 2617.5) {
        // evaluate sums for quadratic fit of e_raw vs. drift time
        s00[0][chan]++;
        s10[0][chan] += drift;
        s20[0][chan] += drift * drift;
        s30[0][chan] += drift * drift * drift;
        s40[0][chan] += drift * drift * drift * drift;
        s01[0][chan] += e_raw;
        s11[0][chan] += e_raw * drift;
        s21[0][chan] += e_raw * drift * drift;
        // evaluate sums for quadratic fit of e_ctc vs. lamda
        if (lamda < -2 || lamda > 4) continue;
        s00[1][chan]++;
        s10[1][chan] += lamda;
        s20[1][chan] += lamda * lamda;
        s30[1][chan] += lamda * lamda * lamda;
        s40[1][chan] += lamda * lamda * lamda * lamda;
        if (FIT_CTC_QLC) {
          s01[1][chan] += e_ctc_adc ;
          s11[1][chan] += e_ctc_adc * lamda;
          s21[1][chan] += e_ctc_adc * lamda * lamda;
        } else {
          s01[1][chan] += e_raw ;
          s11[1][chan] += e_raw * lamda;
          s21[1][chan] += e_raw * lamda * lamda;
        }
        continue;
      }

      //e_qdtc, e_qdtc_adc, e_qlc, e_qlc_adc;
      e_qdtc_adc = e_raw + (e_qdtc_slope[chan] + e_qdtc_quad[chan]*drift) * drift;
      if (FIT_CTC_QLC) {
        e_qlc_adc  = e_ctc_adc + (e_qlc_slope[chan] + e_qlc_quad[chan]*lamda) * lamda;
      } else {
        e_qlc_adc  = e_raw + (e_qlc_slope[chan] + e_qlc_quad[chan]*lamda) * lamda;
      }
      // histogram energies in [ADC] and [0.5 keV] units
      if (step == 7) {
        // histogram energies corrected with new quadratic fit
        if (e_qdtc_adc < 8192) his[600+chan][(int) (e_qdtc_adc + 0.5)]++;
        if (e_qlc_adc < 8192) his[800+chan][(int) (e_qlc_adc + 0.5)]++;
        // histogram energies corrected with new linear fit
        his[2200+chan][(int) (2.0*e_ctc + 0.5)]++;
        his[2400+chan][(int) (2.0*e_lamda + 0.5)]++;
        his[2600+chan][(int) (e_ctc + e_lamda + 0.5)]++;
      } else if (step == 8) {
        e_qdtc = e_qdtc_adc * e_qdtc_gain[chan];
        e_qlc  = e_qlc_adc  * e_qlc_gain[chan];
        if (e_qdtc < 4090) his[1000+chan][(int) (2.0*e_qdtc + 0.5)]++;
        if (e_qdtc < 16380/5) his2[(int) (5.0*e_ctc + 0.5)]++;
        if (e_qlc < 4090) {
          his[1200+chan][(int) (2.0*e_qlc + 0.5)]++;
          his[1400+chan][(int) (e_qdtc + e_lamda + 0.5)]++;
        }
        if (e_qdtc > 1310 && e_qdtc < 2750 && e_qlc > 1310 && e_qlc < 2750) {
          his[1600+chan][(int) (4.0*(e_qdtc-1307.25) + 0.5)]++;
          his[1800+chan][(int) (4.0*(e_qlc-1307.25) + 0.5)]++;
          //his[2000+chan][(int) (2.0*(e_qdtc + e_lamda - 2.0*1307.25) + 0.5)]++;
          his[2000+chan][(int) (2.0*(e_qdtc + e_qlc - 2.0*1307.25) + 0.5)]++;
        }
        if (e_qdtc > 20 && e_qdtc < 2750 && e_qlc > 20 && e_qlc < 2750) {
          his3[    chan][(int) (4.0*(e_ctc)           + 0.5)]++;
          his3[200+chan][(int) (4.0*(e_lamda)         + 0.5)]++;
          his3[400+chan][(int) (2.0*(e_ctc + e_lamda) + 0.5)]++;
          his3[600+chan][(int) (4.0*(e_qdtc)          + 0.5)]++;
          his3[800+chan][(int) (2.0*(e_qdtc + e_qlc)  + 0.5)]++;
        }
      }
      // make file for 2D plots  of E vs CTC
      // roi_elo = DEP_E - 40.0;
      // if (f_out_2d && step == 8 && chan == chan_2d && e_ctc >= roi_elo && e_ctc <= roi_elo+80)
      roi_elo = CAL_E - 40.0;
      if (f_out_2d && step == 8 && chan == chan_2d && e_ctc >= roi_elo && e_ctc <= roi_elo+700) {
        fprintf(f_out_2d, "%4d %9.2f %8.2f %8.2f %8.2f %10.3f %7.3f %7.3f %9.2f\n",
                chan, e_raw*gain, e_ctc, e_qdtc, e_lamda, drift, lamda, dcr, aovere);
      }
    }
    // --------------------------------------------------------------------------
    // end of event processing for this step
    // now process the current histograms to find calibrations, optimal correction factors, etc

    if (step == 5) {
      // extract E_ctc fwhm

      /*
      // extract dcr-vs-lamda slopes
      for (chan=0; chan<200; chan++) {
        if (s00[2][chan] > 500) {
          // numerical errors from dcr*dcr are much larger in s20 than lamda*lamda in s02, so use second set of  numbers
          double a = (s10[2][chan]*s01[2][chan] - s11[2][chan]*s00[2][chan])/(s10[2][chan]*s10[2][chan] - s20[2][chan]*s00[2][chan]);
          double b = (s01[2][chan] - a*s10[2][chan]) / s00[2][chan];
          printf("DCR vs lamda fit:  %3d %7.2f %6.2f   ->  %5.3f %6.3f\n", chan, a, b, 1.0/a, -b/a);
          dcr2lamda[0][chan] = 1.0/a;
          dcr2lamda[1][chan] = -b/a;
        }
      }
      */
      continue;
    }

    if (step == 6) {
      /* analyse fit sums to extract quadratice fits of energies vs drift time and lamda */
      for (chan=0; chan<200; chan++) {
        e_qdtc_slope[chan] = CTC.e_dt_slope[chan];
        e_qlc_slope[chan]  = CTC.e_lamda_slope[chan];
        if (s00[0][chan] > 500) {
          // quadratic fit against drift time:
          double d = (s00[0][chan]*s20[0][chan]*s40[0][chan] - s10[0][chan]*s10[0][chan]*s40[0][chan] -
                      s00[0][chan]*s30[0][chan]*s30[0][chan] + s10[0][chan]*s20[0][chan]*s30[0][chan]*2.0 -
                      s20[0][chan]*s20[0][chan]*s20[0][chan]);
          if (fabs(d) < 0.1) continue;
          e_qdtc_slope[chan] = -(s11[0][chan]*s00[0][chan]*s40[0][chan] - s01[0][chan]*s10[0][chan]*s40[0][chan] +
                                 s01[0][chan]*s20[0][chan]*s30[0][chan] - s21[0][chan]*s00[0][chan]*s30[0][chan] -
                                 s11[0][chan]*s20[0][chan]*s20[0][chan] + s21[0][chan]*s10[0][chan]*s20[0][chan]) / d;
          e_qdtc_quad[chan]  = -(s01[0][chan]*s10[0][chan]*s30[0][chan] - s11[0][chan]*s00[0][chan]*s30[0][chan] -
                                 s01[0][chan]*s20[0][chan]*s20[0][chan] + s11[0][chan]*s10[0][chan]*s20[0][chan] +
                                 s21[0][chan]*s00[0][chan]*s20[0][chan] - s21[0][chan]*s10[0][chan]*s10[0][chan]) / d;
          // linear fit:
          double a = -(s10[0][chan]*s01[0][chan] - s11[0][chan]*s00[0][chan])/
                      (s10[0][chan]*s10[0][chan] - s20[0][chan]*s00[0][chan]);
          printf("%3d %5.1f -> %5.1f  -> %5.1f %6.2f;  ",
                 chan, CTC.e_dt_slope[chan], a, e_qdtc_slope[chan], e_qdtc_quad[chan]);
          CTC.e_dt_slope[chan] = a;

          // quadratic fit against lamda:
          d = (s00[1][chan]*s20[1][chan]*s40[1][chan] - s10[1][chan]*s10[1][chan]*s40[1][chan] -
               s00[1][chan]*s30[1][chan]*s30[1][chan] + s10[1][chan]*s20[1][chan]*s30[1][chan]*2.0 -
               s20[1][chan]*s20[1][chan]*s20[1][chan]);
          if (fabs(d) < 0.1) { printf("\n"); continue; }
          e_qlc_slope[chan] = -(s11[1][chan]*s00[1][chan]*s40[1][chan] - s01[1][chan]*s10[1][chan]*s40[1][chan] +
                                s01[1][chan]*s20[1][chan]*s30[1][chan] - s21[1][chan]*s00[1][chan]*s30[1][chan] -
                                s11[1][chan]*s20[1][chan]*s20[1][chan] + s21[1][chan]*s10[1][chan]*s20[1][chan]) / d;
          e_qlc_quad[chan]  = -(s01[1][chan]*s10[1][chan]*s30[1][chan] - s11[1][chan]*s00[1][chan]*s30[1][chan] -
                                s01[1][chan]*s20[1][chan]*s20[1][chan] + s11[1][chan]*s10[1][chan]*s20[1][chan] +
                                s21[1][chan]*s00[1][chan]*s20[1][chan] - s21[1][chan]*s10[1][chan]*s10[1][chan]) / d;
          // linear fit:
          a = -(s10[1][chan]*s01[1][chan] - s11[1][chan]*s00[1][chan])/
               (s10[1][chan]*s10[1][chan] - s20[1][chan]*s00[1][chan]);
          printf("%5.1f -> %5.1f -> %5.1f %6.2f\n",
                 CTC.e_lamda_slope[chan], a, e_qlc_slope[chan], e_qlc_quad[chan]);
          CTC.e_lamda_slope[chan] = a;
        }
      }

    } else if (step == 7) {
      /* find position of 2614.5-keV peak and use that to compute the calibration gainsc */
      FILE *fp = fopen("fwhm_ctc2.txt", "w");
      fprintf(fp, "#chan   QDTC    QLC\n");
      float fwhm0 = 0;
      for (chan=0; chan<200; chan++) {
        e_qdtc_gain[chan] =  e_qlc_gain[chan] = 1;
        if (s00[0][chan] < 501) continue;
        best_qdt_ql[chan] = 0;
        if (chan%100 >= runInfo.nGe) continue;
        if (chan%100 == 0) printf("\n");
        // first for E_qdtc
        fwhm = 5;
        j = 4000;
        if (runInfo.flashcam == 3) j = 3000;  // HADES data has lower effective gain
        if (chan > 99) {
          fwhm = 3;
          j = 1700;
        }
        if ((pos = autopeak3(his[600+chan], j, 8000, &area, &fwhm))) {
          fwhm0 = fwhm;
          if (chan < 100) {
            e_qdtc_gain[chan] = Dets[chan].HGcalib[2] = CAL_E/pos;
            printf("Ch %3d %6.1f:   E_qdtc_adc  P = %7.1f A = %7.0f  FWHM = %7.2f keV\n",
                   chan, CAL_E, pos, area, fwhm*Dets[chan].HGcalib[2]);
            fprintf(fp, "%4d %7.3f", chan, fwhm*Dets[chan].HGcalib[2]);
          } else {
            e_qdtc_gain[chan] = Dets[chan-100].LGcalib[2] = CAL_E/pos;
            printf("Ch %3d %6.1f:   E_qdtc_adc  P = %7.1f A = %7.0f  FWHM = %7.2f keV\n",
                   chan, CAL_E, pos, area, fwhm*Dets[chan-100].LGcalib[2]);
            fprintf(fp, "%4d %7.3f", chan, fwhm*Dets[chan-100].LGcalib[2]);
          }
        }

        // now for E_qlc
        fwhm = 5;
        j = 4000;
        if (runInfo.flashcam == 3) j = 3000;  // HADES data has lower effective gain
        if (chan > 99) {
          fwhm = 3;
          j = 1700;
        }
        if ((pos = autopeak3(his[800+chan], j, 8000, &area, &fwhm))) {
          if (fwhm < fwhm0 * FWHM_RATIO) best_qdt_ql[chan] = 1;
          if (chan < 100) {
            e_qlc_gain[chan] = Dets[chan].HGcalib[3] = CAL_E/pos;
            printf("                 E_qlc_adc   P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   pos, area, fwhm*Dets[chan].HGcalib[3]);
            fprintf(fp, " %7.3f\n", fwhm*Dets[chan].HGcalib[3]);
          } else {
            e_qlc_gain[chan] = Dets[chan-100].LGcalib[3] = CAL_E/pos;
            printf("                 E_qlc_adc   P = %7.1f A = %7.0f; FWHM = %7.2f keV\n",
                   pos, area, fwhm*Dets[chan-100].LGcalib[3]);
            fprintf(fp, " %7.3f\n", fwhm*Dets[chan-100].LGcalib[3]);
          }
        }
      }
      fclose(fp);
    }
  }

  printf(">>  All done...\n"); fflush(stdout);

  // write energy calibrations to gains.output
  if ((f_out = fopen("gains2.output","w"))) {
    printf("\n Writing energy calibrations to gains.output\n");
    for (i=0; i<runInfo.nGe; i++)
      fprintf(f_out,
              "%3d %10.8lf %10.8lf %s\n",
              i, Dets[i].HGcalib[0], Dets[i].LGcalib[0], Dets[i].StrName);
    fclose(f_out);
  }

  // write charge-trapping data to ctc.output
  // CTC_info_write(&runInfo, &CTC);

  // write out histograms
  f_out = fopen("ctc2.rms", "w");
  for (i=0; i<HIS_COUNT; i++) {
   char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d DT-corrected energy [0.5 keV]", i, i);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d lamda-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d mean linear corrected energy [0.5 keV]", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d QDT-corrected energy [ADC]", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d QL-corrected energy [ADC]", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d QDT-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d QL-corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d mean quad corrected energy [0.5 keV]", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; ch %d QDT-corrected energy [0.25 keV]", i, i%200);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d QL-corrected energy [0.25 keV]", i, i%200);
    } else {
      sprintf(spname, "%d; ch %d mean quad corrected energy [0.25 keV]", i, i%200);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);
  
  f_out = fopen("ctc3.rms", "w");
  for (i=0; i<1000; i++) {
   char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d DT-corrected energy [0.25 keV]", i, i);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d lamda-corrected energy [0.25 keV]", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d mean linear corrected energy [0.25 keV]", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d QDT-corrected energy [0.25 keV]", i, i%200);
    } else {
      sprintf(spname, "%d; ch %d mean quad corrected energy [0.25 keV]", i, i%200);
    }
    write_his(his3[i], 16384, i, spname, f_out);
  }
  fclose(f_out);
  
  f_out = fopen("ctc.dat", "w");
  fwrite(his2, 1, sizeof(his2), f_out);
  fclose(f_out);
  
  return 0;
}

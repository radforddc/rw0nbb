#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define MAKE_2D 1        // make a file with A/E, drift, and energy data for channel CHAN_2D, for 2d plots
#define CHAN_2D 0        // channel of interest for 2d plotting
#define HIS_COUNT 2400   // number of spectra in his[][] array

#define SUBTR_DCR_MEAN (e_adc < 8000 && \
                        ((chan < 100 && SUBTRACT_MEAN_DCR_HG) || \
                         (chan > 99  && SUBTRACT_MEAN_DCR_LG)))

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1, keep_ae_cut = 0, keep_pos[200];


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname_in\n\n", argv[0]);
    return -1;
  }
  /* open skim data file as input */
  while (argn < argc && argv[argn][0] == '-') argn++;
  for (int i=1; i<argc; i++) {
    if (argv[i][0] == '-' && argv[i][1] == 'a') {
      // -a flag: preserve A/E cut relative to A/E position
      keep_ae_cut = 1;
      printf(" %s: Preserving A/E cut relative to A/E position\n", argv[i]);
    }
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

  CTCinfo CTC;
  PSAinfo PSA;
  PZinfo  PZI;

  // data used, stored, and reused in the different steps
  SavedData **sd;
  int      chan;
  float    aovere, drift;
  double   e_raw;
  int      nsd = 0, isd = 0;  // number of saved data, and pointer to saved data id

  float    aovere_norm, a_e_pos[200], a_e_pos1[200];
  float    dcr, lamda;

  double e_ctc, e_adc, gain;
  float  pos, area, fwhm, s1, s2, s3, s4, s5, ppos[200][6] = {{0}};
  int    i, j, k, kk, n, roi_elo;
  int    *his[HIS_COUNT];
  int    *dcr_mean[200], mean_dcr_ready = 0;
  FILE   *f_out, *f_out_2d = 0, *fp;


  /* initialize */
  /* malloc and clear histogram space */
  if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
    printf("ERROR in PSAcal.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1] + 8192;
  if ((dcr_mean[0] = calloc(200*8192, sizeof(int))) == NULL) {
    printf("ERROR in PSAcal.c; cannot malloc dcr_mean!\n");
    exit(-1);
  }
  for (i=1; i<200; i++) dcr_mean[i] = dcr_mean[i-1] + 8192;

  // read saved skim data from f_in
  fread(&nsd, sizeof(int), 1, f_in);
  fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
  fread(&runInfo, sizeof(runInfo), 1, f_in);
  /* malloc space for SavedData */
  if ((sd = malloc(nsd*sizeof(*sd))) == NULL ||
      (sd[0] = malloc(nsd*sizeof(SavedData))) == NULL) {
    printf("ERROR in PSAcal.c; cannot malloc SavedData!\n");
    exit(-1);
  }
  for (i=1; i<nsd; i++) sd[i] = sd[i-1] + 1;
  fread(*sd, sizeof(SavedData), nsd, f_in);
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

  /* read PZ orrection info */
  if (!PZ_info_read(&runInfo, &PZI)) {
    printf("\n ERROR: No initial pole-zero data read. Does PZ.input exist?\n");
    exit(-1);
  }
  /* read energy correction factors from ctc.input */
  if (!CTC_info_read(&runInfo, &CTC)) {
    printf("\n ERROR: No initial charge-trapping correction data read. Does ctc.input exist?\n");
    exit(-1);
  }
  /* read A/E, DCR, and lamda values from psa.input */
  /* also read individual trapezoid values from filters.input (if it exists) */
  if (!PSA_info_read(&runInfo, &PSA)) {
    printf("\n ERROR: No initial PSA data read. Does psa.input exist?\n");
    if (keep_ae_cut) {
      keep_ae_cut = 0;
      printf("Will NOT preserve A/E cut relative to A/E position\n");
    }
  } else {
    for (i=0; i<200; i++) keep_pos[i] = PSA.ae_pos[i];
  }
  /* read mean DCR and lamda values vs. E for residual INL correction */
  if (SUBTRACT_MEAN_DCR_HG || SUBTRACT_MEAN_DCR_LG) {
    if ((fp = fopen("INL_DCR_input.sec", "r"))) {
      fread(dcr_mean[0], 8192*sizeof(int), 200, fp);        // get or skip smoothed values
      if (SUBTRACT_MEAN_DCR_HG > 1) {
        fread(dcr_mean[0], 8192*sizeof(int), 100, fp);      // get unsmoothed HG values
        if (SUBTRACT_MEAN_DCR_LG > 1)
          fread(dcr_mean[100], 8192*sizeof(int), 100, fp);  // get unsmoothed LG values
      } else if (SUBTRACT_MEAN_DCR_LG > 1) {
        fread(dcr_mean[100], 8192*sizeof(int), 100, fp);    // skip unsmoothed HG values
        fread(dcr_mean[100], 8192*sizeof(int), 100, fp);    // get unsmoothed LG values
      }
      fclose(fp);
      mean_dcr_ready = 1;   // all is ready
    } else {
      mean_dcr_ready = -1;  // -1 indicates need to take further action later
    }
  }

  if (MAKE_2D &&
      ((CHAN_2D < 100 && Dets[CHAN_2D].HGChEnabled) ||
       (CHAN_2D >  99 && Dets[CHAN_2D%100].LGChEnabled))) {
    char fname[64];
    sprintf(fname, "PSA_ch%3.3d_2d.dat", CHAN_2D);
    f_out_2d = fopen(fname, "w");
    fprintf(f_out_2d, "#chan    E_ctc     A/E    A/E_raw    DT   DT_corr  DCR   lamda  t90-100\n");
  }

  // end of initialization
  // start loop over events in skim data

  // ---------------------- steps 1-5 ----------------------
  for (int step = 1; step <= 5; step++) {
    printf(" ******************** Step %d of 5 ********************\n", step);

    /*
     * Step 1: Histogram: Raw A/E
     *            to get: Raw A/E pos
     *         Also create averaged DCR and lamda values as a function of E_raw
     * Step 2: Histogram: A/E_ctc vs. DT, DCR vs. DT, lamda vs. DT
     *            to get: A/E_ctc slope and pos, DCR_DT slope, lamda_DT slope
     * Step 3: Histogram: A/E_ctc_e vs. E,  DCR', lamda',  A/E_ctc
     *            to get: A/E_ctc_e slope and pos
     * Step 4: Histogram: A/E_ctc_e gated by E_ctc      Also final A/E_ctc_e
     *            to get: A/E_ctc_e cut value
     * Step 5: Histogram: E_ctc cut by A/E_ctc_e
     *            to get: Final spectra, A/E acceptances
     *         Also make file for 2D plots of A/E|E|CTC if required
     */

    for (isd = 0; isd < nsd; isd++) {
      // if (isd%(nsd/10) == 0) printf(">>  event %7d (%d/10\n", isd, isd*10/nsd);
      chan   = sd[isd]->chan;
      e_raw  = sd[isd]->e;
      drift  = sd[isd]->drift;
      aovere = sd[isd]->a_over_e;
      dcr    = sd[isd]->dcr;
      lamda  = sd[isd]->lamda;
      // if (sd[isd]->t100 - sd[isd]->t90 > 15) continue;

      e_adc = e_raw + drift*CTC.e_dt_slope[chan];
      lamda *= 8.0;                              // scaled to roughly match DCR at 2 MeV
      // lamda *= e_adc/5000.0;                  // energy scaling that would roughly match DCR
      // if (e_adc > 50) dcr *= 5000.0/e_adc;    // energy scaling that would roughly match lamda
 
      if (step == 1) {
        // histogram raw A/E
        if ((j = aovere + 0.5) < 2000 && j > 0) his[chan][j]++;

        // make averaged DCR and lamda, as a function of e_adc
        dcr -= PSA.dcr_dt_slope[chan] * drift;
        lamda -= PSA.lamda_dt_slope[chan] * drift;
        if (e_adc < 8000 && fabs(dcr) < 50.0) {
          his[1800+chan][(int) e_adc/2] += dcr;
          his[2000+chan][(int) e_adc/2]++;
        }
        if (e_adc < 8000 && fabs(lamda) < 50.0) {
          his[1800+chan][4000 + (int) e_adc/2] += lamda;
          his[2000+chan][4000 + (int) e_adc/2]++;
        }
        continue;
      }

      if (chan < 100) {
        gain = Dets[chan].HGcalib[0];
      } else {
        gain = Dets[chan-100].LGcalib[0];
      }
      e_ctc = e_adc * gain;

      float dtc = drift*3.0*2614.0/e_adc;   // in (x)us, for DT- correction to A/E
      if (chan%100 == 50) dtc *= 3.0;       // special hack for detector 50; has especially strong variation
      float ec = e_ctc/1000.0 - 2.0;        // in MeV, for E-correction to A/E
      float aovere1 = aovere + PSA.ae_dt_slope[chan] * dtc;
      float aovere2 = aovere1 + PSA.ae_e_slope[chan] * ec;

      /* find optimum drift-time corrections for A/E and energy */
      roi_elo =  1574;
      if (e_ctc >= roi_elo && e_ctc <= roi_elo+800) {  // 800 keV window covering ROI
        if (step == 2) {
          /* try 20 options to find optimum drift-time correction for A/E_ctc */
          aovere_norm = aovere * 800.0/a_e_pos[chan];
          s1 = aovere_norm - 600.0 - 1.0 * dtc;
          for (k=0; k < 20; k++) {
            if ((j = s1 + 0.5) > 0 && j < 400) his[200+chan][j + 400*k]++;
            s1 += 0.5 * dtc;
          }
        } else if (step == 3) {
          /* try 20 different slopes of A/E_ctc_e against energy to see which one is best
             each step adds a slope of 2 A/E per MeV
          */
          aovere_norm = aovere1 * 800.0/a_e_pos1[chan];
          s1 = aovere_norm - 600.0 - 2.0*1.4 * ec;
          for (k=0; k < 20; k++) {
            if ((j = s1 + 0.5) > 0 && j < 400) his[400+chan][j + 400*k]++;
            s1 += 2.0*0.7 * ec;
          }
        }
      }

      if (step == 4) {
        aovere_norm = aovere * 800.0/a_e_pos[chan];
        if ((j = aovere_norm + 0.5) < 2000 && j > 0) his[chan][1500+j]++;
        aovere_norm = aovere1 * 800.0/a_e_pos1[chan];
        if ((j = aovere_norm + 0.5) < 2000 && j > 0) his[chan][3000+j]++;
        aovere_norm = aovere2 * 800.0/PSA.ae_pos[chan];
        if ((j = aovere_norm + 0.5) < 2000 && j > 0) his[chan][4500+j]++;
      } else if (step == 5) {
        aovere_norm = aovere2 * 800.0/PSA.ae_pos[chan];
        if ((j = aovere_norm + 0.5) < 2000 && j > 0) his[chan][6000+j]++;

        // make file for 2D plots  of A/E|E|CTC
        if (f_out_2d && chan == CHAN_2D && e_ctc >= roi_elo)// && e_ctc <= roi_elo+700)
          fprintf(f_out_2d, "%4d %9.3f %8.2f %8.2f %7.2f %7.2f %7.2f %7.2f %6d\n",
                  chan, e_ctc, aovere_norm, aovere * 800.0/a_e_pos[chan], drift, dtc,
                  sd[isd]->dcr, sd[isd]->lamda/1.3, sd[isd]->t100 - sd[isd]->t90);
      }

      /* find and test A/E cut */
      if (step == 4) {
        // put narrow gate on DE peak, and wider gate on bgnd
        if ((j = aovere2 * 800.0/PSA.ae_pos[chan] + 0.5) < 1500 && j > 0) {
          if (e_ctc >= DEP_E - 2.5 && e_ctc < DEP_E + 2.5) {            // DEP; 5.0 keV
            his[800+chan][j+3000]++;
          } else if ((e_ctc >= DEP_E - 28.0 && e_ctc < DEP_E - 6.0) ||
                     (e_ctc >= DEP_E + 6.0  && e_ctc < DEP_E + 24.0)) { // backgnd; 22+18 keV
            his[800+chan][j+1500]++;
          } else if (e_ctc >= SEP_E - 2.5 && e_ctc < SEP_E + 2.5) {     // SEP; 5 keV
            his[800+chan][j+6000]++;
          } else if ((e_ctc >= SEP_E - 30.0 && e_ctc < SEP_E - 6.0) ||
                     (e_ctc >= SEP_E + 6.0 && e_ctc < SEP_E + 22.0)) { // backgnd; 24+16 keV
            his[800+chan][j+4500]++;
          }
        }

      } else if (step == 5) {
        s2 = aovere2;
        // adjust for energy dependence of cut
        // first correct for series noise: assume 2*sigma cut
        // series noise contribution to A/E sigma = BL_RMS * sqrt(2*rise) * factor / E_raw
        if (AOE_CORRECT_NOISE)
          s2 -= (2.0 * PZI.bl_rms[chan] * sqrt(2.0 * (float) PSA.a_e_rise[chan]) *
                 PSA.a_e_factor[chan] * gain/1593.0 * (1.0 - 1593.0/e_ctc));        // 1593 keV = DEP
        // now deal with energy dependence of variation in A/E due to bremsstrahlung etc
        s2 += AOE_CORRECT_EDEP * (PSA.ae_pos[chan] - PSA.ae_cut[chan]) * (1.0 - 1593.0/e_ctc);  // FIXME: Add limit at low e_ctc

        // count events in narrow gate on DEP and SEP, and wider gates on bgnd
        if (e_ctc >= DEP_E - 2.5 && e_ctc < DEP_E + 2.5) {            // DEP; 5.0 keV
          his[800+chan][1]++;
          if (s2 >= PSA.ae_cut[chan]) his[800+chan][2]++;
        } else if ((e_ctc >= DEP_E - 28.0 && e_ctc < DEP_E - 6.0) ||
                   (e_ctc >= DEP_E + 6.0  && e_ctc < DEP_E + 24.0)) { // backgnd; 22+18 keV
          his[800+chan][3]++;
          if (s2 >= PSA.ae_cut[chan]) his[800+chan][4]++;
        } else if (e_ctc >= SEP_E - 2.5 && e_ctc < SEP_E + 2.5) {     // SEP; 5 keV
          his[800+chan][5]++;
          if (s2 >= PSA.ae_cut[chan]) his[800+chan][6]++;
        } else if ((e_ctc >= SEP_E - 30.0 && e_ctc < SEP_E - 6.0) ||
                   (e_ctc >= SEP_E + 6.0 && e_ctc < SEP_E + 22.0)) { // backgnd; 24+16 keV
          his[800+chan][7]++;
          if (s2 >= PSA.ae_cut[chan]) his[800+chan][8]++;
        } else if ((CAL_E > 2610 && e_ctc >= 2010.0 && e_ctc < 2070.0) ||  // continuum at ROI (Th-228)
                   (CAL_E < 2600 && e_ctc >= 2041.0 && e_ctc < 2081.0))  { // continuum at ROI (Co-56)
          his[800+chan][9]++;
          if (s2 >= PSA.ae_cut[chan]) his[800+chan][10]++;
        }
        if (a_e_pos[chan] > 100 &&
            e_ctc >= 1510 && e_ctc <= 2500) { // DE and SE peaks
          int e2 = (e_ctc + 0.5 - 1500.0);
          his[600+chan][e2]++;
          if (s2 >= PSA.ae_cut[chan])   his[600+chan][1000 + e2]++;
          if (s2 >= PSA.ae_cut[chan]-1) his[600+chan][2000 + e2]++;
          if (s2 >= PSA.ae_cut[chan]+1) his[600+chan][3000 + e2]++;
        }
      }

      if (step == 2) {
        /* Step 2: Histogram: DCR vs. DT in his[1000], lamda vs. DT in his[1200]
                   to get DCR_DT slope, lamda_DT slope */
        if ((s1 = 100.5 + dcr) > 0 && s1 < 200) {
          // try multiple correction factors for DT dependence of DCR
          for (j=1; j<35; j++) {
            if (s1 > 0 && s1 < 200) his[1000+chan][j*200 + (int) s1]++;
            s1 -= drift;
          }
        }
        // lamda is taken from the 1/tau value fitted with frac2 fixed
        if ((s1 = 100.5 + lamda) > 0 && s1 < 200) {
          // try multiple DT correction factors for preamp lamda
          for (j=1; j<35; j++) {
            if (s1 > 0 && s1 < 200) his[1200+chan][j*200 + (int) s1]++;
            s1 -= 0.8 * drift;
          }
        }

      } else if (step == 3) {
        /* Step 3: Histogram: DCR_DT, lamda_DT in his[1600], to get peak positions for alignment */
        if ((s1 = 500.5 + dcr) > 200 && s1 < 1100) his[1600+chan][(int) s1]++;
        s1 -= PSA.dcr_dt_slope[chan] * drift;
        if (s1 > 200 && s1 < 1100) his[1600+chan][1000 + (int) s1]++;
        if (SUBTR_DCR_MEAN) s1 -= (float) dcr_mean[chan][(int) e_adc/2] / 10.0;        // correct for residual INL
        if (s1 > 200 && s1 < 1100) his[1600+chan][2000 + (int) s1]++;

        if ((s1 = 500.5 + lamda) > 200 && s1 < 1100) his[1600+chan][4000 + (int) s1]++;
        s1 -= PSA.lamda_dt_slope[chan] * drift;
        if (s1 > 200 && s1 < 1100) his[1600+chan][5000 + (int) s1]++;
        if (SUBTR_DCR_MEAN) s1 -= (float) dcr_mean[chan][4000 + (int) e_adc/2] / 10.0; // correct for residual INL
        if (s1 > 200 && s1 < 1100) his[1600+chan][6000 + (int) s1]++;

      } else if (step == 4) {
        /* Step 3: Histogram: final DCR, lamda, moved to cut value, in his[1400] */
        // s1 = 500.5 + dcr - PSA.dcr_dt_slope[chan] * drift - PSA.dcr_lim[chan];
        if ((s1 = 500.5 + dcr) > 200 && s1 < 1100) his[1400+chan][(int) (s1 - ppos[chan][0])]++;
        s1 -= PSA.dcr_dt_slope[chan] * drift;
        if (s1 > 200 && s1 < 1100) his[1400+chan][1000 + (int) (s1 - ppos[chan][1])]++;
        if (SUBTR_DCR_MEAN) s1 -= (float) dcr_mean[chan][(int) e_adc/2] / 10.0;        // correct for residual INL
        if (s1 > 200 && s1 < 1100) his[1400+chan][2000 + (int) (s1 - ppos[chan][2])]++;
        if (s1 > 200 && s1 < 1100) his[1400+chan][3000 + (int) (s1 - PSA.dcr_lim[chan])]++;

        // s1 = 500.5 + lamda - PSA.lamda_dt_slope[chan] * drift - PSA.lamda_lim[chan];
        if ((s1 = 500.5 + lamda) > 200 && s1 < 1100) his[1400+chan][4000 + (int) (s1 - ppos[chan][3])]++;
        s1 -= PSA.lamda_dt_slope[chan] * drift;
        if (s1 > 200 && s1 < 1100) his[1400+chan][5000 + (int) (s1 - ppos[chan][4])]++;
        if (SUBTR_DCR_MEAN) s1 -= (float) dcr_mean[chan][4000 + (int) e_adc/2] / 10.0; // correct for residual INL
        if (s1 > 200 && s1 < 1100) his[1400+chan][6000 + (int) (s1 - ppos[chan][5])]++;
        if (s1 > 200 && s1 < 1100) his[1400+chan][7000 + (int) (s1 - PSA.lamda_lim[chan])]++;
      }
    }

    // --------------------------------------------------------------------------
    // end of event processing for this step
    // now process the current histograms to find calibrations, optimal correction factors, etc

    if (step ==1) {
      /* find position of raw A/E peak */
      for (chan = 0; chan < 200; chan++) {
        if (chan%100 >= runInfo.nGe) continue;
        // look for the most counts over a 7-bin interval
        j = 4;
        k = kk = 0;
        for (i=0; i<9; i++) k += his[chan][i];
        for (i=5; i<1990; i++) {
          k += his[chan][i+4] - his[chan][i-5];
          if (kk < k) {
            j = i;
            kk = k;
          }
        }
        if (kk > 100) a_e_pos[chan] = j;

        // calculate mean values of DCR and lamda as a function of energy
        for (i=0; i<8192; i++) {
          if (his[2000+chan][i] > 5) {
            his[1800+chan][i] = lrint(10.0 * (float) his[1800+chan][i] / (float) his[2000+chan][i]);
          } else {
            his[1800+chan][i] = 0;
          }
        }
        // calculate a smoothed value from a running average over a 7-bin interval
        for (i=10; i<8182; i++) {
          double sum = 0;
          for (k = -3; k <4; k++) sum += his[1800+chan][i+k];
          his[2200+chan][i] = lrint(sum/7.0);
        }
      }

      // copy mean DCR and lamda vs. energy to proper array
      if (mean_dcr_ready < 0 && SUBTRACT_MEAN_DCR_HG) {
        for (chan = 0; chan < 100; chan++) {
          for (i=0; i<8192; i++)
            dcr_mean[chan][i] = his[1800 + chan + 400*(2-SUBTRACT_MEAN_DCR_HG)][i];
        }
      }
      if (mean_dcr_ready < 0 && SUBTRACT_MEAN_DCR_LG) {
        for (chan = 100; chan < 200; chan++) {
          for (i=0; i<8192; i++)
            dcr_mean[chan][i] = his[1800 + chan + 400*(2-SUBTRACT_MEAN_DCR_LG)][i];
        }
      }

    } else if (step == 2 || step == 3) {
      int sp_offset;
      if (step == 2) {
        printf(" A/E results:\n Chan DT-slope   pos   (peak)\n");
        sp_offset = 200;  // use data from his[200+chan]
      } else {
        printf(" A/E results:\n Chan  E-slope   pos   (peak)\n");
        sp_offset = 400;  // use data from his[400+chan]
      }
      for (chan = 0; chan < 200; chan++) {
        if (chan%100 >= runInfo.nGe) continue;
        j = 3;
        k = kk = 0;
        // look for test peak with the most counts over a 7-bin interval
        for (i=0; i<7; i++) k += his[sp_offset+chan][i];
        for (i=4; i<8000; i++) {
          k += his[sp_offset+chan][i+3] - his[sp_offset+chan][i-4];
          if (kk < k) {
            j = i;
            kk = k;
          }
        }
        if (kk > 100) {
          if (step == 2) {
            /* found slope of A/E with drift time that optimizes resolution */
            a_e_pos1[chan]        = a_e_pos[chan] * (j%400 + 600) / 800.0;
            PSA.ae_dt_slope[chan] = a_e_pos[chan] * (j/400 - 2) / 800.0 * 0.5;
            printf("%3d %6.1f %4.0f (%4d)\n",
                   chan, PSA.ae_dt_slope[chan], a_e_pos1[chan], j/400);
          } else if (step == 3) {
            /* found slope of A/E with energy that optimizes resolution */
            PSA.ae_pos[chan]     = a_e_pos1[chan] * (j%400 + 600) / 800.0;
            PSA.ae_e_slope[chan] = a_e_pos1[chan] * (j/400 - 2) / 800.0 * 0.7;
            printf("%3d %6.1f %4.0f (%4d)\n",
                   chan, PSA.ae_e_slope[chan], PSA.ae_pos[chan], j/400);
          }
        }
      }
    }

    if (step == 4) {
      /* try to find the 90% DEP cut on A/E */
      double s6=0, s7=0;
      for (chan = 0; chan < 200; chan++) {
        if (chan%100 >= runInfo.nGe) continue;
        s1 = s2 = 0;
        for (i=1; i<1500; i++) {
          s1 += his[800+chan][i+3000];  // DE peak gate
          s2 += his[800+chan][i+1500];  // 8 * backgnd gate
        }
        if (s1 < 100 || s2 < 20 || s2 > 3*s1) continue;
        s3 = 0.1 * (s1 - s2/8.0);  // 10% of DEP counts
        s1 = 0;
        for (i=1; i<1500 && s1 < s3; i++) {
          s2 = s1;
          s1 += his[800+chan][i+3000] - his[800+chan][i+1500]/8.0;  // DEP pk - bkgnd
        }
        s3 = (float) i - 1.5 + (s3-s2)/(s1-s2); //    1.5 is a slight fudge; should be 2.0
        if (!keep_ae_cut) {
          PSA.ae_cut[chan] = s3 * PSA.ae_pos[chan] / 800.0;
          printf(" chan, cut = %3d %8.2f  ;", chan, PSA.ae_cut[chan]);
        }

        // calculate and report cut acceptance values
        s1 = s2 = s4 = s5 = 0;
        for (i=1500; i>(int)s3; i--) {
          s2 = s1;
          s1 += his[800+chan][i+6000] - his[800+chan][i+4500]/8.0;  // SEP pk - bkgnd
          s5 = s4;
          s4 += his[800+chan][i+3000] - his[800+chan][i+1500]/8.0;  // DEP pk - bkgnd
        }
        // derivatives
        s6 += (s1-s2); s7 += (s4-s5);
        s2 += (ceil(s3) - s3) * (s1-s2);
        s5 += (ceil(s3) - s3) * (s4-s5);
        for (; i>1; i--) {
          s1 += his[800+chan][i+6000] - his[800+chan][i+4500]/8.0;  // SEP pk - bkgnd
          s4 += his[800+chan][i+3000] - his[800+chan][i+1500]/8.0;  // SDP pk - bkgnd
        }
        printf(" SEP, DEP acceptance:  %.4f  %.4f\n", s2/s1, s5/s4);
      }
      if (s6 > 0 && s7 > 0)
        printf("\n  >>>>>  Mean ratio of acceptance derivatives: %.4f  %.4f\n\n",
               s6/s7, s7/s6);

      /* try to find better value for final A/E position */
      for (chan = 0; chan < 200; chan++) {
        if (chan%100 >= runInfo.nGe) continue;
        s2 = s3 = -10;
        for (k = 10; k >= -10; k--) {
          s3 = s2;
          s1 = s2 = 0;
          for (i=760; i<831; i++) {
            s1 += his[chan][i+k+4500];
            s2 += (i-793) * his[chan][i+k+4500];
          }
          if (s1 < 100) break;
          s2 /= s1; // centroid of A/E distrbution, minus 793
          if (s2 > 0) break;
        }
        if (s1 < 100) continue;
        PSA.ae_pos[chan] *= (1.0 + (s2/(s2-s3) + (float) k)/800);
        printf("chan %3d |  k, s2 = %3d %6.3f -> %3d %6.3f;  pos adjustment = %6.2f\n",
               chan, k+1, s3, k, s2, s2/(s2-s3) + (float) k);
        if (keep_ae_cut) {
          PSA.ae_cut[chan] += PSA.ae_pos[chan] - keep_pos[chan];
        }
      }
    }

    /* find optimum DT corrections for DCR and lamda */
    if (step == 2) {
      printf("             DCR/DT             lamda/DT\n");
      printf("Chan | slope area  FWHM  | slope area  FWHM |\n");
      for (chan = 0; chan < 200; chan++) {
        if (chan%100 >= runInfo.nGe) continue;
        // first make sure this channel is working by counting up DCR peak
        s1 = 0;
        for (i=200; i<400; i++) s1 += his[1000+chan][i];
        if (s1 < 100) continue;
        printf(" %3d | ", chan);

        /* find best choice (minimum fwhm) for DTC to DCR value */
        for (int dcr_or_lamda = 0; dcr_or_lamda < 2; dcr_or_lamda++) {
          j = k = n = 0;
          s1 = 999;
          for (j=1; j<35; j++) {
            fwhm = 8;
            if (dcr_or_lamda < 1 && chan > 99) fwhm = 3;
            if ((pos = autopeak3(his[1000 + 200*dcr_or_lamda + chan], 200*j, 200+200*j, &area, &fwhm)) &&
                area > 100 && fwhm < s1) {
              s1 = fwhm;
              s2 = pos;
              k = area;
              n = j;
            }
          }
          if (k > 99) {
            printf(" %2d %6d  %5.2f | ", n, k, s1);
            if (dcr_or_lamda < 1) {
              PSA.dcr_dt_slope[chan] = n-1;
              PSA.dcr_lim[chan] = s2 - (float) (100 + 200*n);
            } else {
              PSA.lamda_dt_slope[chan] = 0.8 * (n-1);
              PSA.lamda_lim[chan] = s2 - (float) (100 + 200*n);
            }
          } else {
            printf("                  | ");
          }
        }
        printf("\n");
      }
    }

    /* find DCR and lamda acceptance limits*/
    if (step == 3) {
      printf("\n                DCR                              lamda\n");
      printf("Chan |   pos    area   FWHM      cut  |   pos    area   FWHM      cut  |\n");
      for (chan = 0; chan < 200; chan++) {
        for (j = 0; j < 3; j++) {
          fwhm = 8;
          area = 0;
          if (chan > 99) fwhm = 3;
          pos = autopeak3(his[1600 + chan], 400+1000*j, 700+1000*j, &area, &fwhm) - 500.0 - 1000.0*j;
          if (area > 200) ppos[chan][j] = pos;
        }
        if (area < 200) continue;
        // set cut limit to 1.3 FWHM above centroid
        PSA.dcr_lim[chan] = 1.3 * fwhm + pos;
        printf(" %3d | %5.1f %7.0f %6.2f %8.1f  |", chan, pos, area, fwhm, PSA.dcr_lim[chan]);

        for (j = 3; j < 6; j++) {
          fwhm = 8;
          area = 0;
          pos = autopeak3(his[1600 + chan], 1400+1000*j, 1700+1000*j, &area, &fwhm) - 1500.0 - 1000.0*j;
          if (area > 200) ppos[chan][j] = pos;
        }
        if (area < 200) continue;
        fwhm = 8;
        pos = autopeak3(his[1600 + chan], 6400, 6700, &area, &fwhm) - 6500.0;
        if (area < 200) continue;
        // set cut limit to 1.3 FWHM above centroid
        PSA.lamda_lim[chan] = 1.3 * fwhm + pos;
        printf(" %5.1f %7.0f %6.2f %8.1f  |\n", pos, area, fwhm, PSA.lamda_lim[chan]);
        // printf("   ppos: %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n",
        //        ppos[chan][0], ppos[chan][1], ppos[chan][2], ppos[chan][3], ppos[chan][4], ppos[chan][5]);
      }
    }

    if (step == 5) {
      // re-calculate and report cut acceptance values
      fp = fopen("aoe_eff.txt", "w");
      fprintf(fp, "#chan  SEP   err    DEP   err   Continuum err\n");
      double s5 = 0, s6=0, s7=0, s8=0, s9=0;
      double s0, e0, e1, e2, e3, e4, e5;
      for (chan = 0; chan < 100; chan++) {  // HG channels ony!
        if (chan%100 >= runInfo.nGe) continue;
        if (chan%100 == 30 || chan%100 == 50) continue;  // FIXME!  - dets with bad A/E
        s1 = (double) his[800+chan][1] - (double) his[800+chan][3]/8.0;  // DEP - bknd, all
        s2 = (double) his[800+chan][2] - (double) his[800+chan][4]/8.0;  // DEP - bknd, cut
        s3 = (double) his[800+chan][5] - (double) his[800+chan][7]/8.0;  // SEP - bknd, all
        s4 = (double) his[800+chan][6] - (double) his[800+chan][8]/8.0;  // SEP - bknd, cut
        e1 = (double) his[800+chan][1] + (double) his[800+chan][3]/64.0;
        e2 = (double) his[800+chan][2] + (double) his[800+chan][4]/64.0;
        e3 = (double) his[800+chan][5] + (double) his[800+chan][7]/64.0;
        e4 = (double) his[800+chan][6] + (double) his[800+chan][8]/64.0;
        s0 = (double) (his[800+chan][1] - his[800+chan][2]) -
             (double) (his[800+chan][3] - his[800+chan][4])/8.0;
        e0 = (double) (his[800+chan][1] - his[800+chan][2]) +
             (double) (his[800+chan][3] - his[800+chan][4])/64.0;
        if (his[800+chan][9] > 100)
            s5 = (double) his[800+chan][10] / (double) his[800+chan][9];

        if (s3 < 100 || s1 < 100) continue;
        e0 /= s0*s0;
        e1 /= s1*s1; e2 /= s2*s2;
        e3 /= s3*s3; e4 /= s4*s4;
        e0 = sqrt(e0 + e1) * s0 / s1;
        e1 = sqrt(e1 + e2) * s2 / s1;
        e3 = sqrt(e3 + e4) * s4 / s3;
        e5 = sqrt(1.0/(double) his[800+chan][10] + 1.0/(double) his[800+chan][9]) * s5;
        printf("chan %3d   SEP, DEP, continuum acceptance:  %.4f(%.4f)  %.4f(%.4f)  %.4f(%.4f)\n",
               chan, s4/s3, e3, s2/s1, e0, s5, e5);
        fprintf(fp, "%4d %6.2f %5.2f %6.2f %5.2f %6.2f %5.2f\n",
                chan, s4/s3*100.0, e3*100.0, s2/s1*100.0, e0*100.0, s5*100.0, e5*100.0);
        s6 += s4/s3;
        s7 += s2/s1;
        s8 += s5;
        s9++;
      }
      if (s7 > 0 && s9 > 1) {
        printf("\n     >>> Mean acceptances:  %.4f  %.4f  %.4f\n\n", s6/s9, s7/s9, s8/s9);
        fprintf(fp, "#means: %5.2f %12.2f %12.2f\n",s6/s9*100.0, s7/s9*100.0, s8/s9*100.0);
      }
      fclose(fp);
    }
  }

  printf(">>  All done...\n");

  // write PSA data to psa.output
  PSA_info_write(&runInfo, &PSA);

  // write out histograms
  f_out = fopen("psa.rms", "w");
  for (i=0; i<HIS_COUNT; i++) {
   char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d A/; raw, norm, DT-corrected, DT- and E-correct; run %d", i, i, runInfo.runNumber);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d A/E DT-corrected with test factors", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d A/E_DT E-corrected with test factors", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d A/E-gated energy", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d Energy-gated A/E (BG, DEP, BG, SEP", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d DCR DT-corrected with test factors", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d lamda DT-corrected with test factors", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d DCR and lamda: raw, DT-corrected; shifted to mean or cut", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; DCR and lamda, unshifted", i);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d averaged DCR and lamda vs energy", i, i%200);
    } else if (i < 2200) {
      sprintf(spname, "%d; ch %d energy count for sp. 1800+chan", i, i%200);
    } else if (i < 2400) {
      sprintf(spname, "%d; ch %d averaged and smoothed (over 7 bins) DCR & lamda vs E", i, i%200);
    } else {
      sprintf(spname, "%d;", i);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  /* write spectra 1900 and 2200 for DCR and lamda vs. E correction  */
  f_out = fopen("INL_DCR_output.sec", "w");
  fwrite(his[2200], 8192*sizeof(int), 200, f_out);
  fwrite(his[1800], 8192*sizeof(int), 200, f_out);
  fclose(f_out);
  if (mean_dcr_ready  < 0)
    printf("\n Warning; no INL_DCR_input.sec was read; used internally generated values.\n"
           " Maybe mv INL_DCR_output.sec INL_DCR_input.sec ?\n\n"); 

  return 0;
}

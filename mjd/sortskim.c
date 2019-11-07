#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define MAKE_2D 1        // make a file with A/E, drift, and energy data for channel CHAN_2D, for 2d plots
#define CHAN_2D 0        // channel of interest for 2d plotting
#define HIS_COUNT 2000   // number of spectra in his[][] array

#define SUBTR_DCR_MEAN (e_adc < 8000 && \
                        ((chan < 100 && SUBTRACT_MEAN_DCR_HG) || \
                         (chan > 99  && SUBTRACT_MEAN_DCR_LG)))

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1;


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname_in\n", argv[0]);
    return -1;
  }
  /* open skim data file as input */
  while (argn < argc && argv[argn][0] == '-') argn++;
 
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
  float    aovere, aovere_norm, drift, dcr, lamda, s1, s2, s3;
  double   e_raw;
  int      nsd = 0, isd = 0;  // number of saved data, and pointer to saved data id

  double e_adc, e_ctc, e_lamda, e_fin, gain;
  int    i, j, ie;
  int    *his[HIS_COUNT];
  int    *dcr_mean[200], mean_dcr_ready = 0;
  FILE   *f_out, *f_out_2d = 0, *fp;

  float  elo_2d = 2605.0, ehi_2d = 2625.0;


  /* initialize */
  /* malloc and clear histogram space */
  if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
    printf("ERROR in sortskim.c; cannot malloc his!\n");
    exit(-1);
  }
  for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1] + 8192;
  if ((dcr_mean[0] = calloc(200*8192, sizeof(int))) == NULL) {
    printf("ERROR in sortskim.c; cannot malloc dcr_mean!\n");
    exit(-1);
  }
  for (i=1; i<200; i++) dcr_mean[i] = dcr_mean[i-1] + 8192;

  // read saved skim data from f_in
  fread(&nsd, sizeof(int), 1, f_in);
  fread(&Dets[0], sizeof(Dets[0]), NMJDETS, f_in);
  fread(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_in);
  if (runInfo.idNum == 0) {
    runInfo.flashcam = 1;
    fread(&(runInfo.flashcam), 8*sizeof(int), 1, f_in);
  }
  /* malloc space for SavedData */
  if ((sd = malloc(nsd*sizeof(*sd))) == NULL ||
      (sd[0] = malloc(nsd*sizeof(SavedData))) == NULL) {
    printf("ERROR in sortskim.c; cannot malloc SavedData!\n");
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

  for (isd = 0; isd < nsd; isd++) {
    // if (isd%(nsd/10) == 0) printf(">>  event %7d (%d/10\n", isd, isd*10/nsd);
    chan   = sd[isd]->chan;
    e_raw  = sd[isd]->e;
    drift  = sd[isd]->drift;
    aovere = sd[isd]->a_over_e;
    dcr    = sd[isd]->dcr;
    lamda  = sd[isd]->lamda;
    // if (sd[isd]->t100 - sd[isd]->t90 > 15) continue;
 
    if (chan < 100) {
      gain = Dets[chan].HGcalib[0];
    } else {
      gain = Dets[chan-100].LGcalib[0];
    }

    e_adc = e_raw + drift * CTC.e_dt_slope[chan];
    e_ctc = e_adc * gain;
    e_lamda = (e_raw + lamda * CTC.e_lamda_slope[chan] / 6000.0) * gain * CTC.e_lamda_gain[chan];
    if (CTC.best_dt_lamda[chan]) {
      e_fin = e_lamda;
    } else {
      e_fin = e_ctc;
    }
    ie = e_fin * 2.0 + 0.5;

    /* histogram raw A/E, E_dtc, E_lamda, final E, all with no psa cuts */
    if (e_raw < 8192)                              his[chan][(int)(e_raw + 0.5)]++;
    if (ie < 8192 && ie > 0)                       his[200+chan][ie]++;
    if ((j = e_ctc   * 2.0 + 0.5) < 8192 && j > 0) his[1400+chan][j]++;
    if ((j = e_lamda * 2.0 + 0.5) < 8192 && j > 0) his[1600+chan][j]++;

    /*  now calculate psa variables */
    lamda *= 8.0;      // scaled to roughly match DCR at 2 MeV

    float dtc = drift*3.0*2614.0/e_adc;   // in (x)us, for DT- correction to A/E
    if (chan%100 == 50) dtc *= 3.0;       // special hack for detector 50; has especially strong variation
    float ec = e_ctc/1000.0 - 2.0;        // in MeV, for E-correction to A/E

    float aovere2 = aovere + PSA.ae_dt_slope[chan] * dtc + PSA.ae_e_slope[chan] * ec;
    aovere_norm = aovere2 * 800.0/PSA.ae_pos[chan];
    if ((j = aovere_norm + 0.5) < 2000 && j > 10) his[1800+chan][j]++;

    s2 = aovere2;
    /*
    // adjust for energy dependence of cut
    // first correct for series noise: assume 2*sigma cut
    // series noise contribution to A/E sigma = BL_RMS * sqrt(2*rise) * factor / E_raw
    if (AOE_CORRECT_NOISE)
      s2 -= (2.0 * PZI.bl_rms[chan] * sqrt(2.0 * (float) PSA.a_e_rise[chan]) *
             PSA.a_e_factor[chan] * gain/1593.0 * (1.0 - 1593.0/e_ctc));        // 1593 keV = DEP
    // now deal with energy dependence of variation in A/E due to bremsstrahlung etc

    // linear; FIXME: Add limit at low e_ctc?
    // s2 -= AOE_CORRECT_EDEP * (PSA.ae_pos[chan] - PSA.ae_cut[chan]) * (1.0 - e_ctc/1593.0);
    // quadratic; FIXME: Add limit at low e_ctc?
    s2 -= AOE_CORRECT_EDEP * (PSA.ae_pos[chan] - PSA.ae_cut[chan]) * (1.0 - e_ctc*e_ctc/1593.0/1593.0);
    */
    s2 -= PSA.ae_cut[chan];
    if (PSA.ae_pos[chan] < 100) {
      s2 = -1000;
    } else {
      s2 *= 800.0/PSA.ae_pos[chan];
    }

    s1 = dcr - PSA.dcr_dt_slope[chan] * drift - PSA.dcr_lim[chan];
    if (SUBTR_DCR_MEAN) s1 -= (float) dcr_mean[chan][(int) e_adc/2] / 10.0;   // correct for residual INL

    s3 = lamda - PSA.lamda_dt_slope[chan] * drift - PSA.lamda_lim[chan];
    if (SUBTR_DCR_MEAN) s3 -= (float) dcr_mean[chan][4000 + (int) e_adc/2] / 10.0; // correct for residual INL

    if ((j = lrintf(s2)) > -500 && j < 500) his[1800+chan][3000+j]++;
    if ((j = lrintf(s1)) > -500 && j < 500) his[1800+chan][5000+j]++;
    if ((j = lrintf(s3)) > -500 && j < 500) his[1800+chan][6000+j]++;

    /* do A/E and DCR cuts on energy */
    if (ie < 8192 && ie > 0) {
      if (s2 >= 0)            his[400+chan][ie]++;  // A/E
      if (s1 <= 0)            his[600+chan][ie]++;  // DCR
      if (s2 >= 0 && s1 <= 0) his[800+chan][ie]++;  // A/E + DCR
      if (s3 <=0)             his[1000+chan][ie]++; // lamda
      if (s2 >= 0 && s3 <= 0) his[1200+chan][ie]++; // A/E + lamda
    }

    // make file for 2D plots  of A/E|E|CTC
    if (f_out_2d && chan == CHAN_2D && e_ctc >= elo_2d && e_ctc <= ehi_2d)
      fprintf(f_out_2d, "%4d %9.3f %8.2f %8.2f %7.2f %7.2f %7.2f %7.2f %6d\n",
              chan, e_ctc, aovere_norm, s2, drift, dtc,
              s1, s3, sd[isd]->t100 - sd[isd]->t90);

  }


  printf(">>  All done...\n");

  // write out histograms
  f_out = fopen("sort.rms", "w");
  for (i=0; i<HIS_COUNT; i++) {
   char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; ch %d Raw energy [ADC]; run %d", i, i, runInfo.runNumber);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d CT-corrected energy, no cuts [0.5 keV]", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d CT-corrected energy, A/E cut [0.5 keV]", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d CT-corrected energy, DCR cut [0.5 keV]", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d CT-corrected energy, A/E + DCR cuts [0.5 keV]", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d CT-corrected energy, lamda cut [0.5 keV]", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d CT-corrected energy, A/E + lamda cuts [0.5 keV]", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d Energy, DT-corrected [0.5 keV]", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; ch %d Energy, lamda-corrected [0.5 keV]", i, i%200);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d A/E raw, moved to cut; DCR, lamda moved to cut", i, i%200);
    } else {
      sprintf(spname, "%d;", i);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  if (mean_dcr_ready  < 0)
    printf("\n WARNING; no INL_DCR_input.sec was read; maybe mv INL_DCR_output.sec INL_DCR_input.sec ?\n\n"); 
  return 0;
}

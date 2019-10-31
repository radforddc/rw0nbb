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

int fitter(float *pars, float *sig, int tlo, int thi);

int main(int argc, char **argv) {

  MJDetInfo  Dets[NMJDETS];
  MJRunInfo  runInfo;
  int        argn=1;
  char       *c, data_file_name[256], list_file_name[256];
  FILE       *f_in = NULL, *f_file_list = NULL;


  if (argc < 2) {
    fprintf(stderr, "\nusage: %s fname_in [chnum_lo] [chnum_hi] [e_lo] [e_hi] [-n]\n\n", argv[0]);
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
      return 0;
    }
    if (!fgets(data_file_name, sizeof(data_file_name), f_file_list)) {
      fprintf(stderr, "\n Failed to read input file %s\n", list_file_name);
      return 0;
    }
    for (c = data_file_name + strlen(data_file_name); *c < '0'; c--) *c = 0;  // removing trailing junk
  } else {
    strncpy(data_file_name, argv[argn], sizeof(data_file_name));
  }

  f_in = fopen(data_file_name, "r");
  if (f_in == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", data_file_name);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", data_file_name);
  strncpy(runInfo.filename, data_file_name, 256);
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* read file header */
  int nDets = decode_runfile_header(f_in, Dets, &runInfo);
  if (nDets < 1) return 1;

  if (!runInfo.flashcam) {
    if (runInfo.dataIdGM == 0 && runInfo.dataIdGA == 0)
      printf("\n No data ID found for Gretina4M or 4A data!\n");
    if (runInfo.dataIdGM)
      printf("\n Data ID %d found for Gretina4M data\n", runInfo.dataIdGM);
    if (runInfo.dataIdGA)
      printf("\n Data ID %d found for Gretina4A data\n", runInfo.dataIdGA);
    printf(" Run number: %d in file %s\n", runInfo.runNumber, runInfo.filename);
  }

/* ---------------------------------------------------- */

  int totevts=0, out_evts=0;
  int clo=0, chi, elo=3000, ehi=7200;
  int module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  int det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  int chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  int presum[200] = {0};         // expected presum step (or zero for no presum)
  PZinfo  PZI;
  CTCinfo CTC;
  PSAinfo PSA;

  unsigned int head[2], evtdat[20000];
  short  *signal, sigu[8192], siguc[8192];
  float  fsignal[8192] = {0};

  // data used, stored, and reused in the different steps
  SavedData **sd;
  int    chan;
  double e_raw;
  float  drift, aovere, dcr, lamda;
  int    nsd = 0, isd = 0;  // number of saved data, and pointer to saved data id
  int    sdchunk = 1000000; // number of SavedData events to malloc at one time

  double e_ctc, e_adc;
  int    i, j, k, t90, t100, sig_len, his[8192] = {0};
  FILE   *f_out;


  /* initialize */

  /* malloc initial space for SavedData */
  if ((sd = malloc(sdchunk*sizeof(*sd))) == NULL ||
      (sd[0] = malloc(sdchunk*sizeof(SavedData))) == NULL) {
    printf("ERROR in skim.c; cannot malloc SavedData!\n");
    exit(-1);
  }
  for (i=1; i<sdchunk; i++) sd[i] = sd[i-1] + 1;
  nsd = sdchunk;

  if (ep_init(Dets, &runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
  /* read pole-zero values from PZ.input */
  if (!PZ_info_read(&runInfo, &PZI)) {
    printf("\n ERROR: No initial pole-zero data read. Does PZ.input exist?\n");
    return -1;
  }
  /* read energy correction factors from ctc.input */
  if (!CTC_info_read(&runInfo, &CTC)) {
    printf("\n Warning: No initial charge-trapping correction data read. Does ctc.input exist?\n");
  }
  /* read individual trapezoid values from filters.input (if it exists) */
  if (2 == (PSA_info_read(&runInfo, &PSA))) {
    printf("\n Individual trap values read from filters.input\n");
  }

  // see if channel and energy limits are defined in the command line
  chi=100+runInfo.nGe-1;
  elo = 3000;
  ehi = 7200;
  j = 2;

  if (runInfo.argc > j) clo = atoi(runInfo.argv[j++]);
  if (runInfo.argc > j) chi = atoi(runInfo.argv[j++]);
  if (runInfo.argc > j) elo = atoi(runInfo.argv[j++]);
  if (runInfo.argc > j) ehi = atoi(runInfo.argv[j++]);
  if (clo < 0) clo = 0;
  if (chi > 100+runInfo.nGe) chi = 100+runInfo.nGe;

  printf("\nChs %d to %d, e_trapmax %d to %d\n\n", clo, chi, elo, ehi);

  /* decide if there is presumming; if so, then we also need to
     make sure we have the space to hold expanded signals */
  k = 2008;
  int chan_k = -1;
  for (chan = 0; chan < runInfo.nGe; chan++) {
    if (runInfo.flashcam) continue;
    // HG channels
    if (Dets[chan].type == 0 && Dets[chan].HGChEnabled &&           // GRETINA4M type
        (Dets[chan].HGPostrecnt + Dets[chan].HGPrerecnt) < 2008) {
      presum[chan] = Dets[chan].HGPostrecnt + Dets[chan].HGPrerecnt;
      if (k > presum[chan]) {
        k = presum[chan];
        chan_k = chan;
      }
    }
    if (Dets[chan].type == 1 && Dets[chan].HGPreSumEnabled) {       // GRETINA4A type
      presum[chan] = 1 << Dets[chan].HGdecimationFactor;
      k = 1111;
      chan_k = chan;
    }
    if (VERBOSE > 1 && presum[chan] > 0)
      printf("chan, type, presum = %d %d %d\n", chan, Dets[chan].type, presum[chan]);

    // LG channels
    if (Dets[chan].type == 0 && Dets[chan].LGChEnabled &&           // GRETINA4M type
        (Dets[chan].LGPostrecnt + Dets[chan].LGPrerecnt) < 2008) {
      presum[100+chan] = Dets[chan].LGPostrecnt + Dets[chan].LGPrerecnt;
      if (k > presum[100+chan]) {
        k = presum[100+chan];
        chan_k = 100+chan;
      }
    }
    if (Dets[chan].type == 1 && Dets[chan].LGPreSumEnabled) {       // GRETINA4A type
      presum[100+chan] = 1 << Dets[chan].LGdecimationFactor;
      k = 1111;
      chan_k = 100+chan;
    }
    if (VERBOSE > 1 && presum[100+chan] > 0)
      printf("chan, type, presum = %d %d %d\n", 100+chan, Dets[chan].type, presum[100+chan]);
  }
  if (k < 2008 && k > 1000) {  // yes, at least one channel has presumming enabled
    printf("   Presumming: At least one channel has PS start at or after sample %d\n", k);
    if (k == 1111) { // GRETINA4A type, x16 expansion
      k = presum[chan_k]*2048;   // TEST TEMPORARY
    } else {
      k = 2048 + (2018-k)*3;     // max expanded signal length
    }
    if (k > 8192) {
      printf("ERROR: Expanded signals will be longer than 8192 samples! Code edit required!\n");
      exit(-1);
    }
  } else {
    printf("No presumming\n");
  }


  // end of initialization
  // start loop over reading events from input file

  while (1) {
    int evlen = 0, crate=0, slot=0, board_type = 0, last_read = 0;
    long long int time = 0;
    chan = -1;

    if (runInfo.flashcam) {  // ----- flashcam data
      while ((last_read = fread(&k, sizeof(int), 1, f_in)) == 1 && k < 1200) {
        if (k > 1) {
          if ((last_read = fread(evtdat, k, 1, f_in)) != 1) break;
        }
      }
      evlen = k/4 + 2;
      chan = 0;

    } else {  // ----- not flashcam data
      if ((last_read = fread(head, sizeof(head), 1, f_in)) == 1) {
        board_type = head[0] >> 18;
        evlen = (head[0] & 0x3ffff);

        if (board_type == 0) {  // a new runfile header! file must be corrupt?
          printf("\n >>>> ERROR: DataID = 0; found a second file header??"
                 " Ending scan of this file!\n"
                 " >>>> head = %8.8x %8.8x  evlen = %d\n", head[0], head[1], evlen);
          break;
        }

        /* if we don't want to decode this type of data, just skip forward in the file */
        if (board_type != runInfo.dataIdGM &&
            board_type != runInfo.dataIdGA) {
          if (evlen > 10000) {
            printf("\n >>>> ERROR: Event length too long??  board_type = %d, length = %d\n"
                   " >>>> This file is probably corruped, ending scan!\n", board_type, evlen);
            break;
          }
          fseek(f_in, 4*(evlen-2), SEEK_CUR);
          continue;
        }

        slot  = (head[1] >> 16) & 0x1f;
        crate = (head[1] >> 21) & 0xf;
        if (crate < 0 || crate > NCRATES ||
            slot  < 0 || slot > 20) {
          printf("ERROR: Illegal VME crate or slot number %d %d\n", crate, slot);
          if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) break;
          continue;
        }
      }
    }   // -------- if (flashcam) {} else {

    /* ========== read in the rest of the event data ========== */
    if (last_read != 1 || fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
      /* if no more data, *break* to go to histgram processing
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
      printf(" %8d evts in, %d out, %d saved\n", totevts, out_evts, isd); fflush(stdout);
    }

    if (runInfo.flashcam) {  // ----- flashcam data

      if (chan < clo || chan > chi) continue;
      sig_len = 2*(evlen-2);
      signal = (short *) evtdat;

    } else {  // ----- not flashcam data

      int ch = (evtdat[1] & 0xf);
      if ((j = module_lu[crate][slot]) < 0 || ch >= 10) continue;
      chan = chan_lu[j][ch];
      if (chan > 99 + runInfo.nGe && chan < 100 + runInfo.nGe + runInfo.nPT) continue; // pulser tag channels
      if (chan < 0 || chan > 99 + runInfo.nGe + runInfo.nPT ||
          ((chan < 100 && !Dets[chan].HGChEnabled) ||
           (chan > 99 && !Dets[chan-100].LGChEnabled))) {
        printf("Data from detector not enabled! Chan = %d  crate, slot, j, ch = %d %d %d %d  len = %d\n",
               chan, crate, slot, module_lu[crate][slot], ch, evlen);
        continue;
      }

      if (++totevts % 50000 == 0) {   
        printf(" %8d evts in, %d out, %d saved\n", totevts, out_evts, isd); fflush(stdout);
      }
      if (chan < clo || chan > chi) continue;

      /* --------------- Gretina4M or 4A digitizer data ---------------- */
      time = (evtdat[3] & 0xffff);
      time = time << 32 | evtdat[2];
      if (time < 0) continue;
      sig_len = 2008;

      signal = (short *) evtdat + 28;
      if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
          2020 == decompress_signal((unsigned short *)signal, sigu, 2*(evlen - 2) - 28)) {
        signal = sigu;
        evlen = 1026;
      }

      /* deal with presumming of signal */
      if (presum[chan] > 0) {   // need to correct for presumming
        // FIXME: replace hard-coded factor of 4 with data from header
        int step = presum[chan];
        int cts = (signal[step-5] + signal[step-4] + signal[step-3] + signal[step-2])/4;
        if (cts > 10 || cts < -10) {
          /* --- Find transition to 4x presumming ------- */
          for (step=presum[chan]; step<presum[chan]+2; step++) {
            if (cts > 0 && signal[step] > 2*cts) break;
            if (cts < 0 && signal[step] < 2*cts) break;
          }
          if (step == presum[chan]+2) {
            printf("Hmmm... step not found in %d - %d in chan %3d, counts %5d; discarding event\n",
                   step-1, step, chan, cts);
            continue;
          } else if (VERBOSE > 1) {
            printf("step found at %d in chan %3d\n", step, chan);
          }
        }
        /* --- uncompress presumming --- */
        for (i=0; i<step; i++) siguc[i] = signal[i];
        for (i=step; i<sig_len; i++) {                   // NOTE: This is wrong for signal < 0
          for (j=0; j<4; j++)
            siguc[step + 4*(i-step) + j] = signal[i]/4; // distribute sums over 4 bins each
          for (j=0; j<signal[i]%4; j++)
            siguc[step + 4*(i-step) + j]++;             // and make sure sum is right
        }
        sig_len += 3*(sig_len-step);
        signal = siguc;
        /* --- Done with 4x presumming ------- */
        if (VERBOSE > 1)
          printf("Corrected chan %d for presumming at step %d, sig_len = %d\n",
                 chan, step, sig_len);
      }
    }   // -------- if (flashcam) {} else {

    int e = trap_max(signal, &j, TRAP_RISE, TRAP_FLAT)/TRAP_RISE;
    if (runInfo.flashcam) e /= 2;
    if (chan < 100 && (e < elo || e > ehi)) continue;
    if (chan > 99 && (e < elo/3.4 || e > ehi/3.2)) continue;
    out_evts++;

    /* sticky-bit fix */
    int d = 128;  // basically a sensitivity threshold; max change found will be d/2
    if (chan > 99 && chan < 100+runInfo.nGe) d = 64;
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
        printf("Fixing sticky bit in signal %d, chan %d, t=%d, change %d\n",
               out_evts-1, chan, i, j*d/2);
        signal[i] += j*d/2;
        // break;
      }
    }

    if (sig_len > 2500) sig_len = 2500;   // FIXME

    /* ---------------------- process selected signals ---------------------- */
      
    /* find t100 and t90*/
    t100 = 700;                 // FIXME? arbitrary 700?
    for (i = t100+1; i < 1500; i++)
      if (signal[t100] < signal[i]) t100 = i;
    if (t100 > 1300) continue;  // FIXME ??  - important for cleaning, gets rid of pileup
    /* get mean baseline value */
    int bl = 0;
    for (i=300; i<400; i++) bl += signal[i];
    bl /= 100;
    if ((bl - PZI.baseline[chan]) < -10 || (bl - PZI.baseline[chan]) > 50) continue;   // a little data cleaning
    for (t90 = t100-1; t90 > 500; t90--)
      if ((signal[t90] - bl) <= (signal[t100] - bl)*19/20) break;

    /* do (optional) INL  correction */
    if (DO_INL) {
      if (inl_correct(Dets, &runInfo, signal+10, fsignal+10, sig_len-10, chan)) {
        printf(" >>> inl_correct return error for chan %d\n", chan);
        return -1;
      }
    } else {
      for (i=0; i<sig_len; i++) fsignal[i] = signal[i];
    }

    /* do fitting of pole-zero parameters to get lamda (~ DCR) */

    float chisq, lamda1, frac2;
    int   tlo = t100+50, thi = sig_len - 10;
    if (thi > tlo + 1500) thi = tlo + 1500;  // FIXME; check performance
    chisq = pz_fitter(fsignal, tlo, thi, chan, &PZI, &lamda1, &frac2, &lamda);
    if (chisq < 0.01 || chisq > 10.0) continue;            // fit failed, or bad fit chisq
    lamda = (0.01 / PZI.tau[chan] - lamda) * 1.5e6 + 0.3;  // 0.3 fudge to get mean ~ 0

    /* do (required) PZ correction */
    if (PZ_fcorrect(fsignal+10, sig_len-10, chan, &PZI))
      printf(" >>> PZ_fcorrect return error for chan %d\n", chan);

    //-------------------------

    /* get raw (e_raw) and drift-time-corrected energy (e_adc and e_ctc)
       and effective drift time */
    int tmax, t0;
    e_ctc = get_CTC_energy(fsignal, sig_len, chan, Dets, &PSA,
                           &t0, &e_adc, &e_raw, &drift, CTC.e_dt_slope[chan]);
    if (e_ctc < 0.001) {
      printf("E_ctc = %.1f < 1 eV!\n", e_ctc);
      if (VERBOSE) printf("chan, t0, t100: %d %d %d; e, drift: %.2f %.2f\n",
                          chan, t0, t100, e_raw, drift);
      continue;
    }
    if (runInfo.flashcam) {
      e_ctc /= 2.0;
      e_adc /= 2.0;
      e_raw /= 2.0;
      drift /= 2.0;
    }

    if (elo == 3000 && (E_THRESH > 0 && e_ctc < E_THRESH)) continue;

    /* find A/E */
    if (e_ctc > 50 && t0 > 600 && t0 < 1300) {
      aovere = float_trap_max_range(fsignal, &tmax, PSA.a_e_rise[chan], 0, t0-20, t0+300);
      aovere *= PSA.a_e_factor[chan] / e_raw;
      if (runInfo.flashcam) aovere /= 2.0;
    } else {
      aovere = 0;
    }

    /* ---- This next section calculates the GERDA-style A/E ---- */
    if (PSA.gerda_aoe[chan]) {
      float  ssig[6][2000];
      for (k=100; k<2000; k++) ssig[0][k] = fsignal[k] - fsignal[300];
      for (j=1; j<6; j++) {  // number of cycles
        for (i=200; i < 1500; i++) {
          ssig[j][i] = 0.0;
          for (k=0; k<10; k++) {
            ssig[j][i] += ssig[j-1][i+k-5];
          }
          ssig[j][i] /= 10.0;
        }
      }
      aovere = 0;
      for (k=250; k<1450; k++) {
        if (aovere < ssig[5][k] - ssig[5][k-1]) aovere = ssig[5][k] - ssig[5][k-1];
      }
      aovere *= 12.71*1593/e_raw;
    }
    /* ---- end of GERDA-style A/E ---- */

    /* get DCR from slope of PZ-corrected signal tail */
    if (sig_len < 2450) {
      if (t90 + 760 > sig_len) {
        if (VERBOSE)
          printf("Error getting DCR in skim.c: t90 = %d, sig_len = %d\n", t90, sig_len);
        continue;
      }
      dcr = float_trap_fixed(fsignal, t90+50, 100, 500) / 25.0;
    } else {
      if (t90 + 1180 > sig_len) {
        if (VERBOSE)
          printf("Error getting DCR in skim.c: t90 = %d, sig_len = %d\n", t90, sig_len);
        continue;
      }
      dcr = float_trap_fixed(fsignal, t90+50, 160, 800) / 40.0;
    }
    if (runInfo.flashcam) dcr /= 2.0;

    if (e_raw > 0 && e_raw < 8000) his[(int) e_raw]++;
    // save data for skim file
    sd[isd]->chan     = chan;
    sd[isd]->e        = e_raw;
    sd[isd]->drift    = drift;
    sd[isd]->a_over_e = aovere;
    sd[isd]->dcr      = dcr;
    sd[isd]->lamda    = lamda;
    sd[isd]->t0       = t0;
    sd[isd]->t90      = t90;
    sd[isd]->t100     = t100;
    if (++isd == nsd) { // space allocated for saved data is now full; malloc more
      if ((sd = realloc(sd, (isd + sdchunk)*sizeof(*sd))) == NULL ||
          (sd[isd] = malloc(sdchunk*sizeof(SavedData))) == NULL) {
        printf("ERROR in skim.c; cannot realloc SavedData! nsd = %d\n", nsd);
        exit(-1);
      }
      nsd += sdchunk;
      for (i=isd+1; i<nsd; i++) sd[i] = sd[i-1]+1;
    }
  }
  f_out = fopen("skim.sec", "w");
  fwrite(his, sizeof(int), 8192, f_out);
  fclose(f_out);

  printf(" %8d evts in, %d out, %d saved\n", totevts, out_evts, isd);
  nsd = isd;

  // save skim data (SavedData) to disk
  f_out = fopen("skim.dat", "w");
  fwrite(&nsd, sizeof(int), 1, f_out);
  fwrite(&Dets[0], sizeof(Dets[0]), NMJDETS, f_out);
  if (runInfo.flashcam) {
    runInfo.idNum = 0;
    fwrite(&runInfo, sizeof(runInfo), 1, f_out);  // for backwards compatibility
  } else {
    fwrite(&runInfo, sizeof(runInfo) - 8*sizeof(int), 1, f_out);
  }
  // fwrite(*sd, sizeof(SavedData), nsd, f_out);
  // writing large files all in one go doesn't seem to work?
  for (isd = 0; isd < nsd; isd += sdchunk) {
    if (nsd - isd < sdchunk) {
      fwrite(sd[isd], sizeof(SavedData), nsd - isd, f_out);
    } else {
      fwrite(sd[isd], sizeof(SavedData), sdchunk, f_out);
    }
  }
  fclose(f_out);
  printf("\n Wrote %d skimmed events to skim.dat\n\n", nsd);

  fclose(f_in);
  printf("\n All Done.\n\n");
  return 0;
}

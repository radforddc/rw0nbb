#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define HIS_COUNT 400  // number of spectra in his[][] array

/*  ------------------------------------------------------------ */

void signalselect(FILE *f_in, MJDetInfo *detInfo, MJRunInfo *runInfo, int step);

int main(int argc, char **argv) {

  MJDetInfo  detInfo[NMJDETS];
  MJRunInfo  runInfo;
  int        argn = 1, use_fraction = 0;


  if (argc < 2) {
    fprintf(stderr,
            "\nusage: %s [-a n] fname_in [chnum_lo] [chnum_hi] [e_lo] [e_hi]\n"
            "         option -a: Use only one out of every n events, for faster processing of large datasets\n\n",
            argv[0]);
    return -1;
  }
  if (argc > 3 && !strncmp(argv[1], "-a", 2)) {
    use_fraction = atoi(argv[2]);
    if (use_fraction < 1 || use_fraction > 100) {
      fprintf(stderr, "ERROR in use_fraction = %s; aborting.\n\n", argv[2]);
      return 0;
    }
  }
  /* open raw/presorted data file as input */
  while (argn < argc && argv[argn][0] == '-') argn += 2;
  FILE *f_in = fopen(argv[argn],"r");
  if (f_in == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", argv[argn]);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", argv[argn]);
  strncpy(runInfo.filename, argv[argn], 256);
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* read file header */
  int nDets = decode_runfile_header(f_in, detInfo, &runInfo);
  if (nDets < 1) return 1;

  if (runInfo.flashcam) {
    printf("FlashCam data; assuming only one detector.\n");
    if (runInfo.flashcam < 2)
      printf("\nNOTE: I suggest you first run preflash on (a list of) your files.\n\n");
  } else {
    if (runInfo.dataIdGM == 0 && runInfo.dataIdGA == 0)
      printf("\n No data ID found for Gretina4M or 4A data!\n");
    if (runInfo.dataIdGM)
      printf("\n Data ID %d found for Gretina4M data\n", runInfo.dataIdGM);
    if (runInfo.dataIdGA)
      printf("\n Data ID %d found for Gretina4A data\n", runInfo.dataIdGA);
    printf(" Run number: %d in file %s\n", runInfo.runNumber, runInfo.filename);
  }

  /* read through all the events in the file, to build and process them */
  runInfo.analysisPass = 0;
  signalselect(f_in, detInfo, &runInfo, 1);

  runInfo.analysisPass = 1;
  fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);  // go back to start of data in file
  signalselect(f_in, detInfo, &runInfo, 2);

  fclose(f_in);
  printf("\n All Done.\n\n");
  return 0;
}

/* ========================================================== */

void signalselect(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo, int step) {

  int totevts=0, out_evts=0;
  static int module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs

  unsigned int head[2], evtdat[20000];
  short  *signal, sigu[8192], siguc[8192];
  float  fsignal[8192] = {0};

  int    i, j, k, chan, sig_len, use_fraction = 0, argn=2;
  FILE   *f_out;
  //FILE   *f_mat;
  int    t95, t100, bl;
  float  pos, area, fwhm;

  static int    *his[HIS_COUNT];
  static int    clo=0, chi, elo=3000, ehi=7200;
  static int    presum[200] = {0};     // expected presum step (or zero for no presum)
  static PZinfo PZI;


  if (runInfo->argc > 3 && !strncmp(runInfo->argv[1], "-a", 2)) {
    use_fraction = atoi(runInfo->argv[2]);
    printf("\n Using one of every %d events...\n", use_fraction);
    argn += 2;
  }

  /*
    step = 1 to get baseline mode
    step = 2 to get PZ tau and frac2
  */
  printf(" >>> Step %d of 2 <<<\n", step); fflush(stdout);
  if (step == 1) {
    //f_mat = fopen("s.mat", "w");
    /* initialize */
    /* malloc and clear histogram space */
    if ((his[0] = calloc(HIS_COUNT*8192, sizeof(int))) == NULL) {
      printf("ERROR in PZcal.c signalselect(); cannot malloc his!\n");
      exit(-1);
    }
    for (i=1; i<HIS_COUNT; i++) his[i] = his[i-1]+8192;

    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return;
    if (!(j = PZ_info_read(runInfo, &PZI))) {
      printf("\n WARNING: No initial pole-zero data read.\n"
             " >>>> Starting with first baseline seen for each channl as initial estimate.\n");
      for (chan = 0; chan < 200; chan++) {
        PZI.baseline[chan] = -9999;  // special flag to say that baseline is unknown
      }
    }

    chi=100+runInfo->nGe-1;
    elo = 3000;
    ehi = 7200;
    if (runInfo->argc > argn) clo = atoi(runInfo->argv[argn++]);
    if (runInfo->argc > argn) chi = atoi(runInfo->argv[argn++]);
    if (runInfo->argc > argn) elo = atoi(runInfo->argv[argn++]);
    if (runInfo->argc > argn) ehi = atoi(runInfo->argv[argn++]);
    if (clo < 0) clo = 0;
    printf("\nChs %d to %d, e_trapmax %d to %d\n\n", clo, chi, elo, ehi);

    /* decide if there is presumming; if so, then we also need to
       make sure we have the space to hold expanded signals */
    k = 2008;
    int chan_k = -1;
    for (chan = 0; chan < runInfo->nGe; chan++) {
      if (runInfo->flashcam) continue;
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
      if (VERBOSE && presum[chan] > 0)
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
      if (VERBOSE && presum[100+chan] > 0)
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

  }

  /* start loop over reading events from input file
     ============================================== */

  while (1) {
    int evlen = 0, crate=0, slot=0, board_type = 0;
    long long int time = 0;
    chan = -1;

    if (runInfo->flashcam) {  // ----- flashcam data -----
      chan = 0;
      if (chan < clo || chan > chi) continue;

      if (runInfo->flashcam < 2) {  // not presorted
        int last_read;
        while ((last_read = fread(&k, sizeof(int), 1, f_in)) == 1 && k < 1200) {
          if (k > 1) {
            if (fread(evtdat, k, 1, f_in) != 1) break;
          }
        }
        evlen = k/4 + 2;
        /*  read in the rest of the event data  */
        if (last_read != 1 || fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
          printf("  No more data...\n");
          break;
        }
        sig_len = 2*(evlen-2);
        signal = (short *) evtdat;
        for (i=0; i < sig_len; i++) signal[i] = ((unsigned short) signal[i]) / 4 - 8192;
        // if rising edge comes too late in signal, trapmax will fail
        if (sig_len > 2400) {
          signal  += (sig_len-2400)/2;
          sig_len -= (sig_len-2400)/2;
        }
      } else {        // presorted data, contains only waveforms
        if (fread(sigu,  sizeof(short), 2, f_in) != 2 ||
            fread(siguc, sizeof(short), sigu[0], f_in) != sigu[0]) {
          printf("  No more data...\n");
          break;
        }
        i = sigu[0];            // compressed length
        j = sig_len = sigu[1];  // uncompressed length
        signal = siguc;
        if (i != j) {
          sig_len = decompress_signal((unsigned short *)siguc, sigu, i);
          signal = sigu;
        }
        if (j != sig_len) {
          printf("\nERROR decompressing signal; sig_len %d != %d\n", sig_len, j);
          return;
        }
      }

      if (++totevts % 50000 == 0)
        printf(" %8d evts in, %d out\n", totevts, out_evts); fflush(stdout);

    } else {  // ----- not flashcam data -----
      if (fread(head, sizeof(head), 1, f_in) != 1) break;
      board_type = head[0] >> 18;
      evlen = (head[0] & 0x3ffff);

      if (board_type == 0) {  // a new runfile header! file must be corrupt?
        printf("\n >>>> ERROR: DataID = 0; found a second file header??"
               " Ending scan of this file!\n"
               " >>>> head = %8.8x %8.8x  evlen = %d\n", head[0], head[1], evlen);
        break;
      }

      /* if we don't want to decode this type of data, just skip forward in the file */
      if (board_type != runInfo->dataIdGM &&
          board_type != runInfo->dataIdGA) {
        if (evlen > 10000) {
          printf("\n >>>> ERROR: Event length too long??\n"
                 " >>>> This file is probably corruped, ending scan!\n");
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

      /* ========== read in the rest of the event data ========== */
      if (fread(evtdat, sizeof(int), evlen-2, f_in) != evlen-2) {
        printf("  No more data...\n");
        break;
      }
      int ch = (evtdat[1] & 0xf);
      if ((j = module_lu[crate][slot]) >= 0 && ch < 10) chan = chan_lu[j][ch];
      if (chan > 99 + runInfo->nGe && chan < 100 + runInfo->nGe + runInfo->nPT) continue; // pulser tag channels
      if (chan < 0 || chan > 99 + runInfo->nGe + runInfo->nPT ||
          ((chan < 100 && !Dets[chan].HGChEnabled) ||
           (chan > 99 && !Dets[chan-100].LGChEnabled))) {
        printf("Data from detector not enabled! Chan = %d  crate, slot, j, ch = %d %d %d %d  len = %d\n",
               chan, crate, slot, module_lu[crate][slot], ch, evlen);
        continue;
      }

      if (++totevts % 50000 == 0) {
        printf(" %8d evts in, %d out\n", totevts, out_evts); fflush(stdout);
      }
      if (chan < clo || chan > chi) continue;
      if (use_fraction && totevts % use_fraction != 0) continue;

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
        int sstep = presum[chan];
        int cts = (signal[sstep-5] + signal[sstep-4] + signal[sstep-3] + signal[sstep-2])/4;
        if (cts > 10 || cts < -10) {
          /* --- Find transition to 4x presumming ------- */
          for (sstep=presum[chan]; sstep<presum[chan]+2; sstep++) {
            if (cts > 0 && signal[sstep] > 2*cts) break;
            if (cts < 0 && signal[sstep] < 2*cts) break;
          }
          if (sstep == presum[chan]+2) {
            printf("Hmmm... step not found in %d - %d in chan %3d, counts %5d; discarding event\n",
                   sstep-1, sstep, chan, cts);
            continue;
          } else if (VERBOSE) {
            printf("step found at %d in chan %3d\n", sstep, chan);
          }
        }
        /* --- uncompress presumming --- */
        for (i=0; i<sstep; i++) siguc[i] = signal[i];
        for (i=sstep; i<sig_len; i++) {                   // NOTE: This is wrong for signal < 0
          for (j=0; j<4; j++)
            siguc[sstep + 4*(i-sstep) + j] = signal[i]/4; // distribute sums over 4 bins each
          for (j=0; j<signal[i]%4; j++)
            siguc[sstep + 4*(i-sstep) + j]++;             // and make sure sum is right
        }
        sig_len += 3*(sig_len-sstep);
        signal = siguc;
        /* --- Done with 4x presumming ------- */
        if (VERBOSE)
          printf("Corrected chan %d for presumming at step %d, sig_len = %d\n",
                 chan, sstep, sig_len);
      }

      /* sticky-bit fix */
      int d = 128;  // basically a sensitivity threshold; max change found will be d/2
      if (chan > 99 && chan < 100+runInfo->nGe) d = 64;
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
    }   // -------- if (flashcam) {} else {

    int e = trap_max(signal, &j, TRAP_RISE, TRAP_FLAT)/TRAP_RISE;
    if (runInfo->flashcam) e /= 2;

    his[100][e]++;
    if (chan < 100 && (e < elo || e > ehi)) continue;
    if (chan > 99 && (e < elo/3.4 || e > ehi/3.2)) continue;
    his[110][e]++;
    out_evts++;


    /* ---------------------- process selected signals ---------------------- */
    if (sig_len > 3500) sig_len = 3500;   // FIXME

    if (step == 1) {
      /* just get mean and RMS for baseline value */
      double s1 = 0, s2 = 0;
      for (i=300; i<400; i++) {
        double s = signal[i];
        s1 += s;
        s2 += s*s;
      }
      s1 /= 100.0;                   // mean baseline
      s2 = sqrt(s2/100.0 - s1*s1);   // RMS
      /* histogram mean baseline value (for new PZ.output) */
      if (PZI.baseline[chan] == -9999) {  // baselines unknown! use this first signal from chan as estimate.
        PZI.baseline[chan] = s1;
        printf("Using BL of first signal (%5.0f) as BL guess for channel %3d\n", s1, chan);
      }
      bl = s1 + 1000.5 - PZI.baseline[chan];
      if (bl > 0 && bl < 1900) his[chan][bl]++;
      if (s2 < 50.0) his[200+chan][(int)(10.0*s2 + 0.5)]++;
      continue;
    }

    // now do step 2 stuff
    /* find and histogram t100 and t95 */
    t100 = 700;                 // FIXME? arbitrary 700?
    for (i = t100+1; i < sig_len - 500; i++)
      if (signal[t100] < signal[i]) t100 = i;
    if (t100 > sig_len - 700) continue;  // FIXME ??  - important for cleaning, gets rid of pileup
    /* get mean baseline value */
    bl = 0;
    for (i=300; i<400; i++) bl += signal[i];
    bl /= 100;
    if ((bl - PZI.baseline[chan]) < -10 || (bl - PZI.baseline[chan]) > 50) continue;   // a little data cleaning
    for (t95 = t100-1; t95 > 500; t95--)
      if ((signal[t95]-bl) <= (signal[t100] - bl)*19/20) break;
    his[200+chan][2000 + t95]++;
    his[200+chan][4000 + t100-t95]++;

    /* do (optional) INL  correction */
    if (DO_INL) {
      if (inl_correct(Dets, runInfo, signal+10, fsignal+10, sig_len-10, chan))
        printf(" >>> inl_correct return error for chan %d\n", chan);
    } else {
      for (i=10; i<sig_len; i++) fsignal[i] = signal[i];
    }
    //fwrite(signal, sizeof(short), sig_len, f_mat);

    /* do fitting of pole-zero parameters (decay time and fraction_2) */

    float chisq, lamda1, frac2, lamda;
    int   tlo = t100+50, thi = sig_len - 10;
    if (thi > tlo + 1500) thi = tlo + 1500;   // FIXME; check performance

    chisq = pz_fitter(fsignal, tlo, thi, chan, &PZI, &lamda1, &frac2, &lamda);
    if (chisq > 0.01) {       // fit successful
      if (chisq < 5.0) {      // good fit
        /* histogram fitted parameter values */
        j = lamda1 * 3000000.0;
        if (j > 0 && j < 1500)   his[chan][j + 2000]++;  // 0.01/tau
        j = frac2 * 4000.0;
        if (j > -100 && j < 400) his[chan][j + 4000]++;  // frac2
        j = lamda * 3000000.0;
        if (j > 0 && j < 1500)   his[chan][j + 4500]++;  // lamda
      }
      j = chisq * 100.0;
      if (j < 2000) {
        his[chan][j + 6000]++;
      } else {
        his[chan][1960]++;
        if (chisq > 100) his[chan][1970]++;
      }
    } else {                 // fit unsuccessful
      his[chan][1950]++;
    }
  }

  printf(" %8d evts in, %d out\n", totevts, out_evts);

  if (step == 1) {   // adjust only baseline for PZ
    // find new baseline values;
    printf("Chan  BL ->  BL ;  BL_RMS -> BL_RMS\n");
    //        1  116 -> 116 ;    2.15 ->   2.15
    for (chan=clo; chan<=chi; chan++) {
      j = 10;
      for (i=11; i<1900; i++)
        if (his[chan][j] < his[chan][i]) j = i;  // mode of baseline distribution
      if (his[chan][j] < 200) continue;
      printf("%3d %4.0f ->", chan, PZI.baseline[chan]);
      PZI.baseline[chan] += j-1000;
      printf("%4.0f ; %7.2f ->", PZI.baseline[chan], PZI.bl_rms[chan]);

      double s1 = 0, s2 = 0;
      for (i=0; i<1000; i++) {
        s1 += his[200+chan][i];
        s2 += i * his[200+chan][i];
      }
      if (s1 > 100) {
        s2 /= s1;
        PZI.bl_rms[chan] = s2/10.0;;
        printf(" %7.2f\n", PZI.bl_rms[chan]);
      } else {
        printf("\n");
      }
    }
    printf("\n ------------------------ End of step 1 \n\n");
    return;
  }

  // write out histograms
  if (!(f_out = fopen("pz.rms", "w"))) return;
  for (i=0; i<HIS_COUNT; i++) {
    char spname[256];
    if (i < 200) {
      sprintf(spname, "%d; PZ correction info: BL, 1/tau, frac2, lamda, chisq; run %d", i, runInfo->runNumber);
    } else {
      sprintf(spname, "%d; ch %d baseline RMS; t95; t100-t95", i, i%200);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  // adjust PZI values using fitted values stored in his[chan]
  // find new tau, frac2 values
  printf("Fitted PZ values:\n"
         "Chan   BL ;    tau     1/tau   ->    tau     1/tau  ;   frac2   ->   frac2\n");
  //        0   -61 ;   58.58   0.01707  ->   58.56   0.01708 ;   0.01200 ->  0.01167

  for (chan=clo; chan<=chi; chan++) {
    j = 10;
    for (i=11; i<1900; i++)
      if (his[chan][j] < his[chan][i]) j = i;  // mode of baseline distribution
    if (his[chan][j] < 200) continue;

    printf("%3d %5.0f ;  %6.2f %9.5f  -> ",
           chan, PZI.baseline[chan], PZI.tau[chan], 1.0/PZI.tau[chan]);
    /* first find mode of tau distribution */
    j = 2010;
    for (i=2011; i<3500; i++) {
      if (his[chan][j] < his[chan][i]) j = i;  // mode of tau distribution
    }
    if (his[chan][j] > 100) {
      PZI.tau[chan] = 30000.0 / (j-2000.0);
      fwhm = 5;
      // now integrate over +- 1.5 FWHM
      if ((pos = autopeak4(his[chan], j-20, j+20, 1.5f, &area, &fwhm)) && area > 100) {
	PZI.tau[chan] = 30000.0 / (pos-2000.0);
      }
      // if (chan <= 50) printf("\n chan %d area pos fwhm = %f %f %f  j = %d %d\n", chan, area, pos, fwhm, j, his[chan][j]);
      printf(" %6.2f %9.5f ;  %4d %.2f %.2f ; %8.5f -> ", PZI.tau[chan], 1.0/PZI.tau[chan], j, pos, fwhm, PZI.frac2[chan]);
    } else {
      printf("\n");
      continue;
    }

    fwhm = 8;
    if ((pos = autopeak4(his[chan], 4000, 4500, 1.0f, &area, &fwhm)) && area > 100)
      PZI.frac2[chan] = (pos - 4000.0) / 4000.0;
    // Use mode of distribution instead?
    j = 4000;
    for (i=4001; i<4500; i++)
      if (his[chan][j] < his[chan][i]) j = i;  // mode of frac2 distribution
    if (his[chan][j] > 100) {
      PZI.frac2[chan] = (j - 4000.0) / 4000.0;
      printf("%8.5f\n", PZI.frac2[chan]);
    } else {
      printf("\n");
    }
  }
  //------------------------------------------------

  PZ_info_write(runInfo, &PZI);
  return;
} /* signalselect() */

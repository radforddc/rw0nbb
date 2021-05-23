#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define VERBOSE 0
#define DELTATIME 1  // resolving time for timestamp selection
#define FILTER_10MHZ 0

void signalselect(FILE *f_in, MJDetInfo *detInfo, MJRunInfo *runInfo);

int main(int argc, char **argv) {

  FILE       *f_in, *f_lis=0;
  MJDetInfo  detInfo[NMJDETS];
  int        nDets;
  MJRunInfo  runInfo;
  int        argn=1;
  char       *c, fname[256], line[256];

  if (argc < 6) {
    fprintf(stderr, "\nusage: %s fname_in chnum_lo chnum_hi e_lo e_hi [-nakKdpzilt]\n"
            "     -n: do not subtract initial baseline\n"
            "     -z: subtract final signal, so that shifted value is zero after rise\n"
            "     -a: align signals to t80 = sample 1100\n"
            "     -k: e_lo, e_hi are in keV rather than ADC\n"
            "     -K: e_lo, e_hi are in keV, calculated using PZ correction\n"
            "     -p: do PZ correction\n"
            "     -d: take derivative of signal\n"
            "     -i: do not try to correct for INL\n"
            "     -l: select only events passing A/E and DCR but failing lq cut\n"
            "     -t: read a list of timestamps to get from ts.input\n\n", argv[0]);
    return -1;
  }
  /* open raw data file as input */
  if (strstr(argv[argn], ".lis")) {    // input argument is a list file
    if ((f_lis = fopen(argv[argn],"r")) == NULL) {
      fprintf(stderr, "\n Failed to open list input file %s\n", argv[argn]);
      return 0;
    }
    if (!fgets(fname, sizeof(fname), f_lis)) {
      fprintf(stderr, "\n Failed to read list input file %s\n", argv[argn]);
      return 0;
    }
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
  } else {                                      // use command argument as input file
    strncpy(fname, argv[argn], sizeof(fname));
  }
  while (argn < argc && argv[argn][0] == '-') argn ++;
  if ((f_in = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", fname);
    return 0;
  }
  printf("\n >>> Reading %s\n\n", fname);
  strncpy(runInfo.filename, fname, 256);
  runInfo.argc = argc;
  runInfo.argv = argv;

  /* read file header */
  nDets = decode_runfile_header(f_in, detInfo, &runInfo);
  if (nDets < 1) return 1;

  if (runInfo.dataIdGM == 0 && runInfo.dataIdGA == 0)
    printf("\n No data ID found for Gretina4M or 4A data!\n");
  if (runInfo.dataIdGM)
    printf("\n Data ID %d found for Gretina4M data\n", runInfo.dataIdGM);
  if (runInfo.dataIdGA)
    printf("\n Data ID %d found for Gretina4A data\n", runInfo.dataIdGA);
  printf(" Run number: %d in file %s\n", runInfo.runNumber, runInfo.filename);

  /* loop over all input files */
  while (1) {

    /* read through all the events in the file */
    runInfo.analysisPass = 0;
    signalselect(f_in, detInfo, &runInfo);
    fclose(f_in);
    if (!f_lis) break;

    /* get next input file name */
    if (!fgets(fname, sizeof(fname), f_lis))  // no more lines in the input list file
      strncpy(fname, "0", strlen(fname));     // special string
    for (c = fname + strlen(fname); *c < '0'; c--) *c = 0;  // removing trailing junk
    if (strlen(fname) < 2) break;  // no more files to process

    /* open next file and read header size and run number */
    if ((f_in = fopen(fname,"r")) == NULL) {
      fprintf(stderr, "\n Failed to open input file %s\n", fname);
      return 0;
    }
    strncpy(runInfo.filename, fname, sizeof(runInfo.filename));
    printf("\n >>> Reading %s\n\n", fname);
    fread(&runInfo.fileHeaderLen, 1, sizeof(int), f_in);
    /* loop through the lines of the XML data until we find the run number */
    while (fgets(line, sizeof(line), f_in) && strncmp(line, "</plist>", 8)) {
      if (strstr(line, "<key>RunNumber</key>")) {
        fgets(line, sizeof(line), f_in);
        if (!(c=strstr(line, "<integer>")) ||
            (1 != sscanf(c+9, "%d", &runInfo.runNumber))) {
          fprintf(stderr, "\n ERROR decoding run number:\n %s\n", line);
          return -1;
        }
        break;
      }
    }
    /* position at start of data */
    fseek(f_in, 4*runInfo.fileHeaderLen, SEEK_SET);
  }

  printf("\n All Done.\n\n");
  return 0;
}

/* ========================================================== */

void signalselect(FILE *f_in, MJDetInfo *Dets, MJRunInfo *runInfo) {

  int    i, j, k, crate=0, slot=0, evlen, ch, chan, sig_len, dd0, dd1, d;
  int    board_type, dataIdRun=0;
  unsigned int  head[2], evtdat[20000];
  static int    totevts=0, out_evts=0, clo, chi, elo, ehi;
  static int    doINL = 1, doPZ = 0, doLQ = 0, current_runNumber = 0;
  static long long int time, tsget[2000] = {0};

  static int module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  static PZinfo  PZI;
  static PSAinfo PSA;

  static FILE   *f_out=0, *f_ts = 0;
  static int    argn = 1, subbl = 1, kev = 0, align = 0, ts_sel = 0, erro, ntsget = 0, norm = 0, deriv=0;
  short  *signal, ucsig[8192], sigu[2048];  // ucsig = uncompressed signal (presumming corrected)
  float  fsignal[8192] = {0}, fs2[8192] = {0}, e1;
  float  drift, aovere, dcr, lq;
  double e_raw, e_adc, e_ctc;
  int    t80, t95, t100, bl, step, tsmatch = 0;
  char   line[256];

 
  /* initialize some values */
  if (!f_out) {
    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return;

    for (argn = 1; argn < runInfo->argc; argn++) {
      if (runInfo->argv[argn][0] == '-') {
        if (strstr(runInfo->argv[argn], "n")) subbl = 0;
        if (strstr(runInfo->argv[argn], "z")) norm  = 1;
        if (strstr(runInfo->argv[argn], "a")) align = 1;
        if (strstr(runInfo->argv[argn], "k")) kev   = 1;
        if (strstr(runInfo->argv[argn], "K")) kev   = 2;
        if (strstr(runInfo->argv[argn], "p")) doPZ  = 1;
        if (strstr(runInfo->argv[argn], "d")) deriv = 1;
        if (strstr(runInfo->argv[argn], "i")) doINL = 0;
        if (strstr(runInfo->argv[argn], "l")) doLQ  = 1;
        if (strstr(runInfo->argv[argn], "t") &&
            (f_ts = fopen("ts.input", "r"))) ts_sel = 1;
      }
    }
    argn = 2;
    clo = atoi(runInfo->argv[argn++]);
    chi = atoi(runInfo->argv[argn++]);
    elo = atoi(runInfo->argv[argn++]);
    ehi = atoi(runInfo->argv[argn++]);
    printf("\nChs %d to %d, e_trapmax %d to %d\n\n", clo, chi, elo, ehi);
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "n")) subbl = 0;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "z")) norm  = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "a")) align = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "k")) kev   = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "K")) kev   = 2;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "p")) doPZ  = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "d")) deriv = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "i")) doINL = 0;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "l")) doLQ  = 1;
    if (runInfo->argc > argn && strstr(runInfo->argv[argn], "t") &&
        (f_ts = fopen("ts.input", "r"))) ts_sel = 1;

    if (doLQ) doPZ = 1;

    /* read A/E, DCR, and lamda values from psa.input */
    /* also read individual trapezoid values from filters.input (if it exists) */
    if (!PSA_info_read(runInfo, &PSA)) {
      printf("\n ERROR: No initial PSA data read. Does psa.input exist?\n");
      if (doLQ) {
        printf("   LQ cut requires psa.input; aborting.\n");
        exit(-1);
      }
    }

    f_out = fopen("s.rms", "w");
    if (kev==2 || doPZ) doPZ = PZ_info_read(runInfo, &PZI);
    if (kev==2 && !doPZ) {
      printf("\n ERROR: No initial pole-zero data read, needed for keV. Does PZ.input exist?\n");
      return;
    }

    while (ts_sel &&
           ntsget < 2000 &&
           fgets(line, sizeof(line), f_ts)) {
      if (strlen(line) > 5 && line[0] != '#' &&
          sscanf(line, "%lld", &tsget[ntsget]) == 1) {
        if (tsget[ntsget] > 168 ||
            sscanf(line, "%d %f %lld", &j, &e1, &tsget[ntsget]) == 3) ntsget++;
      }
    }
    if (ts_sel) printf(" >>>> %d timestamps to look for matches\n", ntsget);

  }
  for (i=0; i<runInfo->idNum; i++) {
    if (strstr(runInfo->decoder[i], "ORRunDecoderForRun")) {
      dataIdRun = runInfo->dataId[i];
      printf("dataIdRun = %d %s %d\n", dataIdRun, runInfo->decoder[i], i);
    }
  }
  printf("dataIdRun = %d\n", dataIdRun);

  /* start loop over reading events from input file
     ============================================== */

  while (fread(head, sizeof(head), 1, f_in) == 1) {

    board_type = head[0] >> 18;
    evlen = (head[0] & 0x3ffff);

    if (board_type == 0) {  // a new runfile header! file must be corrupt?
      printf("\n >>>> ERROR: DataID = 0; found a second file header??"
             " Ending scan of this file!\n"
             " >>>> head = %8.8x %8.8x  evlen = %d\n", head[0], head[1], evlen);
      break;
    }

    if (board_type == dataIdRun) {
      fread(evtdat, 8, 1, f_in);
      if (head[1] & 0x21) {
        // printf("------- START Run %d at %d", evtdat[0], evtdat[1]);
        current_runNumber = evtdat[0];
      }
      continue;
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
    if (++totevts % 20000 == 0) {
      printf(" %8d evts in, %d out\n", totevts, out_evts); fflush(stdout);
    }

    /* --------------- Gretina4M digitizer data ---------------- */
    time = (evtdat[3] & 0xffff);
    time = time << 32 | evtdat[2];
    if (time < 0) continue;
    tsmatch = 0;
    for (i=0; i<ntsget; i++) {
      if (time >= tsget[i] - DELTATIME && time <= tsget[i] + DELTATIME) {
        ch = (evtdat[1] & 0xf);
        chan = chan_lu[module_lu[crate][slot]][ch];
        printf(" >>>> timestamp %14lld matches, chan %3d, run %6d\n",
               time, chan, current_runNumber);
        tsmatch = 1;
        break;
      }
    }
      
    ch = (evtdat[1] & 0xf);
    if ((j = module_lu[crate][slot]) >= 0 && ch < 10) {
      chan = chan_lu[j][ch];
      if (chan >= 0 && chan <= 157 &&
          ((chan < 100 && !Dets[chan].HGChEnabled) ||
           (chan > 99 && !Dets[chan-100].LGChEnabled))) {
        printf("Data from detector not enabled! Chan = %d  crate, slot, j, ch = %d %d %d %d\n", chan, crate, slot, j, ch);
        chan = -1;
      }
      if (chan < clo || chan > chi) continue;
      signal = (short *) evtdat + 28;
      if (evlen != 1026 && signal[0] == 2020 && // signal is compressed; decompress it
          2020 == decompress_signal((unsigned short *)signal, sigu, 2*(evlen - 2) - 28)) {
        signal = sigu;
        evlen = 1026;
      }

      sig_len = 2012;
      if (chan <= 99) {
        k = Dets[chan].HGPrerecnt + Dets[chan].HGPostrecnt;  // expected location of presumming start
      } else {
        k = Dets[chan-100].LGPrerecnt + Dets[chan-100].LGPostrecnt;
      }
      erro = 0;
      if (k > 10 && k < sig_len) { // correct for presumming
        // FIXME: replace hard-coded factor of 4 with data from header
        step = k;
        j = (signal[step-5] + signal[step-4] + signal[step-3] + signal[step-2])/4;
        if (j > 10 || j < -10) {
          /* --- Find transition to 4x presumming ------- */
          for (step=k; step<k+2; step++) {
            if (j > 0 && signal[step] > 2*j) break;
            if (j < 0 && signal[step] < 2*j) break;
          }
          if (step == k+2) {
            step = k+1;
            erro = 1;
            printf("Hmmm... step not found in %d - %d in chan %3d, counts %d\n",
                   step-1, step, chan, j);
          }
        }
        
        /* --- uncompress --- */
        for (i=0; i<step; i++) ucsig[i+10] = signal[i];
        for (i=step; i<sig_len; i++) {               // NOTE: This is wrong for signal < 0
          for (j=0; j<4; j++)
            ucsig[step + 10 + 4*(i-step) + j] = signal[i]/4; // dist sums over 4 bins each
          for (j=0; j<signal[i]%4; j++)
            ucsig[step + 10 + 4*(i-step) + j]++;     // and make sure sum is right
        }
        sig_len += 3*(sig_len-step);
        signal = &ucsig[10];
      }
      /* --- Done with 4x presumming ------- */

      if (kev==2) {
        if (PZ_correct(signal, fsignal, sig_len, chan, &PZI))
          printf(" >>> PZ_correct return error for chan %d\n", chan);
        e1 = float_trap_max(fsignal, &j, 401, 250)/401.0;
      } else {
        e1 = trap_max(signal, &j, 401, 200)/401;
      }
      if (kev) {
        if (chan <= 57) {
          e1 *= Dets[chan].HGcalib[0];
        } else if (chan > 99 && chan <= 157) {
          e1 *= Dets[chan-100].LGcalib[0];
        }
      }

      if (// !erro &&
          (e1 < elo || e1 > ehi) &&
          !tsmatch) continue;

      /* sticky-bit fix */
      d = 128;  // basically a sensitivity threshold; max change found will be d/2
      if (chan > 99 && chan < 100+runInfo->nGe) d = 64;
      for (i=20; i<sig_len-10; i++) {
        // calculate second derivatives
        dd0 = abs((int) signal[i+1] - 2*((int) signal[i]) + (int) signal[i-1]);
        dd1 = abs((int) signal[i+2] - 2*((int) signal[i+1]) + (int) signal[i]);
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
          //signal[i] += j*d/2;
          break;
        }
      }

      /* find mean baseline value */
      bl = 0;
      for (i=300; i<400; i++) bl += signal[i];
      bl /= 100;

      // lost events for doLQ=0 between here
      //out_evts++;
      /* find t100, t95, and t80 (t80 will be redone later of -l is specified) */
      t100 = 700;                 // FIXME? arbitrary 700?
      for (i = t100+1; i < sig_len - 500; i++)
        if (signal[t100] < signal[i]) t100 = i;
      if (t100 > sig_len - 700) continue;  // FIXME ??  - important for cleaning, gets rid of pileup
      /* get mean baseline value */
      if (doPZ &&
          ((bl - PZI.baseline[chan]) < -10 || (bl - PZI.baseline[chan]) > 50)) continue;   // a little data cleaning
      i = bl + (signal[t100] - bl)*19/20;
      for (t95 = t100-1; t95 > 500; t95--)
        if (signal[t95] <= i) break;
      i = bl + (signal[t100] - bl)*4/5;
      for (t80 = t100-1; t80 > 500; t80--)
        if (signal[t80] <= i) break;
      // and here
      //out_evts++;

      /* do (optional) INL and PZ corrections */
      if (doINL) {
        if (inl_correct(Dets, runInfo, signal, fsignal, sig_len, chan))
          printf(" >>> inl_correct return error for chan %d\n", chan);
        if (doPZ) {
          if (PZ_fcorrect(fsignal, sig_len, chan, &PZI))
            printf(" >>> PZ_correct return error for chan %d\n", chan);
        }
        for (i=0; i<sig_len; i++) signal[i] = lrintf(fsignal[i]);
      } else if (doPZ) {
        if (PZ_correct(signal, fsignal, sig_len, chan, &PZI))
          printf(" >>> PZ_correct return error for chan %d\n", chan);
        for (i=0; i<sig_len; i++) signal[i] = lrintf(fsignal[i]);
      }
      if (1 && (doINL || doPZ)) bl = PZI.baseline[chan];

      if (doLQ) {
        /* get raw (e_raw) and drift-time-corrected energy (e_adc and e_ctc)
           and effective drift time */
        int tmax, t0;
        e_ctc = get_CTC_energy(fsignal, sig_len, chan, Dets, &PSA,
                               &t0, &e_adc, &e_raw, &drift, 0.0);
        if (e_ctc < 0.001) {
          printf("E_ctc = %.1f < 1 eV!\n", e_ctc);
          if (VERBOSE) printf("chan, t0, t100: %d %d %d; e, drift: %.2f %.2f\n",
                              chan, t0, t100, e_raw, drift);
          continue;
        }
        if (runInfo->flashcam) {
          e_ctc /= 2.0;
          e_adc /= 2.0;
          e_raw /= 2.0;
          drift /= 2.0;
        }

        /* find A/E */
        if (e_ctc > 50 && t0 > 200 && t0 < 1500) {
          float s2 = float_trap_max_range(fsignal, &tmax, PSA.a_e_rise[chan], 0, t0-20, t0+300); // FIXME: hardwired range??
          aovere = s2 * PSA.a_e_factor[chan] / e_raw;
          /* do quadratic fit/interpolation over +- one sample, to improve max A/E determination
             does this help? seems to work well in some cases, in other cases not so much?
          */
          if (AoE_quad_int) {
            float s1, s3, a, b, aoe;
            int rise = PSA.a_e_rise[chan];
            s1 = s2 - (fsignal[tmax  ] - 2.0*fsignal[tmax   + rise] + fsignal[tmax   + 2*rise]);
            s3 = s2 + (fsignal[tmax+1] - 2.0*fsignal[tmax+1 + rise] + fsignal[tmax+1 + 2*rise]);
            b = s2 - (s1+s3)/2.0;
            a = s2 - s1 + b;
            aoe = s1 + a*a/(4.0*b);
            if (aoe > s2) aovere = aoe * PSA.a_e_factor[chan] / e_raw;
          }
          if (runInfo->flashcam) aovere /= 2.0;
        } else {
          aovere = 0;
          if (e_raw > 1000)
            printf("Error getting A/E: chan %d   e_raw, e_ctc = %.0f, %.0f  t0 = %d\n", chan, e_raw, e_ctc, t0);
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
        float dtc = drift*3.0*2614.0/e_adc;   // in (x)us, for DT- correction to A/E
        if (chan%100 == 50) dtc *= 3.0;       // special hack for detector 50; has especially strong variation
        float ec = e_ctc/1000.0 - 2.0;        // in MeV, for E-correction to A/E
        float aovere1 = aovere + PSA.ae_dt_slope[chan] * dtc;
        float aovere2 = aovere1 + PSA.ae_e_slope[chan] * ec;
        //if (AOE_CORRECT_NOISE)
        //  aovere2 -= (2.0 * PZI.bl_rms[chan] * sqrt(2.0 * (float) PSA.a_e_rise[chan]) *
        //         PSA.a_e_factor[chan] * gain/1593.0 * (1.0 - 1593.0/e_ctc));        // 1593 keV = DEP
        aovere2 -= AOE_CORRECT_EDEP * (PSA.ae_pos[chan] - PSA.ae_cut[chan]) * (1.0 - e_ctc*e_ctc/1593.0/1593.0);
        if (PSA.ae_pos[chan] > 100 && aovere2 < PSA.ae_cut[chan]) continue; // fails A/E cut and doLQ == 1

        /* get DCR from slope of PZ-corrected signal tail */
        if (sig_len < 2450) {
          if (t95 + 760 > sig_len) {
            if (VERBOSE)
              printf("Error getting DCR in skim.c: t95 = %d, sig_len = %d\n", t95, sig_len);
            continue;
          }
          dcr = float_trap_fixed(fsignal, t95+50, 100, 500) / 25.0;
        } else {
          if (t95 + 1180 > sig_len) {
            if (VERBOSE)
              printf("Error getting DCR in skim.c: t95 = %d, sig_len = %d\n", t95, sig_len);
            continue;
          }
          dcr = float_trap_fixed(fsignal, t95+50, 160, 800) / 40.0;
        }
        dcr -= PSA.dcr_dt_slope[chan] * drift;
        if (dcr > PSA.dcr_lim[chan]) continue;         // fails DCR cut and doLQ == 1

        /* find bl, t80, lq from PZ-corrected signal */
        /* get mean baseline value */
        float fbl = 0;
        for (i=300; i<400; i++) fbl += fsignal[i];
        fbl /= 100.0;
        float ff = fbl + e_raw * 0.8;
        for (t80 = t95; t80 > 500; t80--)
          if (fsignal[t80] <= ff) break;
        lq = (fsignal[t80+1] - ff) / (fsignal[t80+1] - fsignal[t80]); // floating remainder for t80
        lq *= (e_raw + fbl - (fsignal[t80+1] + ff)/2.0);  // charge (not yet) arriving during that remainder

        /* get late charge (lq) = delay in charge arriving after t80 */
        lq += float_trap_fixed(fsignal, t80+1, 100, 100);
        lq /= e_raw/100.0;

        lq -= PSA.lq_dt_slope[chan] * dtc;  // do DT correction to lq
        if (lq < PSA.lq_lim[chan]) continue; // passes LQ cut and doLQ == 1 selects only LQ-failing events
        if (0)
          printf("chan %3d | A/E, cut: %6.1f %6.1f | DCR, cut: %5.1f %5.1f | lq, cut: %5.1f %5.1f\n",
                 chan, aovere2, PSA.ae_cut[chan], dcr, PSA.dcr_lim[chan], lq, PSA.lq_lim[chan]);
      }  // end   if (doLQ) {

      out_evts++;

      // optionally subtract initial baseline
      if (subbl)
        for (i=0; i<sig_len; i++) signal[i] -= bl;

      // optionally remove 10MHz noise
      if (FILTER_10MHZ) {
        if (!doINL && !doPZ) {
          for (i=0; i<sig_len; i++) fsignal[i] = signal[i];
        } else if (subbl) {
          for (i=0; i<sig_len; i++) fsignal[i] -= bl;
        }
        // 3-2-3 trap filter
        fs2[0] = 0;
        for (i=0; i < 3; i++) fs2[0] += fsignal[i+5] - fsignal[i];
        for (i=1; i < sig_len-7; i++) {
          fs2[i] = fs2[i-1] + fsignal[i+7] - fsignal[i+4] - fsignal[i+2] + fsignal[i-1];
        }
        // for (i=6; i < sig_len; i++) signal[i] = lrintf(fsignal[i] - 0.19*fs2[i-6]);
        for (i=6; i < sig_len-6; i++) signal[i] = lrintf(fsignal[i] - 0.095*fs2[i-6] + 0.095*fs2[i-1]);
      }

      // optionally subtract final signal value, so signal ends at 0 ADC
      if (norm) {
        d = 100;
        for (i=1700; i<1900; i++) d += signal[i];
        d /= 200;
        for (i=0; i<sig_len; i++)
          if (signal[i] > d/2) signal[i] -= d;
      }

      // optionally align t80 to t=1100
      if (align) {
        if (t80 > 1100) {
          j = t80 - 1100;
          for (i=1; i+j < sig_len; i++)
            signal[i] = signal[i+j];
        } else if (t80 < 1100) {
          j = 1100 - t80;
          for (i=sig_len - 1; i-j > 0; i--)
            signal[i] = signal[i-j];
        }
      }

      // optionally take derivative of signal (4-0-4 samples trap filter)
      if (deriv) {
        k = 0;
        for (i=0; i < 4; i++) k += signal[i+5] - signal[i+1];
        signal[0] = k;
        for (i=1; i < sig_len-9; i++) {
          k += signal[i+8] - 2*signal[i+4] + signal[i];
          signal[i] = k;
        }
        for (; i < sig_len; i++) signal[i]=0;
      }

      // optionally add some extra information at the start
      if (1) {
        signal  -= 8; // 8 extra words (shorts)
        sig_len += 8;
        for (i=0; i<8; i++) signal[i] = 0;
        signal[2] = (short) e1;
        signal[4] = (short) ((evtdat[8]>>11)&0xf)*100;  // CHECKME
        signal[6] = (short) chan;
      }
      write_sig(signal, sig_len, chan, f_out);
      // printf("output signal:  chan %d time %lld energy %.1f\n", chan, time, e1);
    }
  }

  //fclose(f_out);
  printf(" %8d evts in, %d out\n", totevts, out_evts);
  return;
} /* signalselect() */

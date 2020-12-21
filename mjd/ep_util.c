/*
   ep_util.c:  utility functions for eventprocess.c

   David Radford   Nov 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define DEBUG      0
#define VERBOSE    0

/*  ------------------------------------------------------------ */

void fillEvent(MJRunInfo *runInfo, int nChData, BdEvent *ChData[],
               int module_lu[][21], int det_lu[][16], int chan_lu[][16]) {
  
  int ievt, tmax;

  for (ievt = 0; ievt < nChData; ievt++) {

    if (ChData[ievt]->orca_type == runInfo->dataIdGM ||
        ChData[ievt]->orca_type == runInfo->dataIdGA) {
      ChData[ievt]->ch = (ChData[ievt]->evbuf[3] & 0xf);
      if ((ChData[ievt]->mod = module_lu[ChData[ievt]->crate][ChData[ievt]->slot]) >= 0 &&
          ChData[ievt]->ch <= 9) {
        ChData[ievt]->chan = chan_lu[ChData[ievt]->mod][ChData[ievt]->ch];
        ChData[ievt]->det  = det_lu[ChData[ievt]->mod][ChData[ievt]->ch];
        if (DEBUG) printf("energy %d %d = %.1f\n", ievt, ChData[ievt]->chan, ChData[ievt]->e);
      } else {
        ChData[ievt]->det = ChData[ievt]->chan = -1;
      }
      if (ChData[ievt]->chan < 0) continue;
      /* ChData[ievt]->sig  = (short *) ChData[ievt]->evbuf + 32; // now done in eventbuild.c */
      ChData[ievt]->e    = (float) trap_max(ChData[ievt]->sig, &tmax, TRAP_RISE, TRAP_FLAT) / (double) TRAP_RISE;
      // some pulser-tag channels have too high an energy, due to presumming. here is where I correct that.
      if (ChData[ievt]->chan > 99 + runInfo->nGe) {
        //if (ChData[ievt]->e < -10) ChData[ievt]->e = -ChData[ievt]->e;
        while (ChData[ievt]->e > 4000)
          ChData[ievt]->e /= 4.0;
      }
    }
  }
}

/*  ------------------------------------------------------------ */

int checkForPulserEvent(MJDetInfo *Dets, MJRunInfo *runInfo,
                        int nChData, BdEvent *ChData[], PTag *pt) {
  
  int      i, j, k, pulser, idet, ievt, chan;
  int      chan_list[200], hitchan[200];  // list of detectors/channels, and those meeting energy gate
  int      e_list[200];  // hit energy
  int      nchan, ndet, nCC[10];
  float    energy;
  long long int  time=0, dt;
  static int n_error_msgs = 0, n_error2_msgs = 0;


  for (i=0; i<10; i++) nCC[i] = 0;
  pulser = nchan = ndet = 0;
  
  /* LOOP over channel-events in the built-event once to identify pulser events */

  for (ievt = 0; ievt < nChData; ievt++) {
    chan_list[ievt] = -1;
    if (ChData[ievt]->orca_type == runInfo->dataIdGM ||
        ChData[ievt]->orca_type == runInfo->dataIdGA) {
      if (ChData[ievt]->chan < 0 ||
          ChData[ievt]->chan >= 100 + runInfo->nPT + runInfo->nGe) continue;
      chan_list[ievt] = chan = ChData[ievt]->chan;
      time = ChData[ievt]->time;
      e_list[ievt] = energy = ChData[ievt]->e;
      if (DEBUG) printf("energy %d %d = %.1f\n", ievt, chan, energy);

      if (energy >= pt->elo[chan] && energy <= pt->ehi[chan] && pt->pdt[chan] > 0) {
        // signal meets pulser energy gate; add it to a list of hit detectors
        for (i=0; i<nchan; i++)
          if (hitchan[i]%100 == chan%100) break;  // aleady saw this detector
        if (i == nchan) {   // this detector _not_ yet seen
          ndet++;
          /* if at least two different detectors on the same controller card
             are hit with the correct energy, then this must really be a pulser event...
             but we don't fully trust special detector 16 */
          if ((idet = chan%100) < runInfo->nGe && idet != 16 &&
              ++nCC[Dets[idet].CCnum] > 1) pulser = 1;
        }
        hitchan[nchan++] = chan;
        /* if a pulser tag channel is hit with the correct energy,
           then this must really be a pulser event */
        if (chan >= 100 + runInfo->nGe) pulser = 1;
      }
      /* if a pulser tag channel is hit even with a negative energy,
         then this must really be a pulser event */
      if (chan >= 100 + runInfo->nGe && energy < -20) pulser = 1;
    }
  }

  /* if at least two different detectors on the same controller card,
     or a pulser tag channel, are hit with the correct energy,
     then this must really be a pulser event */
  if (pulser) {
    for (j=0; j<nchan; j++) {
      chan = hitchan[j];
      idet = chan%100;
      // pt->pt0[chan] = time;      // use this as the most recent pulser time for hit channels
      if (idet < runInfo->nGe) pt->cct0[Dets[idet].CCnum] = time;    // and hit CC's
    }
  }

  /* else, if a detector is hit with the correct energy,
     and is the correct time since a previous pulser event from its controller card,
     then this must really be a pulser event*/
  if (!pulser) {
    for (j=0; j<nchan; j++) {
      chan = hitchan[j];
      idet = chan%100;
      time  = ChData[i]->time;
      if (idet >= runInfo->nGe) { // pulser-tag channel
        pulser = 1;
      } else if (idet != 16) {    // don't flag as pulser event based solely on det 16
        dt = (time - pt->cct0[Dets[idet].CCnum])/100;     // time since last pulser, 1 us/bin
        if (dt < -6) {  // FIXME: Why is dt ~ -4 sometimes???
          if (n_error2_msgs++ < 5)
            printf("Ackk... in pulser tag check, timestamps out-of order?\n"
                   "    ... chan %3d dt is %12lld for CC number %d\n",
                   chan, dt, Dets[idet].CCnum);
          continue;
        }
        k = (dt+100) / pt->ccdt[Dets[idet].CCnum];       // number of cycles
        if (DEBUG) { printf(" dt = %lld   k = %d\n", dt, k); fflush(stdout); }
        if (dt < 8 ||
            (dt+2*k+8) % pt->ccdt[Dets[idet].CCnum] <= 4*k+16) { // allow dt to slip by 2 us per cycle,
          // plus 8 us for the resolving time
          pulser = 1;
        } else {
          if (1 && ++n_error_msgs < 100) {
            printf("Hmm.. pulser tag check: chan %3d   "
                   "dt = %8lld   pdt = %8lld  remainder(%d) %8lld\n",
                   chan, dt, pt->ccdt[Dets[idet].CCnum], k,
                   dt % pt->ccdt[Dets[idet].CCnum]);
          } else if (n_error_msgs == 100) {
            printf("******* Futher pulser tag error messages suppressed! *******\n");
          }
        }
      }
    }
  }

  /* now count the pulser events for each channel, for dead-time evaluation
     ptag.nevts[chan][0]: energy-ungated HG pulser cts
     ptag.nevts[chan][0]: energy-ungated LG pulser cts
     ptag.nevts[chan][0]: energy-ungated HG&&LG pulser cts
     ptag.nevts[chan][0]: energy-gated HG pulser cts
     ptag.nevts[chan][0]: energy-gated LG pulser cts
     ptag.nevts[chan][0]: energy-gated HG&&LG pulser cts
   */
  if (pulser) {
    for (j=0; j<nchan; j++) {
      chan = hitchan[j];
      if (chan < 0) continue;
      idet = chan%100;
      if (idet >= runInfo->nGe ||
          nCC[Dets[idet].CCnum] > 0) {
        pt->nevts[idet][3+chan/100]++;     // energy-gated pulser count
        for (i = 0; i < j; i++)
          if (hitchan[i]%100 == idet) pt->nevts[idet][5]++;  // both LG and HG channels see this pulser signal
      }
    }
    for (j = 0; j < nChData; j++) {
      if (chan_list[j] >= 0 && chan_list[j] < 200 &&
          e_list[j] > pt->elo[chan_list[j]]/2) {     // coarse check that energy is at least elo/2
        idet = chan_list[j]%100;
        if (idet >= runInfo->nGe ||
            nCC[Dets[idet].CCnum] > 0) {
          pt->nevts[idet][chan_list[j]/100]++;         // energy-ungated pulser count
          for (i = 0; i < j; i++)
            if (chan_list[i] >= 0 && chan_list[i] < 200 &&
                e_list[i] > pt->elo[chan_list[i]]/2 &&           // coarse check that energy is at least elo/2
                chan_list[i]%100 == idet) pt->nevts[idet][2]++;  // both LG and HG channels see this pulser signal
        }
      }
    }
  }

  /* for (j = 0; j < nChData; j++)
     if (chan_list[j] == chan) this_ch_is_pulser[ievt] = 1; */

  return pulser;
}

/*  ------------------------------------------------------------ */

int ep_init(MJDetInfo *Dets, MJRunInfo *runInfo, int module_lu[NCRATES+1][21],
            int det_lu[NBDS][16], int chan_lu[NBDS][16]) {

  int    i, j, k, n;
  int    nm = 0;      // number of VME modules defined
  FILE   *f_in;
  double g1, g2;
  char   line[256];

  /* --------------- initialize --------------- */
  for (i=0; i<=NCRATES; i++)
    for (j=0; j<21; j++) module_lu[i][j] = -1;
  for (i=0; i<NBDS; i++)
    for (j=0; j<16; j++) det_lu[i][j] = chan_lu[i][j] = -1;

  // some default values
  for (i=0; i<runInfo->nGe; i++) {
    Dets[i].HGcalib[0] = 0.5;
    Dets[i].LGcalib[0] = 1.5;
  }
  /* try to read energy calibration coefficients */
  if ((f_in = fopen(ECAL_FILENAME,"r"))) {
    printf("\n Reading energy calibrations from %s\n", ECAL_FILENAME);
    while (fgets(line, sizeof(line), f_in)) {
      if (*line != '#' &&
          sscanf(line, "%d %lf %lf", &i, &g1, &g2) == 3 &&
          i >=0 && i < runInfo->nGe) {
        Dets[i].HGcalib[0] = g1;
        Dets[i].LGcalib[0] = g2;
      }
    }
  } else {
    printf("\n file %s does not exist. runInfo->nGe = %d\n", ECAL_FILENAME, runInfo->nGe);
  }

  /* set up lookup table to map VME crate&slot onto VME modules */
  nm = 0;
  for (i=0; i<runInfo->nGD; i++) {
    if (runInfo->GDcrate[i] < 0 || runInfo->GDcrate[i] > NCRATES) {
      printf("ERROR: Illegal crate number %d > %d!\n\n",
             runInfo->GDcrate[i], NCRATES);
      return -1;
    }
    module_lu[runInfo->GDcrate[i]][runInfo->GDslot[i]] = nm++;
  }
  for (i=0; i<runInfo->nV; i++) {
    if (runInfo->Vcrate[i] < 0 || runInfo->Vcrate[i] > NCRATES) {
      printf("ERROR: Illegal crate number %d > %d!\n\n",
             runInfo->Vcrate[i], NCRATES);
      return -1;
    }
    module_lu[runInfo->Vcrate[i]][runInfo->Vslot[i]] = nm++;
  }
  if (nm > NBDS) {
    printf("ERROR: Too many VME boards (%d)!\n\n", nm);
    return -1;
  }

  /* set up lookup table to map modules and channels onto detector IDs */
  for (i = 0; i < runInfo->nGe; i++) {
    if ((j = module_lu[Dets[i].crate][Dets[i].slot]) < 0) {
      printf("ERROR; Ge det %d has unknown crate and slot (%d %d)!\n\n",
             i, Dets[i].crate, Dets[i].slot);
    } else {
      det_lu[j][Dets[i].chanHi] = i;
      det_lu[j][Dets[i].chanLo] = i;
      chan_lu[j][Dets[i].chanHi] = i;
      chan_lu[j][Dets[i].chanLo] = i+100;
      if (0) printf("Ge %2d crate, slot, ch = %d %2d %d -> chan_lu %2d %3d\n",
                    i, Dets[i].crate, Dets[i].slot, Dets[i].chanHi, i, i+100);
    }
  }
  for (i=0; i<runInfo->nPT; i++) {
    if ((j = module_lu[runInfo->PTcrate[i]][runInfo->PTslot[i]]) < 0) {
      printf("ERROR; Pulser tag %d has unknown crate and slot (%d %d)!\n\n",
             i, runInfo->PTcrate[i], runInfo->PTslot[i]);
    } else {
      det_lu[j][runInfo->PTchan[i]] = i+100+runInfo->nGe;
      chan_lu[j][runInfo->PTchan[i]] = i+100+runInfo->nGe;
      if (0) printf("PT %2d crate, slot, ch = %d %2d %d -> chan_lu %3d\n",
                    i, runInfo->PTcrate[i], runInfo->PTslot[i],
                    runInfo->PTchan[i], i+100+runInfo->nGe);
    }
  }
  n = 200;  // = starting detector ID for veto segments
  for (i=0; i<runInfo->nV; i++) {
    if ((j = module_lu[runInfo->Vcrate[i]][runInfo->Vslot[i]]) < 0) {
      printf("ERROR; Veto module %d has unknown crate and slot (%d %d)!\n\n",
             i, runInfo->Vcrate[i], runInfo->Vslot[i]);
    } else {
      for (k=0; k<16; k++) {
        det_lu[j][k] = n;
        chan_lu[j][k] = n;
        if (0) printf("Veto %2d crate, slot, ch = %d %2d %d -> chan_lu %3d\n",
                    i, runInfo->Vcrate[i], runInfo->Vslot[i], k, n);
      }
    }
  }
  if (nm > NBDS) {
    printf("ERROR: Too many VME boards (%d)!\n\n", nm);
    return -1;
  }

  /* double-check for errors in chan_lu tables */
  for (i=0; i<NBDS; i++) {
    for (j=0; j<16; j++) {
      if (chan_lu[i][j] < -1 || chan_lu[i][j] > 200)
        printf(" >>>>>>>>  ERROR: chan_lu[%d][%d] = %d!\n", i, j, chan_lu[i][j]);
    }
  }

  return  nm;
} /* ep_init() */

/* ---------------------------------------- */

int pulser_tag_info_write(MJRunInfo *runInfo, PTag *pt, int tell_results) {

  int   i;
  char  fname[256], *c;
  FILE  *f_out = 0;

  c = runInfo->filename;
  while (strchr(c, '/')) c = strchr(c, '/') + 1;  // remove initial path
  if (!strncmp(c, "DS", 2)) {                     // input data file starts with "DS"
    sprintf(fname, "%s.pdt", c);
  } else {
    sprintf(fname, "run%d.pdt", runInfo->runNumber);
  }
  if (0 && (f_out = fopen(fname, "r"))) {  // file already exists; optionally modify fname
    fclose(f_out);
    strncat(fname, "_new", sizeof(fname)-5);
  }
  if (!(f_out = fopen(fname, "w"))) {
    printf("ERROR in pulser_tag_info_write(): Cannot open output file %s\n", fname);
    return 1;
  }

  fprintf(f_out,
          "# Pulser tag data, run %d\n"
          "# DetID   pdt (us)         t0 (10ns)  eloHG  ehiHG     eloLG  ehiLG\n",
          runInfo->runNumber);
  if (tell_results)
    printf("  DetID   pdt (us)         t0 (10ns)  eloHG  ehiHG     eloLG  ehiLG\n");
    
  for (i=0; i<69; i++) {
    if (pt->pdt[i] == 0) {
      pt->pt0[i] = pt->pt0[i+100];
      pt->pdt[i] = pt->pdt[i+100];
    }
    if (pt->elo[i] <= pt->ehi[i] || pt->elo[i+100] <= pt->ehi[i+100]) {
      fprintf(f_out,
              "%5d %11lld %15lld %9d %6d %9d %6d\n",
              i, pt->pdt[i], pt->pt0[i], pt->elo[i], pt->ehi[i],
              pt->elo[i+100], pt->ehi[i+100]);
      if (tell_results)
        printf("%5d %11lld %15lld %9d %6d %9d %6d\n",
               i, pt->pdt[i], pt->pt0[i], pt->elo[i], pt->ehi[i],
               pt->elo[i+100], pt->ehi[i+100]);
    }
  }

  fclose(f_out);
  return 0;
}

/* ---------------------------------------- */

int pulser_tag_info_read(MJRunInfo *runInfo, PTag *pt) {

  int   i, n=0;
  char  fname[256], fname2[256], line[256], *c;
  FILE  *f_in = 0;

  /*
    read data pulser tag info from file
    look first for <data_filename>.pdt, then run<number>.pdt, then default.pdt
  */
  c = runInfo->filename;
  while (strchr(c, '/')) c = strchr(c, '/') + 1;  // remove initial path
  sprintf(fname, "%s.pdt", c);
  if ((f_in = fopen(fname, "r"))) {
    printf("\n Initializing pulser tag from file %s\n\n", fname);
  } else {
    sprintf(fname, "run%d.pdt", runInfo->runNumber);
    sprintf(fname2, "run%d.pdt", runInfo->firstRunNumber);
    if ((f_in = fopen(fname, "r"))) {
      printf("\n Initializing pulser tag from file %s\n", fname);
    } else if ((f_in = fopen(fname2, "r"))) {
      printf("\n Initializing pulser tag from file %s\n", fname2);
    } else if ((f_in = fopen("default.pdt", "r"))) {
      printf("\n Initializing pulser tag from file default.pdt\n");
    } else {
      printf("\n No pulser tag info file; tried %s.pdt, %s, %s, and default.pdt\n",
             c, fname, fname2);
      return 1;
    }
  }

  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#') continue;
    if (sscanf(line, "%5d", &i) != 1 ||
        i < 0 || i > 99 || 
        sscanf(line+6, "%11lld %15lld %9d %6d %9d %6d",
               &pt->pdt[i], &pt->pt0[i], &pt->elo[i],
               &pt->ehi[i], &pt->elo[i+100], &pt->ehi[i+100]) != 6) {
      printf(" Error; reading of pulser tag file failed!\n  %s\n", line);
      fclose(f_in);
      return 1;
    }
    pt->pdt[i+100] = pt->pdt[i];
    pt->pt0[i+100] = pt->pt0[i];
    n++;
  }

  if (VERBOSE) {
    printf("  ... %d detectors and tags initialized.\n", n);
    printf("  DetID   pdt (us)             t0     eloHG  ehiHG     eloLG  ehiLG\n");
    for (i=0; i<69; i++) {
      if (pt->elo[i] <= pt->ehi[i] || pt->elo[i+100] <= pt->ehi[i+100])
        printf("%5d %11lld %15lld %9d %6d %9d %6d\n",
               i, pt->pdt[i], pt->pt0[i], pt->elo[i], pt->ehi[i],
               pt->elo[i+100], pt->ehi[i+100]);
    }
  }

  fclose(f_in);
  return 0;
}

/* ---------------------------------------- */

int pulser_tag_init(MJDetInfo *Dets, MJRunInfo *runInfo, PTag *pt) {

  int   i;

  /* --------------- initialize pulser tag --------------- */
  for (i=0; i<200; i++) {
    pt->pdt[i] = pt->pt0[i] = pt->ehi[i] = 0;
    pt->elo[i] = 1;
  }
  memset(pt->nevts[0], 0, sizeof(pt->nevts));

  if (pulser_tag_info_read(runInfo, pt)) return 0;

  if (SPECIAL_DET_16) Dets[16].CCnum = 1;
  /* not the true CC number, but the one that gives the largest pulser signal */

  /* set up common pulser tag times, using controller cards */
  for (i=0; i<8; i++)
    pt->ccdt[i] = pt->cct0[i] = -1;
  for (i=0; i < runInfo->nGe; i++) {
    if (pt->elo[i] <= pt->ehi[i]) {
      pt->ccdt[Dets[i].CCnum] = pt->pdt[i];
      pt->cct0[Dets[i].CCnum] = pt->pt0[i];
    }
  }
  return 1;
} /* pulser_tag_init() */

/*  ------------------------------------------------------------ */

/* malloc INL storage and read complete (unformatted) INL data file */
int read_inl_unformatted(float **inl, int chan, FILE *f_in) {
  int k = -1;

  if (!(inl[0] = calloc(16384, sizeof(float))) ||
      !(inl[1] = calloc(16384, sizeof(float)))) {
    printf("ERROR in read_inl_unformatted; cannot malloc inl!\n");
    exit(-1);
  }
  while (k != chan) {
    if (1 != fread(&k, sizeof(k), 1, f_in) ||
        1 != fread(inl[0], sizeof(float)*16384, 1, f_in) ||
        1 != fread(inl[1], sizeof(float)*16384, 1, f_in)) {
      free(inl[0]);
      free(inl[1]);
      inl[0] = inl[1] = 0;
      return 0;
    }
    //printf(" k= %d, chan = %d...\n", k, chan);
  }
  return 1;
}

/* malloc INL storage and read individual (formatted) INL data files */
int read_inl_formatted(float **inl, int chan, int crate, int slot, int ch) {
  char  fname[256];
  int k;
  float f;
  FILE  *f_in;

  if (!(inl[0] = calloc(16384, sizeof(float))) ||
      !(inl[1] = calloc(16384, sizeof(float)))) {
    printf("ERROR in read_inl_formatted; cannot malloc inl!\n");
    exit(-1);
  }

  sprintf(fname, "%s/Crate%d_GRET%d_Ch%d_part1a.dat",
          PATH_TO_NONLIN_DATA, crate, slot, ch);
  if (VERBOSE) printf("Nonlin File %s\n", fname);
  if (!(f_in = fopen(fname, "r"))) {
    printf("No INL data file! %s\n\n", fname);
    exit(-1);
  }
  while(fscanf(f_in, "%d %f\n", &k, &f) == 2) inl[0][k+8192] = f;
  fclose(f_in);

  sprintf(fname, "%s/Crate%d_GRET%d_Ch%d_part2a.dat",
          PATH_TO_NONLIN_DATA, crate, slot, ch);
  if (!(f_in = fopen(fname, "r"))) {
    printf("No INL data file! %s\n\n", fname);
    exit(-1);
  }
  while(fscanf(f_in, "%d %f\n", &k, &f) == 2) inl[1][k+8192] = f;
  fclose(f_in);

  return 0;
}

int inl_correct(MJDetInfo *Dets, MJRunInfo *runInfo,
                short *signal_in, float *signal_out, int len, int chan) {

  static int   first = 1;
  static float *inl[200][2] = {{0}};
  // static float inl_tau = 0;  // set inl_tau to zero to use time-independent correction
  float  inl_tau = 140;
  float  current_inl, inl_factor = 1.0;
  int    i, j;

  
  if (first) {
    /* ---------- initialize; malloc INL storage and read values ---------- */
    /* read in all required INL data */
    char  fname[256];
    FILE  *f_in, *f_out;
    int   n = 0;

    /* first try to read unformatted file with data for all channels */
    sprintf(fname, "%s/required_inl_data.unformatted", PATH_TO_NONLIN_DATA);
    if ((f_in = fopen(fname, "r"))) {
      printf(" Attempting to read INL data from file %s\n", fname);
      for (i=0; i < runInfo->nGe; i++) {
        inl[i][0] = inl[i+100][0] = 0;
        if (Dets[i].HGChEnabled) {  // high gain channel
          n++;
          if (!read_inl_unformatted(inl[i], i, f_in)) break;
        }
        if (Dets[i].LGChEnabled) {  // low gain channel
          n++;
          if (!read_inl_unformatted(inl[i+100], i+100, f_in)) break;
        }
      }
      fclose(f_in);
      if (i < runInfo->nGe) {
        n = 0;
        printf(" Error, that doesn't meet our needs... data for detector %d is missing?\n"
               "  ... Reading original files\n", i);
        for (j=0; j < i; j++) {
          if (inl[j][0]) free(inl[j][0]);
          if (inl[j][1]) free(inl[j][1]);          // FIXME - this is not tested
          if (inl[j+100][0]) free(inl[j+100][0]);
          if (inl[j+100][1]) free(inl[j+100][1]);  // FIXME - this is not tested
          inl[j][0] = inl[j][1] = 0;
          inl[j+100][0] = inl[j+100][1] = 0;
        }
      }
    }

    if (n == 0) {
      /* that didn't work... try to read individual formatted files */
      for (i=0; i < runInfo->nGe; i++) {
        inl[i][0] = inl[i+100][0] = 0;
        if (Dets[i].HGChEnabled) {  // high gain channel
          read_inl_formatted(inl[i], i, Dets[i].crate, Dets[i].slot, Dets[i].chanHi);
          n++;
        }
        if (Dets[i].LGChEnabled) {  // low gain channel
          read_inl_formatted(inl[i+100], i+100, Dets[i].crate, Dets[i].slot, Dets[i].chanLo);
          n++;
        }
      }
      /* that worked... now try to save a new summary file */
      sprintf(fname, "%s/required_inl_data.unformatted", PATH_TO_NONLIN_DATA);
      if ((f_out = fopen(fname, "w"))) {
        for (i=0; i < runInfo->nGe; i++) {
          j = i;  // high gain channel
          if (inl[j][0]) {
            fwrite(&j, sizeof(j), 1, f_out);
            fwrite(inl[j][0], sizeof(float)*16384, 1, f_out);
            fwrite(inl[j][1], sizeof(float)*16384, 1, f_out);
          }
          j = i+100;  // low gain channel
          if (inl[j][0]) {
            fwrite(&j, sizeof(j), 1, f_out);
            fwrite(inl[j][0], sizeof(float)*16384, 1, f_out);
            fwrite(inl[j][1], sizeof(float)*16384, 1, f_out);
          }
        }
        fclose(f_out);
        printf(" INL data for %3d channels saved in %s\n", n, fname);
      }
    }
    printf(" ... INL data read for %3d channels\n\n", n);
    first = 0;
  }

  /* check we have everything we need */
  if (len < 200 || signal_in == 0 || signal_out == 0 ||
      chan < 0 || chan > 199) return 1;
  if (chan >= 100 + runInfo->nGe) {
    for (i = 0; i < len; i++) signal_out[i] = (float) signal_in[i];
    return 0;
  }
  if (inl[chan][0] == 0 || inl[chan][1] == 0) {
    printf("\n ERROR in inl_correct(): inl data undefined for channel %d!\n", chan);
    for (i = 0; i < len; i++) signal_out[i] = (float) signal_in[i];
    return 1;
  }

  /* do INL correction
     inl[chan][0][adc] = course-step INL values, that show ~1.9 us delay
     inl[chan][1][adc] = remaining fine-step INL values, no apparent delay
   */
  /*
  j = 50;
  for (i = 10; i < 110; i++) j += 8192 + signal_in[i];
  current_inl = inl[chan][0][j/100]; // INL for mean ADC value over start of signal
  */
  if (inl_tau) {
    // use time constant in INL correction
    // first calculate mean INL  over start of signal
    current_inl = 0;
    for (i = 10; i < 110; i++) current_inl += inl[chan][0][signal_in[i] + 8192]/100.0;

    for (i = 0; i < 10; i++) // the first 10 samples can be presummed; treat them as special
      signal_out[i] = (float) signal_in[i] - current_inl;
    for (; i < len; i++) {
      current_inl += (inl[chan][0][signal_in[i] + 8192] - current_inl) / inl_tau;
      // inl_tau divisor => time constant ~ 1.4 - 1.9 us
      signal_out[i] = (float) signal_in[i] -
        inl_factor * (inl[chan][1][signal_in[i] + 8192] + current_inl);
    }
  } else {  // no time constant
    for (i = 0; i < 10; i++) // the first 10 samples can be presummed; treat them as special
      signal_out[i] = (float) signal_in[i];
    for (; i < len; i++) {
      signal_out[i] = (float) signal_in[i] -
        inl_factor * (inl[chan][0][signal_in[i] + 8192] + inl[chan][1][signal_in[i] + 8192]);
    }
  }

  return  0;
} /* inl_correct() */

/* ---------------------------------------- */

float float_trap_fixed(float *signal, int t0, int rise, int flat) {
  // FIXME: check against signal length
  int   i;
  float e = 0;

  for (i=0; i<rise; i++)
    e += signal[i+t0+rise+flat] - signal[i+t0];

  return e;
}

float float_trap_max(float *signal, int *tmax, int rise, int flat) {

  int   i, t=50;
  float e = 0, emax;

  for (i=0; i<rise; i++)
    e += signal[i+t+rise+flat] - signal[i+t];
  emax = e;
  *tmax = t-1;
  for (; t < 1950-2*rise-flat; t++) {  // FIXME: hard-coded signal length
    e += signal[t+2*rise+flat] - signal[t+rise+flat] - signal[t+rise] + signal[t];
    if (emax < e) {
      emax = e;
      *tmax = t;
    }
  }

  return emax;
}

float float_trap_max_range(float *signal, int *tmax, int rise, int flat, int tlo, int thi) {
   // FIXME: check against signal length
  int   i, t;
  float e = 0, emax;

  t = tlo;
  // if (thi + 2*rise+flat > 2000) thi = 2000 - 2*rise - flat;

  for (i=0; i<rise; i++)
    e += signal[i+t+rise+flat] - signal[i+t];
  emax = e;
  *tmax = t-1;
  for (; t < thi; t++) {
    e += signal[t+2*rise+flat] - signal[t+rise+flat] - signal[t+rise] + signal[t];
    if (emax < e) {
      emax = e;
      *tmax = t;
    }
  }

  return emax;
}

int trap_fixed(short *signal, int t0, int rise, int flat) {
  // FIXME: check against signal length
  int i, e = 0;

  for (i=0; i<rise; i++)
    e += signal[i+t0+rise+flat] - signal[i+t0];

  return e;
}

int trap_max(short *signal, int *tmax, int rise, int flat) {

  int i, e = 0, emax, t=50;

  for (i=0; i<rise; i++)
    e += signal[i+t+rise+flat] - signal[i+t];
  emax = e;
  *tmax = t;
  for (; t < 1950-2*rise-flat; t++) {  // FIXME: hard-coded signal length
    e += signal[t+2*rise+flat] - signal[t+rise+flat] - signal[t+rise] + signal[t];
    if (emax < e) {
      emax = e;
      *tmax = t;
    }
  }

  return emax;
}

int trap_max_range(short *signal, int *tmax, int rise, int flat, int tlo, int thi) {
  // FIXME: check against signal length
  int i, e = 0, emax, t;

  t = tlo;
  // if (thi + 2*rise+flat > 2000) thi = 2000 - 2*rise - flat;

  for (i=0; i<rise; i++)
    e += signal[i+t+rise+flat] - signal[i+t];
  emax = e;
  *tmax = t;
  for (; t < thi; t++) {
    e += signal[t+2*rise+flat] - signal[t+rise+flat] - signal[t+rise] + signal[t];
    if (emax < e) {
      emax = e;
      *tmax = t;
    }
  }

  return emax;
}

/* CFD: find timestamp when signal exceeds some fraction of its amplitude */

float sig_frac_time(short *signal, float fraction, int width, int tstart){

  int   e = 0, emax, t, tmax = tstart;
  float e0 = 0.0, e1, thresh;

  /* find average signal baseline */
  for (t = 50; t < 50 + 4*width; t++) e0 += (float) signal[t];
  e0 /= 4.0;
  /* find signal maximum */
  for (t = tstart; t < tstart+width; t++) e += signal[t];
  emax = e;
  for (t = tstart; t < 1950-width; t++) {
    e += signal[t+width] - signal[t];
    if (emax < e) {
      emax = e;
      tmax = t;
    }
  }

  /* find signal threshold */
  e = emax;
  thresh = e0 + (((float) emax) - e0) * fraction;
  for (t = tmax; t > 1; t--) {
    e += signal[t] - signal[t+width];
    if (e < thresh) break;
  }

  if (t < 2) return 0.0;
  e0 = e;
  e1 = e + signal[t+width] - signal[t];
  return t + (thresh-e0)/(e1-e0);
}

/* ---------------------------------------- */

int get_sig_t0(float *fsignal, int len, int chan, PSAinfo *PSA){

  /*
    input: fsignal = PZ-corrected signal waveform, of length len, from channel chan
           e = amplitude of signal in ADC units
    returns: t0 (starting sample of signal in 10 ns units)
             corrected for rise and flat times of trapezoid
  */

  int    j, t0, t1, tmax;
  float  ttrap[len], ttmax = 0;          // FIXME: variable len here could slow down code? use max_len?
  double s2, s3;
  int    dt = PSA->t0_rise[chan] + T0_TRAP_FLAT + T0_TRAP_FALL;

  t1 = 20 + PSA->e_ctc_rise[chan] - dt;      // start of range for calculation of asymmetric trap
  if (t1 < 20) t1= 20;
  len -=  PSA->e_ctc_rise[chan] + dt + 10;   // end of range for calculation of asymmetric trap
  //if (len > max_len) len = maxlen;         // FIXME

  // calculate time-pick-off asymmetric trapezoid
  s2 = s3 = 0;
  t0 = tmax = t1;
  ttrap[t0-1] = 0;
  for (j=0; j<PSA->t0_rise[chan]; j++) s2 += fsignal[t0 + j + T0_TRAP_FLAT + T0_TRAP_FALL];
  for (j=0; j<T0_TRAP_FALL; j++) s3 += fsignal[t0 + j];
  for (; t0 < len; t0++) {
    ttrap[t0] = s2 / (double) PSA->t0_rise[chan] - s3 / (double) T0_TRAP_FALL;
    if (ttmax < ttrap[t0]) {
      ttmax = ttrap[t0];
      tmax = t0;
    }
    if (ttrap[t0] > 100) break;
    s2 += fsignal[t0 + PSA->t0_rise[chan] + T0_TRAP_FLAT + T0_TRAP_FALL] -
          fsignal[t0 + T0_TRAP_FLAT + T0_TRAP_FALL];
    s3 += fsignal[t0 + T0_TRAP_FALL] - fsignal[t0];
  }

  // find threshold time
  if (DEBUG && ttmax < 100)
    printf("< t1, t0, ttmax, tmax: %d %d %.1f %d  ttrap[t1] : %.1f ", t1, t0, ttmax, tmax, ttrap[t1]);
  if (ttmax > PSA->t0_thresh[chan]) {
    for (t0 = tmax; t0 > t1 && ttrap[t0-1] > PSA->t0_thresh[chan]; t0--) ;
  } else {
    t0 = tmax;
  }
  if (DEBUG && ttmax < 100) printf(" << final t0: %d\n", t0); fflush(stdout);
  if (t0 <= t1) return -1;
  return t0 + dt;
} /* get_sig_t0() */

/* ---------------------------------------- */

double get_CTC_energy(float *fsignal, int len, int chan, MJDetInfo *Dets, PSAinfo *PSA,
                      int *t0, double *e_adc, double *e_raw, float *drift, float ctc_factor) {

  /*
    extract charge-trapping-corrected, PZ-corrected, and INL-corrected energy from signal

    input:
      fsignal: PZ-corrected and INL-corrected signal of length len, from channel chan
      Dets: MJ detector info  data structure
      PSA:  contains filter params to use for trapezoids
      CTC_factor: the value used in the correction, usually CTC.e_dt_slope[chan]
    outputs:
      returned value: energy in keV, or -1.0f in case of error
      t0: start time of drift/signal
      e_adc: energy in ADC units
      e_raw: uncorrected energy in 0.001 ADC units
      drift: charge trapping value (drift time * charge)
             to be used for optimizing correction, in ADC units
             CTC correction = drift*ctc_factor[chan]
  */
  double s, e_ctc;

  /* find t0 */
  *t0 = get_sig_t0(fsignal, len, chan, PSA);
  /* check t0 lower limit to make sure we have enough signal baseline to calculate the fixed trap */
  if (*t0 < PSA->e_ctc_rise[chan] - 31) return -1.0f;
  /* check t0 upper limit to make sure we have enough remaining signal to calculate the fixed trap */
  if (*t0 + PSA->e_ctc_rise[chan] + PSA->e_ctc_flat[chan - 10] >= len) return -1.0f;

  /* get new energy estimation based on *fixed* trap */
  *e_raw = float_trap_fixed(fsignal, *t0 - PSA->e_ctc_rise[chan] - 10,
                            PSA->e_ctc_rise[chan], PSA->e_ctc_flat[chan]) / (double) PSA->e_ctc_rise[chan];

  /* find time-related value for trapping correction */
  s = float_trap_fixed(fsignal, *t0-210, 200, 10)/200.0;  // FIXME?? increase 200 to TRAP_RISE or TRAP_FLAT?
  *drift = (0.7 * *e_raw - s) / 1000.0;   // drift time * charge for use in charge-trapping correction,
                                          //with some offset to get average near zero

  /* do optimum charge-trapping correction */
  *e_adc = *e_raw + *drift * ctc_factor;
  if (chan < 100) {
    e_ctc = *e_adc * Dets[chan].HGcalib[0];
  } else {
    e_ctc = *e_adc * Dets[chan-100].LGcalib[0];
  }

  if (DEBUG)
    printf("e_raw, e_adc, t0, s1, s2, drift: %6.1f %6.1f %6d %6.0f %6.2f  E_ctc = %.0f\n",
           *e_raw, *e_adc, *t0, s, *drift, e_ctc);

  return e_ctc;
}

/* ---------------------------------------- */

int peak_find(int *his, int lo, int hi) {    // find highest-energy peak between chs lo and hi
  // extremely simple algorithm, requires zero background counts; could be much improved
  int e, emax = hi;

  for (e=hi-1; e>lo; e--) {
    if (his[e] > his[emax]) emax = e;
    if (his[e] <= 0 && his[emax] > 5) return emax;
  }
  return 0;
}

/* ======================================================================= */
int find_cent(int *his, int loch, int hich,
              float *area_ret, float *cent_ret, float *fwhm_ret) {

  //static float factor = 1.5f;
  static float factor = 2.5f;
  float x, y, d, sxy, sx2y, sigma;
  float oldcen, xm, flo, fhi, xhi, xlo;
  float area=0, cent=0, fwhm=0;
  int   maxw2, i, iteration, hi, lo;

  /* integrate peak over Cent +/- Factor*FWHM
     assume zero background
  */

  *area_ret = 0;
  *cent_ret = 0;
  *fwhm_ret = 0;

  lo = loch;
  hi = hich;
  oldcen = (float) (lo + hi) / 2;
  maxw2 = (hi-lo) * (hi-lo) / 10;
  if (hi - lo < 2) {
    lo = oldcen - 3;
    hi = lo + 6;
  }

  flo = 1;
  fhi = 1;
  for (iteration = 1; iteration < 20; ++iteration) {
    xm = oldcen;
    xlo = (float) lo + (1.0 - flo)/2.0 - xm;
    xhi = (float) hi - (1.0 - fhi)/2.0 - xm;

    area = his[lo] * flo + his[hi] * fhi;
    sxy = his[lo] * flo * xlo +
          his[hi] * fhi * xhi;
    sx2y = his[lo] * flo * xlo*xlo +
           his[hi] * fhi * xhi*xhi;

    for (i = lo + 1; i < hi; ++i) {
      x = (float) i - xm;
      y = his[i];
      area += y;
      sxy += y * x;
      sx2y += y * x*x;
    }
    if (area < 10) {
      //printf(" << %d lo, hi, area: %d %d %.0f\n", iteration, lo, hi, area);
      return 0;
    }

    cent = sxy / area;
    sigma = sx2y / area - cent*cent;
    cent += xm;
    if (sigma <= 0) sigma = 0.01;
    if (sigma > (float) maxw2) {
      if (VERBOSE)
        printf(" << Sigma too large: iteration %d; lo, hi, cent, sigma: %d %d %.1f %.1f\n",
               iteration, lo, hi, cent, sigma);
      return 0;
    }

    sigma = sqrtf(sigma);
    fwhm = sigma * 2.355;
    if (fabs(cent - oldcen) < fwhm/1000.0) break;

    oldcen = cent;
    d = factor * fwhm;
    if (d < 2) d = 2;
    lo = (int) (cent - d);
    hi = (int) (cent + d);
    if (lo < 1) lo = 1;
    flo = d - cent + (float) (lo + 1);
    fhi = d + cent - (float) hi;
  }

  if (VERBOSE) printf("find_cent: %d iterations\n", iteration);
  *area_ret = area;
  *cent_ret = cent;
  *fwhm_ret = fwhm;
  return 0;
} /* find_cent */

/* =======================================================================
autopeak()
finds and returns centroid of peak positioned around channel e in histogram his
          also calculates area and fwhm of the peak
 */

float autopeak(int *his, int e, float *area_ret, float *fwhm_ret) {

  float area, cent = 0.0f, fwhm, fwhm0, y;
  int   i, w2, hi, lo;

  *fwhm_ret = *area_ret = 0;
  fwhm0 = 2.5; // first rough guess
  
  i = e;
  y = his[i];
  if (y < 4) {
    printf("AutoPeak: Sorry, that peak is too weak! y = %.0f\n", y);
    return 0.0f;
  }
  w2 = lrint(fwhm0 * 3.0);
  lo = i - w2;
  hi = i + w2;
  if (lo < 1) lo = 1;
  find_cent(his, lo, hi, &area, &cent, &fwhm);

  if (VERBOSE) {
    printf("ap: %d; %d to %d, h = %.0f\n", e, lo, hi, y);
    printf("  >> pos, area, fwhm = %.1f, %.0f, %.2f\n", cent, area, fwhm);
  }

  if (area < 10) {
    if (VERBOSE) printf("AutoPeak: Sorry, could not identify the peak!\n");
    return 0.0f;
  }
  *area_ret = area;
  *fwhm_ret = fwhm;
  return cent;
} /* autopeak */

float autopeak1(int *his, int lo, int hi, float *area, float *fwhm) {

  float cent = 0.0f, fwhm0;
  int   i=lo, j, max=0, w2, lo1, hi1;

  *area = 0;
  fwhm0 = *fwhm; // first rough guess
  // first find mode for use as a first guess for the position
  for (j=lo; j<hi; j++) {
    if (his[j] > max) {
      max = his[j];
      i  = j;
    }
  }
  if (max < 4) {
    printf("AutoPeak: Sorry, that peak is too weak! max = %d\n", max);
    return 0.0f;
  }
  w2 = lrint(fwhm0 * 2.0);
  lo1 = i - w2;
  hi1 = i + w2;
  if (lo1 < lo) lo1 = lo;
  if (hi1 > hi) hi1 = hi;
  find_cent(his, lo1, hi1, area, &cent, fwhm);

  if (VERBOSE) {
    printf("autopk1: %d; %d to %d, h = %d\n", i, lo1, hi1, max);
    printf("  >> pos, area, fwhm = %.1f, %.0f, %.2f\n", cent, *area, *fwhm);
  }

  if (*area < 10) {
    if (VERBOSE) printf("AutoPeak: Sorry, could not identify the peak!\n");
    return 0.0f;
  }
  return cent;
} /* autopeak1 */


float autopeak2(int *his, int lo, int hi, float *area_ret, float *fwhm_ret) {

  double area=0, cent=0, fwhm=0;
  int    i;

  for (i=lo; i<hi; i++) {
    area += his[i];
    cent += his[i] * i;
  }
  if (area < 10) {
    if (VERBOSE) printf("AutoPeak2: Sorry, that peak is too weak! Area = %.0lf\n", area);
    return 0.0;
  }
  cent /= area;

  for (i=lo; i<hi; i++)
    fwhm += his[i] * (((double) i) - cent) * (((double) i) - cent);
  fwhm = sqrt(fwhm/area) * 2.355;

  *area_ret = area;
  *fwhm_ret = fwhm;
  return (float) cent;
} /* autopeak2 */

/* ---------------------------------------------------------------- */

/* ----- write signals to file, to look at with radware (gf3) ----- */
int write_sig(short *sig, int nsamples, int chan, FILE *file) {
  /* create four different sets of data structures, so that we can handle
     writing to up to four different files  */
  static struct {
    int nout, dptr, sptr, nchan[200];
    FILE *fp;
  } ws[4];            // info about each of the 4 different files
  static struct {
    int num;
    int nextptr;
    int dir[100][2];
  } dir[4];           // directory of spectra in the file
  static struct {
    int id;
    int mode;
    int xlen;
    int expand;
    char txt[64];
    int x0;
    int nx;
  } shead[4];         // spectrum header
  int i, id, fnum;
  static int first = 1;


  if (first) {
    for (i=0; i<4; i++) ws[i].fp = 0;
    first = 0;
  }

  /* look for data strutures matching this file
     if there is none, create one 
   */
  for (fnum=0; fnum<4; fnum++) {
    if (ws[fnum].fp == file) break;
    if (ws[fnum].fp == 0) {  // no match to existing structure; make a new one
      ws[fnum].fp = file;
      ws[fnum].nout = ws[fnum].dptr = ws[fnum].sptr = 0;
      memset(ws[fnum].nchan, 0, sizeof(ws[fnum].nchan));  // zero sigs out so far, for all chans
      dir[fnum].num = 0;                        // up to 100 sig IDs per directory
      dir[fnum].nextptr = -1;                   // pointer to next directory
      for (i=0; i<100; i++) dir[fnum].dir[i][0] = dir[fnum].dir[i][1] = -1;
      fwrite(&dir[fnum], sizeof(dir[fnum]), 1, file);
      ws[fnum].sptr = sizeof(dir[fnum]);        // pointer to where to write next spectrum
      break;
    }
  }
  if (fnum > 4) {
    printf("\nERROR in write_sig: Too many files used!\n\n");
    exit(-1);
  }

  ws[fnum].nout++;
  if (ws[fnum].nchan[chan] > 999) return 1;     // already wrote 1000 sigs for this channel
  id = 1000*chan + ws[fnum].nchan[chan]++;      // spectrum id depends on channel
  dir[fnum].dir[dir[fnum].num][0] = id;         // give the spectrum that complex id 
  dir[fnum].dir[dir[fnum].num][1] = ws[fnum].sptr;
  dir[fnum].num++;
  dir[fnum].dir[dir[fnum].num][0] = ws[fnum].nout-1; // but also give it a simple incremental id
  dir[fnum].dir[dir[fnum].num][1] = ws[fnum].sptr;
  fseek(file, ws[fnum].sptr, SEEK_SET);        // go to where we should write the spectrum
  shead[fnum].id = id;
  shead[fnum].mode = 1;                        // shorts
  shead[fnum].xlen = nsamples;
  shead[fnum].expand = 0;
  sprintf(shead[fnum].txt, "%.5d Ch %3d %d-sample signal", id, chan, nsamples);
  shead[fnum].x0 = 0;
  shead[fnum].nx = nsamples;
  if (VERBOSE) printf("write_sig %s to %d\n", shead[fnum].txt, ws[fnum].sptr);
  fwrite(&shead[fnum], sizeof(shead[fnum]), 1, file);
  fwrite(sig, sizeof(short), nsamples, file);
  ws[fnum].sptr = ftell(file);                // pointer to where to write next spectrum

  /* now re-write the directory */
  fseek(file, ws[fnum].dptr, SEEK_SET);       // go to the start of the directory
  dir[fnum].num++;                            // up to 100 sig IDs per directory
  dir[fnum].nextptr = -1;                     // pointer to next directory
  if (dir[fnum].num == 100)
    dir[fnum].nextptr = ws[fnum].dptr = ws[fnum].sptr;
  fwrite(&dir[fnum], sizeof(dir[fnum]), 1, file); // do wre-write

  if (dir[fnum].num == 100) {
    fseek(file, ws[fnum].dptr, SEEK_SET);     // go to the start of NEW directory
    dir[fnum].num = 0;
    dir[fnum].nextptr = -1;                   // pointer to next directory
    for (i=0; i<100; i++) dir[fnum].dir[i][0] = dir[fnum].dir[i][1] = -1;
    fwrite(&dir[fnum], sizeof(dir[fnum]), 1, file);
    ws[fnum].sptr += sizeof(dir[fnum]);       // pointer to where to write next spectrum
  }
  return 0;
} /* write_sig */


/* ----- write histograms to file, to look at with radware (gf3) ----- */
int write_his(int *his, int nch, int sp_id, char *sp_name, FILE *file) {

  static int nout = 0, dptr = 0, sptr = 0, id;
  static FILE *fsave = 0;
  static struct {
    int num;
    int nextptr;
    int dir[100][2];
  } dir;
  static struct {
    int id;
    int mode;
    int xlen;
    int expand;
    char txt[64];
  } shead;
  int i, j=0, k, dd[2];


  if (file != fsave || sp_id == 0) {  // new file; initialize
    nout = dptr = 0;
    dir.num = 0;                      // up to 100 sp IDs per directory
    dir.nextptr = -1;                 // pointer to next directory
    for (i=0; i<100; i++) dir.dir[i][0] = dir.dir[i][1] = -1;
    fwrite(&dir, sizeof(dir), 1, file);
    sptr = sizeof(dir);               // pointer to where to write next spectrum
  }
  nout++;
  fsave = file;

  /* find first non-zero channel */
  for (k=0; k<nch; k++) if (his[k] != 0) break;
  id = sp_id;
  dir.dir[dir.num][0] = id;
  dir.dir[dir.num][1] = sptr;

  fseek(file, sptr, SEEK_SET);        // go to where we should write the spectrum
  shead.id = id;
  shead.mode = 2;                     // ints
  shead.xlen = nch;
  shead.expand = 0;
  strncpy(shead.txt, sp_name, 64);
  // printf("write_his %s to %d\n", shead.txt, sptr);
  fwrite(&shead, sizeof(shead), 1, file);

  while (k < nch) {
    dd[0] = k;  // first non-zero channel
    dd[1] = nch-k;
    /* find last non-zero ch before a string of zeroes */
    for (i=k+1; i<nch; i++) {
      if (his[i-1] != 0 && his[i] == 0) j = i;
      if (his[i] == 0 && i >= j+5) {
        dd[1] = j-k;
        break;
      }
    }
    fwrite(dd, sizeof(int), 2, file);
    fwrite(his+k, sizeof(int), dd[1], file);
    /* find next non-zero channel */
    for (k=i; k<nch; k++) if (his[k] != 0) break;
  }
  dd[0] = dd[1] = -1;
  fwrite(dd, sizeof(int), 2, file);
        
  sptr = ftell(file);                 // pointer to where to write next spectrum


  /* now re-write the directory */
  fseek(file, dptr, SEEK_SET);        // go to the start of the directory
  dir.num++;                          // up to 100 sp IDs per directory
  dir.nextptr = -1;                   // pointer to next directory
  if (dir.num == 100)
    dir.nextptr = dptr = sptr;
  fwrite(&dir, sizeof(dir), 1, file); // do re-write

  if (dir.num == 100) {
    fseek(file, dptr, SEEK_SET);      // go to the start of NEW directory
    dir.num = 0;
    dir.nextptr = -1;                 // pointer to next directory
    for (i=0; i<100; i++) dir.dir[i][0] = dir.dir[i][1] = -1;
    fwrite(&dir, sizeof(dir), 1, file);
    sptr += sizeof(dir);              // pointer to where to write next spectrum
  }
  return 0;
} /* write_his */

/* ---------------------------------------- */

int data_clean_info_write(MJRunInfo *runInfo, DataClean *dcInfo) {

  int   i;
  char  fname[256], *c;
  FILE  *f_out = 0;

  c = runInfo->filename;
  while (strchr(c, '/')) c = strchr(c, '/') + 1;  // remove initial path
  if (!strncmp(c, "DS", 2)) {                     // input data file starts with "DS"
    sprintf(fname, "%s.dcl", c);
  } else {
    sprintf(fname, "run%d.dcl", runInfo->runNumber);
  }
  if ((f_out = fopen(fname, "r"))) {  // file already exists; modify fname
    fclose(f_out);
    strncat(fname, "_new", sizeof(fname)-5);
  }
  if (!(f_out = fopen(fname, "w"))) {
    printf("ERROR in data_clean_info_write(): Cannot open output file %s\n", fname);
    return 1;
  }

  fprintf(f_out,"#        (0.3ADC)  (0.3ADC) (0.06ADC)      (0.1ADC)  (0.1ADC) (0.02ADC)\n");
  fprintf(f_out,"#DetID  BLValue_HG  RMS_HG  slope_HG      BLValue_LG  RMS_LG  slope_LG\n");
  for (i=0; i < runInfo->nGe; i++) {
    if (dcInfo->blsl_lo[i] != -99 || dcInfo->blsl_lo[i+100] != -99)
      fprintf(f_out,
              "%4d %7d %4d %6d %6d %4d %9d %4d %6d %6d %4d %12.0f %6.0f\n",
              i, dcInfo->bl_lo[i], dcInfo->bl_hi[i], dcInfo->blrms_hi[i],
              dcInfo->blsl_lo[i], dcInfo->blsl_hi[i], dcInfo->bl_lo[i+100],
              dcInfo->bl_hi[i+100], dcInfo->blrms_hi[i+100], dcInfo->blsl_lo[i+100],
              dcInfo->blsl_hi[i+100], dcInfo->bl[i], dcInfo->bl[i+100]);
  }

  fclose(f_out);
  return 0;
}

/* ---------------------------------------- */

int data_clean_info_read(MJRunInfo *runInfo, DataClean *dcInfo) {

  int   i, default_read=0;
  char  fname[256], line[256], *c;
  FILE  *f_in = 0;

  /*
    read data cleaning limits from file
    first read from default.dcl, if it exists
      // then overwrite those values with data from a run-specific file
      then, if not, look for a run-specific file
        -look first for <data_filename>.dcl, then run<number>.dcl
  */

  if ((f_in = fopen("default.dcl", "r"))) {
    printf("\n Initializing data cleaning from file default.dcl\n\n");
    while (fgets(line, sizeof(line), f_in)) {
      if (line[0] == '#') continue;
      if (sscanf(line, "%5d", &i) != 1 ||
          i < 0 || i > 99 || 
          sscanf(line+5, "%d %d %d %d %d %d %d %d %d %d",
                 &dcInfo->bl_lo[i], &dcInfo->bl_hi[i], &dcInfo->blrms_hi[i],
                 &dcInfo->blsl_lo[i], &dcInfo->blsl_hi[i], &dcInfo->bl_lo[i+100],
                 &dcInfo->bl_hi[i+100], &dcInfo->blrms_hi[i+100], &dcInfo->blsl_lo[i+100],
                 &dcInfo->blsl_hi[i+100]) != 10) {
        printf(" ERROR; reading of data cleaning limits from file default.dcl failed!\n\n");
        fclose(f_in);
        return 1;
      }
    }
    fclose(f_in);
    default_read = 1;
    return 0;             // comment this line out to overwrite default.dcl values
  }

  // now look for run-specific dcl file
  c = runInfo->filename;
  while (strchr(c, '/')) c = strchr(c, '/') + 1;  // remove initial path
  sprintf(fname, "%s.dcl", c);
  if ((f_in = fopen(fname, "r"))) {
    printf("\n Initializing data cleaning from file %s\n\n", fname);
  } else {
    sprintf(fname, "run%d.dcl", runInfo->runNumber);
    if ((f_in = fopen(fname, "r"))) {
      printf("\n Initializing data cleaning from file %s\n\n", fname);
    } else {
      if (!default_read) printf("\n No data cleaning limits file; tried %s.dcl, %s, and default.dcl\n",
                       c, fname);
      return 1-default_read;
    }
  }

  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#') continue;
    if (sscanf(line, "%5d", &i) != 1 ||
        i < 0 || i > 99 || 
        sscanf(line+5, "%d %d %d %d %d %d %d %d %d %d",
               &dcInfo->bl_lo[i], &dcInfo->bl_hi[i], &dcInfo->blrms_hi[i],
               &dcInfo->blsl_lo[i], &dcInfo->blsl_hi[i], &dcInfo->bl_lo[i+100],
               &dcInfo->bl_hi[i+100], &dcInfo->blrms_hi[i+100], &dcInfo->blsl_lo[i+100],
               &dcInfo->blsl_hi[i+100]) != 10) {
    
      printf(" ERROR; reading of data cleaning limit file failed!\n\n");
      fclose(f_in);
      return 1;
    }
  }

  fclose(f_in);
  return 0;
} /* data_clean_info_read */

/* ---------------------------------------------------------------- */

int data_clean_init(MJRunInfo *runInfo, DataClean *dcInfo) {

  int i, j;

  /* initialize data-cleaning values to defaults */
  if (0) { // loose cuts
    for (i=0; i<100; i++) {
      dcInfo->bl_lo[i]    = -99;   //
      dcInfo->bl_hi[i]    = 699;   // these seem like conservative (loose) limits
      dcInfo->blrms_hi[i] =  79;   //      for HG chs with no pulser signals
      dcInfo->blsl_lo[i]  = -99;   //
      dcInfo->blsl_hi[i]  =  99;   //
      dcInfo->modified[i] =  0;
    }
    for (i=100; i<200; i++) {
      dcInfo->bl_lo[i]    = -399;  //
      dcInfo->bl_hi[i]    = 899;   // these seem like conservative (loose) limits
      dcInfo->blrms_hi[i] =  69;   //      for LG chs with no pulser signals
      dcInfo->blsl_lo[i]  = -99;   //
      dcInfo->blsl_hi[i]  =  99;   //
      dcInfo->modified[i] =  0;
    }
  } else { // tighter cuts against noise
    for (i=0; i<100; i++) {
      dcInfo->bl_lo[i]    = -99;   // 0.3 ADC
      dcInfo->bl_hi[i]    = 699;   // 0.3 ADC
      dcInfo->blrms_hi[i] =  20;   // 0.3 ADC
      dcInfo->blsl_lo[i]  = -20;   // 0.06 ADC
      dcInfo->blsl_hi[i]  =  20;   // 0.06 ADC
      dcInfo->modified[i] =  0;
    }
    for (i=100; i<200; i++) {
      dcInfo->bl_lo[i]    = -399;  // 0.1 ADC
      dcInfo->bl_hi[i]    = 899;   // 0.1 ADC
      dcInfo->blrms_hi[i] =  28;   // 0.1 ADC
      dcInfo->blsl_lo[i]  = -24;   // 0.02 ADC
      dcInfo->blsl_hi[i]  =  24;   // 0.02 ADC
      dcInfo->modified[i] =  0;
    }
  }

  j = data_clean_info_read(runInfo, dcInfo);
  if (VERBOSE) {
    printf(" DetID  BLValue_HG  RMS_HG  slope_HG      BLValue_LG  RMS_LG  slope_LG\n");
    for (i=0; i < runInfo->nGe; i++) {
      if (dcInfo->blsl_lo[i] != -99 || dcInfo->blsl_lo[i+100] != -99)
        printf("%4d %7d %4d %6d %6d %4d %9d %4d %6d %6d %4d\n",
               i, dcInfo->bl_lo[i], dcInfo->bl_hi[i], dcInfo->blrms_hi[i],
               dcInfo->blsl_lo[i], dcInfo->blsl_hi[i], dcInfo->bl_lo[i+100],
               dcInfo->bl_hi[i+100], dcInfo->blrms_hi[i+100], dcInfo->blsl_lo[i+100],
               dcInfo->blsl_hi[i+100]);
    }
  }
  return 1-j;
} /* data_clean_init */

/* ---------------------------------------------------------------- */

/*
  returns 0 for clean signals, or sum of:
          1 for bad mean baseline value (offset)
          2 for noisy RMS of baseline
          4 for bad initial baseline slope
          8 for bad signal tail slope
         16 for saturated signal
*/

int data_clean(short *signal, int chan, DataClean *dcInfo) {

  int    i, sl, dirty=0;
  double s1=0, s2=0;

  if (dcInfo->blsl_lo[chan] == -99) return 0;

#ifdef SHORT_BASELINE
  // Note: Need only 620 samples of signal baseline for this...
  for (i=20; i < 620; i++) {
    s1 += signal[i];
    s2 += (int) signal[i] * (int) signal[i];
  }
  s2 = sqrt((600.0*s2 - s1*s1))/60.0; // RMS, 0.1 ADC units
  s1 /= 60.0;                         // mean, 0.1 ADC units
  sl = (trap_fixed(signal, 20, 200, 200))/4;  // slope, 0.02-ADC units
#else
  // Note: Need at least 820 samples of signal baseline for this...
  for (i=20; i < 820; i++) {
    s1 += signal[i];
    s2 += (int) signal[i] * (int) signal[i];
  }
  s2 = sqrt((800.0*s2 - s1*s1))/80.0; // RMS, 0.1 ADC units
  s1 /= 80.0;                         // mean, 0.1 ADC units
  sl = (trap_fixed(signal, 20, 300, 200))/6;  // slope, 0.02-ADC units
#endif
  if (chan<100) {   // HG
    s1 /= 3.0;      // 0.3 ADC units
    s2 /= 3.0;
    sl /= 3;        // 0.06 ADC units
  }

  if (DEBUG &&
      signal[1500] > 4800 && signal[1500] < 5150 &&
      (chan == 23 || chan == 27 || chan == 51) &&
      (sl < dcInfo->blsl_lo[chan] || sl > dcInfo->blsl_hi[chan]))
    printf("++++ chan %2d  sl = %d  h = %d\n", chan, sl, signal[1500]);

  if (s1 < dcInfo->bl_lo[chan] || s1 > dcInfo->bl_hi[chan]) dirty = 1;
  if (s2 > dcInfo->blrms_hi[chan]) dirty += 2;
  if (sl < dcInfo->blsl_lo[chan] || sl > dcInfo->blsl_hi[chan]) dirty += 4;
  if (0)  printf("s1 s2 sl %.0lf %.0lf %d lims %d %d %d %d %d >> %d\n", s1,s2,sl, 
                 dcInfo->bl_lo[chan], dcInfo->bl_hi[chan], dcInfo->blrms_hi[chan],
                 dcInfo->blsl_lo[chan], dcInfo->blsl_hi[chan], dirty);

  if (signal[1300] > 8100 || signal[1400] > 8100 ||
      signal[1600] > 8100 || signal[1800] > 8100 ) { // overflow event
    dirty += 16;
  } else {
    /* make sure end of signal is sloping down */
    s1 = (trap_fixed(signal,  350, 300, 600))/6;  // ~ signal height
    s2 = (trap_fixed(signal, 1250, 300, 140))/6;  // final slope, 1250-1990
    if (chan < 100) {   // LG
      s1 /= 3.0;        // 0.3 ADC units
      s2 /= 3.0;
    }
    // FIXME: Factor 0.04 could be better optimized?
    if (s2 + 0.04*s1 > dcInfo->blsl_hi[chan]) dirty += 8;
  }

  return dirty;
} /* data_clean */

float autopeak3a(int *his, int lo, int hi, float *area_ret, float *fwhm_ret) {

  double area=0, cent=0, fwhm=0, factor = 2.6;
  int    i, j, w, imax=lo, amax=0;


  w = *fwhm_ret * factor + 0.5;  // integrate over +- 1.3 * fwhm ~ 3 sigma
  if (w > hi - lo) w = hi-lo;    //   or from lo to hi

  /* first find the w-bin-wide region within lo-hi with max # of counts */
  for (i = lo; i <= hi-w; i++) {
    area = 0;
    for (j=i; j<i+w; j++) area += his[j];
    if (area > 50 && area >= amax) {
      imax = i;
      amax = area;
    }
    if (amax && area < amax/2) break;
   }
  if (!amax) return 0.0;

  /* then integrate to get area, centroid, and FWHM = 2.355*RMS */
  area = 0;
  for (j=imax; j<imax+w; j++) {
    area += his[j];
    cent += his[j] * j;
  }
  cent /= area;
  for (j=imax; j<imax+w; j++)
    fwhm += his[j] * (((double) j) - cent) * (((double) j) - cent);
  fwhm = sqrt(fwhm/area) * 2.355;

  *area_ret = area;
  *fwhm_ret = fwhm;
  return (float) cent;
} /* autopeak3a */

float autopeak3(int *his, int lo, int hi, float *area, float *fwhm) {

  double pos = 0, fwhm0=0;
  int    j, max = 0, lo_save=lo, hi_save=hi;

  //printf("ap3: lo, hi, fwhm: %d %d %.0f", lo, hi, *fwhm);
  // first find maximum over search range and use that to restrict range a little
  for (j=lo; j<hi; j++) {
    if (his[j] > max) {
      max = his[j];
      pos = j;
    }
  }
  lo = pos - 4 * *fwhm;
  hi = pos + 4 * *fwhm;
  if (lo < lo_save) lo = lo_save;
  if (hi > hi_save) hi = hi_save;
  //printf(" >>  lo, hi: %d %d\n", lo, hi);

  // now do search using iterative integration
  for (j=0; j<10; j++) {
    //printf(" > j, lo, hi, fwhm: %d %d %d %f\n", j, lo, hi, *fwhm);
    if ((pos = autopeak3a(his, lo, hi, area, fwhm)) < 1) return 0.0;
    if (*fwhm > hi_save - lo_save) return 0.0;
    if (*fwhm == fwhm0) return pos;
    fwhm0 = *fwhm;
    lo = pos - 2.0f * *fwhm + 0.5;
    hi = pos + 2.0f * *fwhm + 0.5;
    if (lo < lo_save) lo = lo_save;
    if (hi > hi_save) hi = hi_save;
  }
 
  return pos;
} /* autopeak3 */

float autopeak4(int *his, int lo, int hi, float input_width, float *area, float *fwhm) {

  // this is like autopeak3, but allows user to specify final integration range
  // in terms of FWHM (using input_width)

  double pos = 0, fwhm0=0;
  int    j, max = 0;

  // first find maximum over search range and use that to restrict range a little
  for (j=lo; j<hi; j++) {
    if (his[j] > max) {
      max = his[j];
      pos = j;
    }
  }
  lo = pos - 1.5 * input_width * *fwhm;  // factor of 1.5 is dropped below
  hi = pos + 1.5 * input_width * *fwhm;

  // now do search using iterative integration
  for (j=0; j<10; j++) {
    //printf(" > j, lo, hi, fwhm: %d %d %d %f\n", j, lo, hi, *fwhm);
    if ((pos = autopeak3a(his, lo, hi, area, fwhm)) < 1) return 0.0;
    if (*fwhm == fwhm0) return pos;
    fwhm0 = *fwhm;
    lo = pos - input_width * *fwhm + 0.5;
    hi = pos + input_width * *fwhm + 0.5;
  }
 
  return pos;
} /* autopeak4 */

/* ---------------------------------------- */

int PZ_correct(short *signal, float *fsignal, int len, int chan, PZinfo *PZI){

  double decay, decay2;
  float  e1, e2=0, frac2;
  int    i, base=0;

  /* check we have everything we need */
  if (len < 200 || fsignal == 0 || chan < 0 || chan > 199) return 1;
  if (PZI->tau[chan] < 1.0) return 0;

  /* start at beginning to do PZ correction */
  decay  = 0.01 / PZI->tau[chan];        // fractional drop in height over one step
  decay2 = 0.01 / PZI->tau2[chan];
  for (i=20; i<40; i++) base += signal[i];
  base /= 20;
  if (base > PZI->baseline[chan]) base = PZI->baseline[chan];
  frac2 = PZI->frac2[chan];

  fsignal[0] = e1 = signal[0];
  for (i = 1; i < len; i++) {
    e1 += signal[i] - signal[i-1] + (signal[i-1]-base)*decay;
    e2 += signal[i] - signal[i-1] - e2*decay2;
    fsignal[i] = e1 - frac2*e2;
  }

  return  0;
} /* PZ_correct() */

/* ---------------------------------------- */

int PZ_fcorrect(float *fsignal, int len, int chan, PZinfo *PZI){

  double decay, decay2;
  float  base=0, e1, e2, e3=0, frac2;
  int    i;

  /* check we have everything we need */
  if (len < 200 || fsignal == 0 || chan < 0 || chan > 199) return 1;
  if (PZI->tau[chan] < 1.0) return 0;

  /* start at beginning to do PZ correction */
  decay  = 0.01 / PZI->tau[chan];        // fractional drop in height over one step
  decay2 = 0.01 / PZI->tau2[chan];
  for (i=20; i<40; i++) base += fsignal[i];
  base /= 20.0;
  if (base > PZI->baseline[chan]) base = PZI->baseline[chan];
  frac2 = PZI->frac2[chan];

  e1 = e2 = fsignal[0];
  for (i = 1; i < len; i++) {
    e1 += fsignal[i] - e2 + (e2-base)*decay;
    e3 += fsignal[i] - e2 - e3*decay2;
    e2 = fsignal[i];
    fsignal[i] = e1 - frac2*e3;
  }

  return  0;
} /* PZ_fcorrect() */

/* ---------------------------------------- */

int PZ_info_read(MJRunInfo *runInfo, PZinfo *PZI) {

  int    i, n;
  char   line[256];
  float  a, b, c, d, e;
  FILE   *f_in;

  /*
    read pole-zero correction data from file PZ.input
  */
  for (i=0; i<200; i++) {
    PZI->baseline[i] = 0;
    PZI->tau[i] = 72.5;     // this and three following values are typical values for MJD
    PZI->tau2[i] = 2.1;
    PZI->frac2[i] = 0.007;
    PZI->bl_rms[i] = 4;
  }

  if (!(f_in = fopen("PZ.input", "r"))) return 0;
  printf("\n Initializing pole-zero values from file PZ.input\n\n");
  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#' || strlen(line) < 4) continue;
    e = 4;
    if ((n=sscanf(line,
                  "%d %f %f %f %f %f",
                  &i,&a,&b,&c,&d,&e)) < 5 ||
        i < 0 || i > 199) {
      printf(" ERROR; reading of PZ input file failed! line: %s\n", line);
      fclose(f_in);
      return 1;
    }
    PZI->baseline[i] = a;
    PZI->tau[i]      = b;
    PZI->tau2[i]     = c;
    PZI->frac2[i]    = d;
    PZI->bl_rms[i]   = e;
  }
  fclose(f_in);

  if (VERBOSE) {
    printf(" DetID       BL    tau     tau2   frac2 BL_RMS"
           "         BL    tau     tau2   frac2 BL_RMS\n");
    for (i=0; i < runInfo->nGe; i++) {
      if (PZI->tau[i] > 1 || PZI->tau[100+i] > 1)
        printf("%4d %10.0f %7.2f %7.2f %7.4f %7.2f %10.0f %7.2f %7.2f %7.4f %7.2f\n",
               i, PZI->baseline[i], PZI->tau[i], PZI->tau2[i], PZI->frac2[i], PZI->bl_rms[i],
               PZI->baseline[100+i], PZI->tau[100+i], PZI->tau2[100+i], PZI->frac2[100+i], PZI->bl_rms[100+i]);
    }
  }

  return 1;
} /* PZ_info_read */

/* ---------------------------------------- */

int PZ_info_write(MJRunInfo *runInfo, PZinfo *PZI) {

  int   i;
  FILE  *f_out;

  /*
    write pole-zero correction data to file PZ.output
  */
  if (!(f_out = fopen("PZ.output", "w"))) {
    printf("\n ERROR: Cannot open file PZ.output for writing!\n");
    return 1;
  }
  printf("\n Writing pole-zero values to file PZ.output\n");

  for (i=0; i<200; i++) {
    if (i%100 == 0) fprintf(f_out, "# ch baseline  PZ-tau PZ-tau2 PZ-frac2  BaselineRMS\n");
    // # ch baseline  PZ-tau PZ-tau2 PZ-frac2  BaselineRMS
    // %3d %8.0fxxx %8.4fxxx %6.2fx %9.5fxxxx %9.2fxxxx
    //  11      -61  58.5892   2.10   0.01190      2.70

    if (i%100 < runInfo->nGe && PZI->tau[i] > 2)
      fprintf(f_out, "%3d %8.0f %8.4f %6.2f %9.5f %9.2f\n",
              i, PZI->baseline[i], PZI->tau[i], PZI->tau2[i], PZI->frac2[i], PZI->bl_rms[i]);
  }
  fclose(f_out);
  return 0;
} /* PZ_info_write */

/*  ------------------------------------------------------------ */

int CTC_info_read(MJRunInfo *runInfo, CTCinfo *CTC) {

  int    i, b, n;
  char   line[256];
  float  a, e;
  double f;
  FILE   *f_in;

  /*
    read charge-trapping correction data from file ctc.input
  */
  for (i=0; i<200; i++) {
    CTC->e_dt_slope[i]    = 1.0;
    CTC->e_lamda_slope[i] = 1.0;
    CTC->e_lamda_gain[i]  = 1.0;
    CTC->best_dt_lamda[i]  = 0;
  }
  if (strlen(CTC->ctc_fname) > 0) printf("Trying to open ctc_fname: %s\n", CTC->ctc_fname);
  if (strlen(CTC->ctc_fname) == 0 || !(f_in = fopen(CTC->ctc_fname, "r"))) {
    strncpy(CTC->ctc_fname, "ctc.input", sizeof(CTC->ctc_fname));
    f_in = fopen(CTC->ctc_fname, "r");
  }
  if (!f_in) return 0;
  printf("\n Initializing charge-trapping correction values from file %s\n\n", CTC->ctc_fname);
  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#' || strlen(line) < 4) continue;
    e = 0; f = 1; b = 1;
    if ((n=sscanf(line,
                  "%d %f %f %lf %d",
                  &i,&a,&e,&f,&b)) < 2 ||
        i < 0 || i > 199) {
      printf(" ERROR; reading of ctc.input file failed! line: %s\n", line);
      fclose(f_in);
      return 1;
    }
    CTC->e_dt_slope[i]    = a;
    CTC->e_lamda_slope[i] = e;
    CTC->e_lamda_gain[i]  = f;
    CTC->best_dt_lamda[i] = b;
  }
  fclose(f_in);

  if (VERBOSE) {
    for (i=0; i<200; i++) {
      if (i%100 == 0)
        printf("Chan  e_dt_slope | e_lamda_slope e_lamda_gain best_option\n");
      //           1        7.50 |          0.42    1.0023012          0
      //        %4dx  %10.2fxxxx | %13.2fxxxxxxx %12.8lfxxxxx %10dxxxxxx
      if (i%100 < runInfo->nGe)
        printf("%4d  %12.2f | %13.2f %12.8lf %10d\n",
               i, CTC->e_dt_slope[i], CTC->e_lamda_slope[i], CTC->e_lamda_gain[i], CTC->best_dt_lamda[i]);
    }
  }

  return 1;
} /* CTC_info_read */

/* ---------------------------------------- */

int CTC_info_write(MJRunInfo *runInfo, CTCinfo *CTC) {

  int   i;
  FILE  *f_out;

  /*
    write charge-trapping correction data to file ctc.output
  */
  if (!(f_out = fopen("ctc.output", "w"))) {
    printf("\n ERROR: Cannot open file ctc.output for writing!\n");
    return 1;
  }
  printf("\n Writing charge-trapping correction values to file ctc.output\n");

  for (i=0; i<200; i++) {
    if (i%100 == 0)
      fprintf(f_out,
              "#Chan  e_dt_slope e_lamda_slope e_lamda_gain best_option\n");
    //             1        7.50          0.42    1.0023012          0
    //         %5dxx  %10.2fxxxx %13.2fxxxxxxx %12.8lfxxxxx %10dxxxxxx
    if (i%100 < runInfo->nGe)
      fprintf(f_out, "%5d  %10.2f %13.2f %12.8lf %10d\n",
              i, CTC->e_dt_slope[i], CTC->e_lamda_slope[i], CTC->e_lamda_gain[i], CTC->best_dt_lamda[i]);
  }

  fclose(f_out);
  return 0;
} /* CTC_info_write */

/*  ------------------------------------------------------------ */

int PSA_info_read(MJRunInfo *runInfo, PSAinfo *PSA) {

  int    i, n, ret_value=0;
  char   line[256];
  float  a, b, c, d, e, f, g, h, s, q, t;
  int    j, k, l, m, o, lo, hi;
  unsigned int t0 = 0;
  FILE   *f_in;

  /*
    read A/E and DCR data from file psa.input
  */
  for (i=0; i<200; i++) {
    PSA->ae_dt_slope[i]    = 0;
    PSA->ae_e_slope[i]     = 0;
    PSA->ae_pos[i]         = 500;
    PSA->ae_cut[i]         = 500;
    PSA->dcr_dt_slope[i]   = 0;
    PSA->dcr_lim[i]        = 0;
    PSA->lamda_dt_slope[i] = 0;
    PSA->lamda_lim[i]      = 0;
    PSA->lq_dt_slope[i]    = 0;
    PSA->lq_lim[i]         = -1;
    PSA->ae_t0[i]          = 0;
    PSA->ae_t_slope[i]     = 0;
  }

  if ((f_in = fopen("psa.input", "r"))) {
    ret_value = 1;
    printf("\n Initializing PSA values from file psa.input\n\n");
    while (fgets(line, sizeof(line), f_in)) {
      if (line[0] == '#' || strlen(line) < 4) continue;
      e = 0; f = 1; s = 0; q = -1; t0=0;
      if ((n=sscanf(line,
                    "%d %f %f %f %f %f %f %f %f %f %f %u %e",
                    &i,&a,&b,&c,&d,&e,&f,&g,&h,&s,&q, &t0,&t)) < 9 ||
          i < 0 || i > 199) {
        printf(" ERROR; reading of psa.input file failed! line: %s\n", line);
        fclose(f_in);
        return 1;
      }
      PSA->ae_dt_slope[i]    = a;
      PSA->ae_e_slope[i]     = b;
      PSA->ae_pos[i]         = c;
      PSA->ae_cut[i]         = d;
      PSA->dcr_dt_slope[i]   = e;
      PSA->dcr_lim[i]        = f;
      PSA->lamda_dt_slope[i] = g;
      PSA->lamda_lim[i]      = h;
      PSA->lq_dt_slope[i]    = s;
      PSA->lq_lim[i]         = q;
      if (n > 12 && t0 > 0) {
        PSA->ae_t0[i]        = t0;
        PSA->ae_t_slope[i]   = t;
      }
    }
    fclose(f_in);

    if (VERBOSE) {
      for (i=0; i<200; i++) {
        if (i%100 == 0)
          printf("Chan  ae_dt_slope ae_e_slope  ae_pos  ae_cut | dcr_dt_slope dcr_lim | lamda_dt_slope lamda_lim | lq_dt_slope  lq_lim\n");
        //           1         7.50      10.70 1007.50 1005.60 |         0.42    1.00 |           0.42      1.00 |   -1.0
        //        %4dx  %11.2fxxxxx %10.2fxxxx %7.2fxx %7.2fxx | %12.2fxxxxxx %7.2fxx | %14.2fxxxxxxxx %9.2fxxxx | %6.1fx
        if (i%100 < runInfo->nGe)
          printf("%4d  %11.2f %10.2f %7.2f %7.2f | %12.2f %7.2f | %14.2f %9.2f | %12.2f %6.1f\n",
                 i, PSA->ae_dt_slope[i], PSA->ae_e_slope[i], PSA->ae_pos[i], PSA->ae_cut[i],
                 PSA->dcr_dt_slope[i], PSA->dcr_lim[i], PSA->lamda_dt_slope[i], PSA->lamda_lim[i], PSA->lq_dt_slope[i], PSA->lq_lim[i]);
      }
    }
  }

  /*
    read individual trap filter parameters from data from file filters.input
  */
  for (i=0; i<200; i++) {
    PSA->e_ctc_rise[i] = TRAP_RISE;
    PSA->e_ctc_flat[i] = TRAP_FLAT;
    PSA->t0_rise[i]    = T0_TRAP_RISE;
    PSA->t0_thresh[i]  = T0_TRAP_THRESH;
    PSA->a_e_rise[i]   = A_E_RISE;
    PSA->a_e_factor[i] = A_E_FACTOR;
    PSA->gerda_aoe[i]  = GERDA_AoE;
  }

  if (!(f_in = fopen("filters.input", "r"))) return ret_value;
  printf("\n Initializing trapezoid values from file filters.input\n\n");
  while (fgets(line, sizeof(line), f_in)) {
    if (line[0] == '#' || strlen(line) < 4) continue;
    o = 0;
    if ((sscanf(line,
                  "%d %d %d %d %d %d %d %f %d",
                  &lo,&hi,&j,&k,&l,&m,&n,&a,&o)) < 8 ||
        lo < 0 || hi > 199 || lo > hi) {
      printf(" ERROR; reading of filters.input file failed! line: %s\n", line);
      fclose(f_in);
      return ret_value;
    }
    for (i=lo; i<=hi; i++) {
      PSA->e_ctc_rise[i] = j;
      PSA->e_ctc_flat[i] = k;
      PSA->t0_rise[i]    = l;
      PSA->t0_thresh[i]  = m;
      PSA->a_e_rise[i]   = n;
      PSA->a_e_factor[i] = a;
      PSA->gerda_aoe[i]  = o;
    }
    if (hi < 100) {
      for (i=lo+100; i<=hi+100; i++) {    // also use these values as the defaults for LG channels
        PSA->e_ctc_rise[i] = j;
        PSA->e_ctc_flat[i] = k;
        PSA->t0_rise[i]    = l;
        PSA->t0_thresh[i]  = m;
        PSA->a_e_rise[i]   = n;
        PSA->a_e_factor[i] = a;
        PSA->gerda_aoe[i]  = o;
      }
    }
  }
  fclose(f_in);

  if (VERBOSE) {
    printf("Chan | e_ctc_rise e_ctc_flat  t0_rise t0_thresh ｜ a_e_rise a_e_factor  gerda_aoe\n");
    for (i=0; i<runInfo->nGe; i++) {
      printf("%4d ｜ %10d %10d  %7d %9d | %8d %10.1f  %9d\n",
             i, PSA->e_ctc_rise[i], PSA->e_ctc_flat[i], PSA->t0_rise[i], PSA->t0_thresh[i],
             PSA->a_e_rise[i], PSA->a_e_factor[i], PSA->gerda_aoe[i]);
    }
  }

  return 2;
} /* PSA_info_read */

/* ---------------------------------------- */

int PSA_info_write(MJRunInfo *runInfo, PSAinfo *PSA) {

  int   i;
  FILE  *f_out;

  /*
    write A/E and DCR data to file psa.output
  */
  if (!(f_out = fopen("psa.output", "w"))) {
    printf("\n ERROR: Cannot open file psa.output for writing!\n");
    return 1;
  }
  printf("\n Writing PSA values to file psa.output\n");

  for (i=0; i<200; i++) {
    if (i%100 == 0)
      fprintf(f_out,
              "#Chan  ae_dt_slope ae_e_slope  ae_pos  ae_cut  dcr_dt_slope dcr_lim  lamda_dt_slope lamda_lim  lq_dt_slope lq_lim\n");
      //           1         7.50      10.70 1007.50 1005.60          0.42    1.00            0.42      1.00         3.25    123.0
      //       %5dxx  %11.2fxxxxx %10.2fxxxx %7.2fxx %7.2fxx  %12.2fxxxxxx %7.2fxx  %14.2fxxxxxxxx %9.2fxxxx %12.2fxxxxxxx %6.1fxx
    if (i%100 < runInfo->nGe) {
      fprintf(f_out, "%5d  %11.2f %10.2f %7.2f %7.2f  %12.2f %7.2f  %14.2f %9.2f %12.2f %6.1f",
              i, PSA->ae_dt_slope[i], PSA->ae_e_slope[i], PSA->ae_pos[i], PSA->ae_cut[i],
              PSA->dcr_dt_slope[i], PSA->dcr_lim[i], PSA->lamda_dt_slope[i], PSA->lamda_lim[i], PSA->lq_dt_slope[i], PSA->lq_lim[i]);
      if (PSA->ae_t0[i] && fabs(PSA->ae_t_slope[i]) > 0.00001) {
        fprintf(f_out, " %12u %10.3e\n", PSA->ae_t0[i], PSA->ae_t_slope[i]);
      } else {
        fprintf(f_out, "\n");
      }
    }
  }

  fclose(f_out);
  return 0;
} /* PSA_info_write */

/*  ------------------------------------------------------------ */

int checkGranularity(MJDetInfo *Dets, MJRunInfo *runInfo,
                     int nChData, BdEvent *ChData[]) {
  
  int      i, ievt, chan;
  int      chan_list[200];     // list of hit detectors/channels in the event
  int      nchan=0, ndet=0;
  
  /* LOOP over channel-events in the built-event */
  for (ievt = 0; ievt < nChData; ievt++) {
    chan_list[nchan] = -1;
    if (ChData[ievt]->orca_type == runInfo->dataIdGM ||
        ChData[ievt]->orca_type == runInfo->dataIdGA) {
      if (ChData[ievt]->chan < 0 ||
          ChData[ievt]->chan >= 100 + runInfo->nGe) continue;
      chan_list[nchan++] = chan = ChData[ievt]->chan;

      for (i=0; i<nchan-1; i++)
        if (chan_list[i]%100 == chan%100) break;  // already saw this detector
      if (i == nchan-1) ndet++;   // this detector _not_ yet seen
    }
  }
  // if (ndet > 2) printf(" nChData, ndet: %d %d\n", nChData, ndet);
  if (ndet <= 0) return 0;
  return ndet-1;
} /* checkGranularity */

/*  ------------------------------------------------------------ */

int compress_signal(short *sig_in, unsigned short *sig_out, int sig_len_in) {

  int   i, j, max1, max2, min1, min2, ds, nb1, nb2;
  int   iso, nw, bp, dd1, dd2;
  unsigned short db[2];
  unsigned int   *dd = (unsigned int *) db;
  static unsigned short mask[17] = {0, 1,3,7,15, 31,63,127,255,
                                    511,1023,2047,4095, 8191,16383,32767,65535};

  //static int len[17] = {4096, 2048,512,256,128, 128,128,128,128,
  //                      128,128,128,128, 48,48,48,48};
  /* ------------ do compression of signal ------------ */
  j = iso = bp = 0;

  sig_out[iso++] = sig_len_in;     // signal length
  while (j < sig_len_in) {         // j = starting index of section of signal
    // find optimal method and length for compression of next section of signal 
    max1 = min1 = sig_in[j];
    max2 = -16000;
    min2 = 16000;
    nb1 = nb2 = 2;
    nw = 1;
    for (i=j+1; i < sig_len_in && i < j+48; i++) { // FIXME; # 48 could be tuned better?
      if (max1 < sig_in[i]) max1 = sig_in[i];
      if (min1 > sig_in[i]) min1 = sig_in[i];
      ds = sig_in[i] - sig_in[i-1];
      if (max2 < ds) max2 = ds;
      if (min2 > ds) min2 = ds;
        nw++;
    }
    if (max1-min1 <= max2-min2) { // use absolute values
      nb2 = 99;
      while (max1 - min1 > mask[nb1]) nb1++;
      //for (; i < sig_len_in && i < j+len[nb1]; i++) {
      for (; i < sig_len_in && i < j+128; i++) { // FIXME; # 128 could be tuned better?
        if (max1 < sig_in[i]) max1 = sig_in[i];
        dd1 = max1 - min1;
        if (min1 > sig_in[i]) dd1 = max1 - sig_in[i];
        if (dd1 > mask[nb1]) break;
        if (min1 > sig_in[i]) min1 = sig_in[i];
        nw++;
      }
    } else {                      // use difference values
      nb1 = 99;
      while (max2 - min2 > mask[nb2]) nb2++;
      //for (; i < sig_len_in && i < j+len[nb1]; i++) {
      for (; i < sig_len_in && i < j+128; i++) { // FIXME; # 128 could be tuned better?
        ds = sig_in[i] - sig_in[i-1];
        if (max2 < ds) max2 = ds;
        dd2 = max2 - min2;
        if (min2 > ds) dd2 = max2 - ds;
        if (dd2 > mask[nb2]) break;
        if (min2 > ds) min2 = ds;
        nw++;
      }
    }

    if (bp > 0) iso++;
    /*  -----  do actual compression  -----  */
    sig_out[iso++] = nw;  // compressed signal data, first byte = # samples
    bp = 0;               // bit pointer
    if (nb1 <= nb2) {
      /*  -----  encode absolute values  -----  */
      sig_out[iso++] = nb1;                    // # bits used for encoding
      sig_out[iso++] = (unsigned short) min1;  // min value used for encoding
      for (i = iso; i <= iso + nw*nb1/16; i++) sig_out[i] = 0;
      for (i = j; i < j + nw; i++) {
        dd[0] = sig_in[i] - min1;              // value to encode
        dd[0] = dd[0] << (32 - bp - nb1);
        sig_out[iso] |= db[1];
        bp += nb1;
        if (bp > 15) {
          sig_out[++iso] = db[0];
          bp -= 16;
        }
      }

    } else {
      /*  -----  encode derivative / difference values  -----  */
      sig_out[iso++] = nb2 + 32;  // # bits used for encoding, plus flag
      sig_out[iso++] = (unsigned short) sig_in[j];  // starting signal value
      sig_out[iso++] = (unsigned short) min2;       // min value used for encoding
      for (i = iso; i <= iso + nw*nb2/16; i++) sig_out[i] = 0;
      for (i = j+1; i < j + nw; i++) {
        dd[0] = sig_in[i] - sig_in[i-1] - min2;     // value to encode
        dd[0]= dd[0] << (32 - bp - nb2);
        sig_out[iso] |= db[1];
        bp += nb2;
        if (bp > 15) {
          sig_out[++iso] = db[0];
          bp -= 16;
        }
      }
    }
    j += nw;
  }

  if (bp > 0) iso++;
  if (iso%2) iso++;     // make sure iso is even for 4-byte padding
  return iso;           // number of shorts in compressed signal data

} /* compress_signal */

/*  ------------------------------------------------------------ */

int decompress_signal(unsigned short *sig_in, short *sig_out, int sig_len_in) {

  int   i, j, min, nb, isi, iso, nw, bp, siglen;
  unsigned short db[2];
  unsigned int   *dd = (unsigned int *) db;
  static unsigned short mask[17] = {0, 1,3,7,15, 31,63,127,255,
                                    511,1023,2047,4095, 8191,16383,32767,65535};

  /* ------------ do decompression of signal ------------ */
  j = isi = iso = bp = 0;
  siglen = (short) sig_in[isi++];  // signal length
  //printf("<<< siglen = %d\n", siglen);
  for (i=0; i<2048; i++) sig_out[i] = 0;
  while (isi < sig_len_in && iso < siglen) {
    if (bp > 0) isi++;
    bp = 0;              // bit pointer
    nw = sig_in[isi++];  // number of samples encoded in this chunk
    nb = sig_in[isi++];  // number of bits used in compression

    if (nb < 32) {
      /*  -----  decode absolute values  -----  */
      min = (short) sig_in[isi++];  // min value used for encoding
      db[0] = sig_in[isi];
      for (i = 0; i < nw && iso < siglen; i++) {
        if (bp+nb > 15) {
          bp -= 16;
          db[1] = sig_in[isi++];
          db[0] = sig_in[isi];
          dd[0] = dd[0] << (bp+nb);
        } else {
          dd[0] = dd[0] << nb;
        }
        sig_out[iso++] = (db[1] & mask[nb]) + min;
        bp += nb;
      }

    } else {
      nb -= 32;
      /*  -----  decode derivative / difference values  -----  */
      sig_out[iso++] = (short) sig_in[isi++];  // starting signal value
      min = (short) sig_in[isi++];             // min value used for encoding
      db[0] = sig_in[isi];
      for (i = 1; i < nw && iso < siglen; i++) {
        if (bp+nb > 15) {
          bp -= 16;
          db[1] = sig_in[isi++];
          db[0] = sig_in[isi];
          dd[0] = dd[0] << (bp+nb);
        } else {
          dd[0] = dd[0] << nb;
        }
        sig_out[iso] = (db[1] & mask[nb]) + min + sig_out[iso-1]; iso++;
        bp += nb;
      }
    }
    j += nw;
  }

  if (siglen != iso) {
    printf("ERROR in decompress_signal: iso (%d ) != siglen (%d)!\n",
           iso, siglen);
  }
  return siglen;       // number of shorts in decompressed signal data

} /* decompress_signal */

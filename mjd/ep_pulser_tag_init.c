/*
   special version of eventprocess()

   event processing code for MJD
   - up to 16 digitizer modules for the Ge detectors, plus
   - QDCs etc for the veto

   This version just for extracting the pulser peak amplitudes and times from a background run.
   It uses two passes through the data:
      pass 1: extract pulser amplitudes
      pass 2: extract pulser period and t0

   David Radford   Dec 2016
*/
/* Process the global events that have been built in eventbuild()
   - extract energy, PSA cuts, etc
   - create histograms
   returns: -2 after finalization of first pass, to ask for second pass through data
            -1 on error,
             1 on pulser event,
             2 on bad event (bad = cannot be processed; not dirty)
             0 otherwise (good event, can be clean or dirty)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MJDSort.h"

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define DEBUG      0
#define VERBOSE    0

/* set USE_CC_TIMING to 1 to use common controller card pulser times for all
     detectors on that that card, or 0 to use individual detector times */
#define USE_CC_TIMING 1

/*  ------------------------------------------------------------ */


int eventprocess(MJDetInfo *Dets, MJRunInfo *runInfo, int nChData, BdEvent *ChData[]) {

  int    i, j, k, t, bad=0;
  int    ievt, idet, ch, crate, slot, chan=0;
  long long int  time;

  static int   module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int   det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int   chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  static int   totevts=0, builtevts=0, badevts=0, goodevts = 0;
  static int   *his[400];
  static int   pulser_period[200];
  static float area[200], d;
  // static float pos2[200], area2[200], fwhm2[200];
  static long long int last_pulser_t[200], dt;
  static int   hitchan[200];  // list of detectors/channels meeting energy gate

  static int first = 1;

  int    e, e_trapmax, CCnum;
  char   spname[256];
  FILE   *f_out;
  int    pass = runInfo->analysisPass;
  int    nchan = 0, ndet = 0;
  PTag   pt;


  /* --------------- initialization --------------- */
  if (first) {
    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
    /* malloc and clear histogram space */
    if ((his[0] = calloc(400*8192, sizeof(int))) == NULL) {
      printf("ERROR in eventprocess(); cannot malloc his!\n");
      exit(-1);
    }
    for (i=1; i<200; i++) {
      pt.pt0[i] = 0;
      pt.pdt[i] = 0;
      pt.pos[i] = pt.fwhm[i] = area[i] = 0.0;
      // pos2[i] = area2[i] = fwhm2[i] = 0.0;
      pulser_period[i] = Dets[0].pulseHighTime + Dets[0].pulseLowTime;  // default
    }
    for (i=1; i<runInfo->nGe; i++) {
      pulser_period[i] = pulser_period[i + 100] =
        Dets[i].pulseHighTime + Dets[i].pulseLowTime;
      pulser_period[Dets[i].CCnum+80] = pulser_period[i];  // common CC channels
    }
    for (i=1; i<400; i++)
      his[i] = his[i-1]+8192;
    // runInfo->minChTimeDiff = 200000;  // set to 2 ms to remove some noise signals
    first = 0;
  }

  /* ------------ finalization / end-of-run handling ------------- */
  if (nChData == -99) {  // special value as flag to trigger finalization
    if (pass == 0) {
      /* -1-1-1-1-1-1-1-1- first pass, to get pulser amplitudes -1-1-1-1-1-1-1-1- */

      /*  extract pulser peak positions  */
      if (VERBOSE) printf(" Chan       pos   area   fwhm\n");
      for (chan=0; chan < 100 + runInfo->nPT + runInfo->nGe; chan++) {
        if (chan > runInfo->nGe && chan < 100) continue;  // empty channel numbers

        /* find highest-energy peak */
        if (!(e = peak_find(his[chan], 45, 8000)))  continue; // no peak found
        pt.pos[chan] = autopeak(his[chan], e, &area[chan], &pt.fwhm[chan]);

        /*  find next-highest-energy peak  */
        /*
        if (!(e = peak_find(his[chan], 3, e/2)) ||
            (pos2[chan] = autopeak(his[chan], e, &area2[chan], &fwhm2[chan])) < 1.0) {
          // no second peak found; report results only for first peak
          printf("%5d %9.1f %6.0f %6.1f\n", chan, pt.pos[chan], area[chan], pt.fwhm[chan]);
          continue;
        }
        */
        /*  report full results  */
        if (VERBOSE) printf("%5d %9.1f %6.0f %6.1f\n",
                            chan, pt.pos[chan], area[chan], pt.fwhm[chan]);
      }
      if (VERBOSE) {
        printf("\n Chan  HG_area LG_area   diff\n");
        for (chan=0; chan<100; chan++) {
          if (area[chan] > 100)
            printf("%5d %8.1f %7.1f %6.1f\n",
                   chan, area[chan], area[chan+100], area[chan]-area[chan+100]);
        }
      }
      printf("\n");

      /* do some initialization to get ready for pass 2 */
      totevts = builtevts = badevts = goodevts = 0;
      for (i=0; i<200; i++) {
        last_pulser_t[i] = -1;
        if (pt.pos[i] > 3) {
          d = 4.0*pt.fwhm[i] + pt.pos[i]/1000;  // max delta-e allowed to meet energy gate
          if (pt.pos[i] > 40 && d < 4.0) d = 4;
          pt.elo[i] = pt.pos[i] - d;
          pt.ehi[i] = pt.pos[i] + d + 0.5;
        } else {
          pt.elo[i] = 1;
          pt.ehi[i] = 0;  // pt.ehi < pt.elo cuts all events
        }
      }

      return -2; // return -2 to ask for next pass through data
    }
    /* -1-1-1-1-1-1- END of first pass, to get pulser amplitudes -1-1-1-1-1-1- */
    /* -2-2-2-2-2-2-2-2- second pass, to get pulser timing -2-2-2-2-2-2-2-2- */

    /* save histograms */
    if ((f_out = fopen("ptag.rms", "w"))) {
      for (i=0; i<400; i++) {
        if (i < 200) {
          sprintf(spname, "%d; PT_init run %d; ch %d E_trap_max [ADC]",
                  i, runInfo->runNumber, i);
        } else if (i > 279 && i < 290) {
          sprintf(spname, "%d; PT_init run %d; CCard %d delta-t [us, 2 windows]",
                  i, runInfo->runNumber, i-280);
        } else if (i == 290) {
          sprintf(spname, "%d; PT_init run %d; PTagChs delta-t [us, 2 windows]",
                  i, runInfo->runNumber);
        } else {
          //sprintf(spname, "%d; PT_init run %d; chan %d delta-time [1.25 ms]",
          sprintf(spname, "%d; PT_init run %d; ch %d delta-t [us, 2 windows]",
                  i, runInfo->runNumber, i-200);
        }
        write_his(his[i], 8192, i, spname, f_out);
      }
      fclose(f_out);
    }

    /*  extract pulser time-peak positions  */
    for (chan=0; chan<200; chan++) {
      if (pt.elo[chan] > pt.ehi[chan] &&   // no pulser known for those chs
          (chan<80 || chan>87)) continue;  // but include common controller card times (chs 80-87)
      /* find highest peak */
      if ((t = peak_find(his[chan+200], 100, 8000))) {  // peak found
        pt.pdt[chan] = t + pulser_period[chan]*2 - 6000;
        pt.pt0[chan] -= 200*pt.pdt[chan]; // give a little extra margin (2 periods) for early pulser events
      }
    }
    if (USE_CC_TIMING) {  // use common controller card pulser times where possible
      for (CCnum = 0; CCnum<8; CCnum++) {
        if (pt.pdt[CCnum+80] > 0) {
          for (i=0; i<runInfo->nGe; i++) {
            if (Dets[i].CCnum == CCnum) {
              pt.pdt[i]     = pt.pdt[CCnum+80];
              pt.pdt[i+100] = pt.pdt[CCnum+80];
              pt.pt0[i]     = pt.pt0[CCnum+80];
              pt.pt0[i+100] = pt.pt0[CCnum+80];
            }
          }
        }
      }
    }

    if (!DS0 && SPECIAL_DET_16 &&
        pt.elo[16] > pt.ehi[16] && pt.elo[116] > pt.ehi[116]) {
      // special hack for very-low-amplitude pulser in ch 16 & 116; cross-talk from CC # 1
      Dets[16].CCnum = 1; // not the true CC number, but the one that gives the largest pulser signal
      pt.pos[16]  = 14.7;
      pt.fwhm[16] = 0.6;
      pt.pos[116]  = 4.5;
      pt.fwhm[116] = 0.3;
      area[16] = area[116] = 999;
      pt.elo[16] = 12;
      pt.ehi[16] = 17;
      pt.elo[116] = 3;
      pt.ehi[116] = 5;
      pt.pdt[16] = pt.pdt[116] = pt.pdt[81];
      pt.pt0[16] = pt.pt0[116] = pt.pt0[81];
    }

    /* report results and save pulser-tag results in a file */
    pulser_tag_info_write(runInfo, &pt, 1);

    /* report any problem channels */
    for (chan=0; chan<200; chan++) {
      if (pt.elo[chan] < pt.ehi[chan] && pt.pdt[chan] < 1) {
        if (chan%100 < runInfo->nGe)
          printf("\n>>>> NOTE that pulser period not found for chan %d (%s %s)\n"
                 ">>>>  This channel should be excluded from this data subset!\n",
                 chan, Dets[chan%100].DetName, Dets[chan%100].StrName);
        else 
          printf("\n>>>> NOTE that pulser period not found for pulser tag chan %d\n"
                 ">>>>  This channel should be excluded from this data subset!\n",
                 chan);
      }
    }

    return 0;
    /* -2-2-2-2-2-2-2- END of second pass, to get pulser timing -2-2-2-2-2-2-2- */
  }

  /* ------------- process input events  ------------- */

  builtevts++;
  totevts += nChData;

  // fillEvent: fill out entris in BdEvent ChData; chan, det, e, sig...
  fillEvent(runInfo, nChData, ChData, module_lu, det_lu, chan_lu);

  /* LOOP over channel-events in the built-event again to handle data */
  for (ievt = 0; ievt < nChData; ievt++) {
    crate = ChData[ievt]->crate;
    slot  = ChData[ievt]->slot;
    ch    = ChData[ievt]->ch;
    time  = ChData[ievt]->time;
    idet  = ChData[ievt]->det;
    chan  = ChData[ievt]->chan;
    if (DEBUG) {
      printf("ievt, len, cr, slot, time: %d %d %d %d %lld\n",
             ievt, ChData[ievt]->evlen, crate, slot, time);
      fflush(stdout);
    }
    if (ChData[ievt]->orca_type != runInfo->dataIdGM &&
        ChData[ievt]->orca_type != runInfo->dataIdGA) continue;  // only interested in Ge data

    /* identify det and chan IDs, and skip some types of bad events */
    if (ChData[ievt]->mod < 0 || ch > 9) {  // unknown module
      // printf("bad module %d %d\n", crate, slot);
      bad = 1;
    } else {
      if (chan > 199) continue;      // not an entirely bad event, but not data we want
      if (idet < 0 || chan < 0 ||                       //  unknown detector or channel
          (idet < runInfo->nGe &&
           ((chan < 99 && !Dets[idet].HGChEnabled)   || //  disabled channel
            (chan > 99 && !Dets[idet].LGChEnabled))) || //  disabled channel
          time < 0) {                                   //  timestamp out of range
        // printf("bad evt %d %d %lld\n", idet, chan, time);
        bad = 1;
      }
    }
    if (bad) {
      badevts++;
      return 2;
    }

    /* ---------------- Ge detector data processing ----------------- */
    /* at this point, we have good Gretina4M data */

    e_trapmax = ChData[ievt]->e;            // offline trap_max energy

    if (pass == 0) {   // ------------- first analysis pass
      if (e_trapmax > 0 && e_trapmax < 8192) his[chan][e_trapmax]++;

    } else {           // ------------- second analysis pass
      k = chan%100;    // detector number
      CCnum = 10;                                   // all pulser tag chs
      if (k < runInfo->nGe) CCnum = Dets[k].CCnum;  // Ge CC number

      /* histogram time since last big pulser event */
      if (last_pulser_t[chan] >= 0) {
        dt = (time - last_pulser_t[chan])/100;   // 1 us / bin

        if (0) {               // histogram time correlations between pulser tags
          if ((j = dt/1250) < 8192) {
            his[200+chan][j]++;      // 1.25 ms/bin
            his[280+CCnum][j]++;
          }

        } else if (0 &&       // histogram _only_ time between two real pulses
                   e_trapmax >= pt.elo[chan] &&
                   e_trapmax <= pt.ehi[chan]) {
          if ((j = dt/1250) < 8192) {
            his[200+chan][j]++;      // 1.25 ms/bin
            his[280+CCnum][j]++;
          }


        } else {          // histogram time since last pulser, but only parts of range, at 1 us/bin 
          j = 2000 + dt - pulser_period[chan];
          if (j > 0 && j < 4000) {
            his[200+chan][j]++;
            his[280+CCnum][j]++;  // common CC chs
          }
          j = 6000 + dt - pulser_period[chan] *2;
          if (j > 4000 && j < 8000) {
            his[200+chan][j]++;
            his[280+CCnum][j]++;  // common CC chs
          }
        }
      }

      if (e_trapmax >= pt.elo[chan] && e_trapmax <= pt.ehi[chan]) {  // new big pulser event
        // signal meets pulser energy gate
        last_pulser_t[chan] = time;
        // add it to a list of hit detectors
        for (i=0; i<nchan; i++)
          if (hitchan[i]%100 == chan%100) break;  // aleady saw this detector
        if (i == nchan) ndet++;
        hitchan[nchan++] = chan;        // save chan in the list of chs that met energy gate
      }
    }
    /* ------------- end of Ge detector data processing -------------- */

  } // END of loop over channel-events

  if (DS0) {
    for (i=1; i< nchan; i++) {
      for (j=0; j<i; j++) {
        his[50+hitchan[i]][hitchan[j]]++;
        his[50+hitchan[j]][hitchan[i]]++;
      }
    }
  }

  /* if at least two different detectors, or a detector and a pulser tag channel,
     are hit with the correct energy, then this must really be a pulser event */
  if (pass > 0 &&   // second pass
      ndet > 1) {   // multiple detectors met energy gate
    for (j=0; j<nchan; j++) {
      //save absolute time stamps of first pulser event for each detector, chan, CC
      chan = hitchan[j];
      idet = chan%100;                                    // detector number
      CCnum = 10;                                         // all pulser tag chs
      if (idet < runInfo->nGe) CCnum = Dets[idet].CCnum;  // Ge CC number
      if (pt.pt0[idet] == 0) pt.pt0[idet] = time;
      if (pt.pt0[chan] == 0) pt.pt0[chan] = time;
      if (pt.pt0[CCnum+80] == 0) pt.pt0[CCnum+80] = time;
    }
    return 1;
  }

  goodevts++;
  if (totevts % 20000 < (totevts-nChData) % 20000) {
    printf(" Processed %d builtevts in, %d %d good/bad evts out (%d%%)\n",
           builtevts, goodevts, badevts, (badevts*100)/builtevts);
  }
  return  0;
} /* eventprocess() */

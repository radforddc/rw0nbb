/*
   special version of eventprocess()

   event processing code for MJD
   - up to 16 digitizer modules for the Ge detectors, plus
   - QDCs etc for the veto

   This version compares hi-gain and low-gain channels to look for ADC nonlinearity.

   David Radford   June 2017
*/
/* Process the global events that have been built in eventbuild()
   - extract ratios
   - create histograms
   returns: -1 on error,
             2 on bad event (bad = cannot be processed; not dirty)
             0 otherwise (good event, can be clean or dirty)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define DEBUG      0
#define VERBOSE    0

/*  ------------------------------------------------------------ */


int eventprocess(MJDetInfo *Dets, MJRunInfo *runInfo, int nChData, BdEvent *ChData[]) {

  int    i, j, k;
  int    ievt, ievt2, ievt_hg, ievt_lg, idet, chan=0;
  long long int  time;

  static int   module_lu[NCRATES+1][21];  // lookup table to map VME crate&slot into module IDs
  static int   det_lu[NBDS][16];          // lookup table to map module&chan into detector IDs
  static int   chan_lu[NBDS][16];         // lookup table to map module&chan into parameter IDs
  static int   totevts=0, builtevts=0, goodevts = 0;
  static int   *his[400];
  static double *diff[300];   // HG/LG

  static int first = 1;

  char   spname[256];
  FILE   *f_out;
  float  fsig_lg[8192], fsig_hg[8192];

  /* --------------- initialization --------------- */
  if (first) {
    if (ep_init(Dets, runInfo, module_lu, det_lu, chan_lu) < 0) return -1;
    /* malloc and clear histogram space */
    if ((his[0] = calloc(400*8192, sizeof(int))) == NULL) {
      printf("ERROR in eventprocess(); cannot malloc his!\n");
      exit(-1);
    }
    if ((diff[0] = calloc(300*8192, sizeof(double))) == NULL) {
      printf("ERROR in eventprocess(); cannot malloc diff!\n");
      exit(-1);
    }
    for (i=1; i<400; i++) his[i] = his[i-1]+8192;
    for (i=1; i<300; i++)  diff[i] = diff[i-1]+8192;

    first = 0;
  }

  /* ------------ finalization / end-of-run handling ------------- */
  if (nChData == -99) {  // special value as flag to trigger finalization

    for (idet=0; idet < 58; idet++) {
      for (j=0; j<8192; j++) {
        if (his[idet][j] > 20) {
          his[idet+100][j] = lrint(100.0 * diff[idet][j] / (double) his[idet][j]);
          his[idet+200][j] = lrint(100.0 * diff[idet+100][j] / (double) his[idet][j]);
          his[idet+300][j] = lrint(100.0 * diff[idet+200][j] / (double) his[idet][j]);
        }
      }
      // remove main slope
      if (1) {
        i = his[idet+200][4200];
        k = his[idet+200][5700] - i;
        if (i != 0 && k != 0) {
          for (j=0; j<8192; j++) {
            if (his[idet][j] > 20) {
              his[idet+100][j] -= i + k * (j-4200) / 1500;
              his[idet+200][j] -= i + k * (j-4200) / 1500;
              his[idet+300][j] -= i + k * (j-4200) / 1500;
            }
          }
        }
      }
    }
    printf("\n");

    /* save histograms */
    if ((f_out = fopen("inlcheck.rms", "w"))) {
      for (i=0; i<358; i++) {
        if (i%100 > 57) continue;
        if (i < 100) {
          sprintf(spname, "%d; LG ADC freq, run %d; ch %d",
                  i, runInfo->runNumber, i);
        } else if (i < 200) {
          sprintf(spname, "%d; HG-LG ADC*100 uncorrected, run %d; ch %d",
                  i, runInfo->runNumber, i-100);
        } else if (i < 300) {
          sprintf(spname, "%d; HG-LG ADC*100 corrected, run %d; ch %d",
                  i, runInfo->runNumber, i-200);
        } else {
          sprintf(spname, "%d; HG-LG ADC*100 HG-uncorr, run %d; ch %d",
                  i, runInfo->runNumber, i-300);
        }
        write_his(his[i], 8192, i, spname, f_out);
      }
      fclose(f_out);
    }

    return 0;
  }

  /* ------------- process input events  ------------- */

  builtevts++;
  totevts += nChData;

  // fillEvent: fill out entris in BdEvent ChData; chan, det, e, sig...
  fillEvent(runInfo, nChData, ChData, module_lu, det_lu, chan_lu);
  
  /* LOOP over channel-events in the built-event again to handle data */
  for (ievt = 1; ievt < nChData; ievt++) {
    if (ChData[ievt]->orca_type != runInfo->dataIdGM &&
        ChData[ievt]->orca_type != runInfo->dataIdGA) continue;  // only interested in Ge data
    time  = ChData[ievt]->time;
    idet  = ChData[ievt]->det;
    chan  = ChData[ievt]->chan;
    if (idet < 0 || idet > 57 || chan < 0 || chan > 199) continue;

    for (ievt2 = 0; ievt2 <ievt; ievt2++) {
      if (ChData[ievt2]->orca_type != runInfo->dataIdGM &&
          ChData[ievt2]->orca_type != runInfo->dataIdGA) continue;  // only interested in Ge data
      if (idet  != ChData[ievt2]->det ||
          (chan != ChData[ievt2]->chan+100 && chan != ChData[ievt2]->chan-100) ||
          time < ChData[ievt]->time - 6 ||
          time > ChData[ievt]->time + 6) continue;

      // have now found a high-gain/low-gain pair
      if (chan == ChData[ievt2]->chan+100) {
        ievt_hg = ievt2;
        ievt_lg = ievt;
      } else {
        ievt_hg = ievt;
        ievt_lg = ievt2;
      }

      // check that energies of HG/LG signals match
      double s1=0;
      double s2=0;
      double gain = Dets[idet].HGcalib[0]/Dets[idet].LGcalib[0];
      for (j=0; j<10; j++) {
        s1 += ChData[ievt_hg]->sig[1300+j] - ChData[ievt_hg]->sig[800+j];
        s2 += ChData[ievt_lg]->sig[1300+j] - ChData[ievt_lg]->sig[800+j];
      }
      s1 *= gain;
      if (idet != 8 && (s1 < s2-400 || s1 > s2+400))
        printf(">>>> HG-LG energy mismatch (%5.0f %5.0f); det %2d, timestamp %lld\n",
               s1/10.0, s2/10.0, idet, time);
      if (s1 < s2-200 || s1 > s2+200) continue;  // discard pairs without a good match

      // do INL correction
      if (inl_correct(Dets, runInfo, ChData[ievt_hg]->sig, fsig_hg, 2008, ChData[ievt_hg]->chan) ||
          inl_correct(Dets, runInfo, ChData[ievt_lg]->sig, fsig_lg, 2008, ChData[ievt_lg]->chan)) {
        printf(" >>> inl_correct return error for det %d!\n", idet);
        continue;
      }

      // store hg-lg differences
      for (j=1300; j<2000; j++) {
        k = ChData[ievt_lg]->sig[j] + 4096;
        if (k < 0 || k > 8000) continue;
        his[idet][k]++;
        diff[idet][k] += 4096.0  + (double) ChData[ievt_hg]->sig[j]*gain - (double) k;
        diff[idet+100][k] += fsig_hg[j]*gain - fsig_lg[j];
        diff[idet+200][k] += (double) ChData[ievt_hg]->sig[j]*gain - fsig_lg[j];
        //diff[idet+100][k] += ((double) (fsig_hg[j] - ChData[ievt_hg]->sig[j]));

      }
    }
  /* ------------- end of Ge detector data processing -------------- */

  } // END of loop over channel-events

  goodevts++;
  if (totevts % 20000 < (totevts-nChData) % 20000) {
    printf(" Processed %d builtevts in, %d good evts out\n",
           builtevts, goodevts);
  }
  return  0;
} /* eventprocess() */

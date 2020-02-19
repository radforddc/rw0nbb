/*
   ep_finalize.c

   event processing code for MJD
   - up to 16 digitizer modules for the Ge detectors, plus
   - QDCs etc for the veto

   David Radford   Nov 2016
*/
/* Close out processing of the events that have been handled by eventprocess()
   - write histograms to file
   - if there are enough data, also:
      - extract information on on-board trap offsets and thresholds, and tmax
      - extract information on baseline (E=0) trapfixed position and FWHM
      - extract data-cleaning limits (on resting baseline mean, RMS, and slope)
                and save the limits in a *.dcl file
      - extract live times from pulser counts

   returns: -1 on error, 0 otherwise
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"
#define HIS_COUNT    2200     // number of spectra in his[][] array

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */

#define VERBOSE    0
/*  ------------------------------------------------------------ */

float autopeak3(int *his, int lo, int hi, float *area_ret, float *fwhm_ret);

void find_limits(int *his, int lo, int hi, int *lolim, int *hilim) {
  // find limiting non-zero channels between chs lo and hi
  int i;
  for (i=lo; i<=hi && his[i]==0; i++);
  if (i > hi) return;  // all zeros
  *lolim = *hilim = i;
  for (i++; i<=hi; i++) if (his[i] > 0) *hilim = i;
  return;
}

/*  ------------------------------------------------------------ */

int ep_finalize(MJDetInfo *Dets, MJRunInfo *runInfo, int **his,
                int *on_bd_rise, int totevts, PTag  *pt, DataClean *dcInfo) {

  int     i, j, k, n, t, chan, e, s1, s2, flag;
  float   pos, area, fwhm, area2;
  float   ol, oh, tl, th, fl, th2, tl2, fh, rl, rh;
  char    spname[256], *c;
  FILE    *f_out;
  PTag    pt2;

  static char  cg[2][4] = {"HG", "LG"};
  static float e0pos[200], e0fwhm[200], bl[200], t1[200], obt_offset[200], blrms[200];
  static int   bl_lo[200], bl_hi[200], blrms_hi[200], blsl_lo[200], blsl_hi[200];

  /*
    detector channels:
         0 -  57:  Ge high gain
       100 - 157:  Ge low gain
       158 - 168:  pulser tag
    his channels:
         0 - 199   pulser-removed trapmax (4-2-4)
       200 - 399   pulser-tagged trapmax (4-2-4)
       400 - 599   dirty energy (e_ctc or trapmax) after data cleaning; ADC units
       600 - 799   clean energy (e_ctc or trapmax) after data cleaning; ADC units
       800 - 999   dirty energy (e_ctc or trapmax) after data cleaning; 0.5 keV/ch
      1000 - 1199  clean energy (e_ctc or trapmax) after data cleaning; 0.5 keV/ch
      1200 - 1399  A/E histograms
      1400 - 1599  DCR and lamda histograms
      1600 - 1799  baseline mean, RMS, and slope (+4000 when trapmax is < 50)
      1800 - 1999  tmax from trapmax (4-2-4) and 2000+t_1 (1% time)
                     and 6000+baseline_trap (E=0, 4-1.5-4) for FWHM, 0.05 ADC units
      2000 - 2199  trap difference (on-board - trapmax) for trapmax > 100
                     in bins of 0.01 (+1000) and 0.2 (+3000) ADC
                     and (t95 - t0) time difference for different A/E and DCR cuts 
   */


  f_out = fopen("his.rms", "w");
  for (i=0; i<HIS_COUNT && his[i] != 0; i++) {
    if (i == 95) {
      sprintf(spname, "%d; no-pulser built-event times [s]", i);
    } else if (i == 96) {
      sprintf(spname, "%d; pulser built-event times [s]", i);
    } else if (i == 97) {
      sprintf(spname, "%d; no-pulser ch-event times [s]", i);
    } else if (i == 98) {
      sprintf(spname, "%d; pulser ch-event times [s]", i);
    } else if (i == 99) {
      sprintf(spname, "%d; channel and crate ID stats", i);
    } else if (i < 200) {
      sprintf(spname, "%d; ch %d pulser-removed energy [ADC]; run %d", i, i, runInfo->runNumber);
    } else if (i < 400) {
      sprintf(spname, "%d; ch %d pulser-tagged energy [ADC]", i, i%200);
    } else if (i < 600) {
      sprintf(spname, "%d; ch %d dirty energy (e_ctc or trapmax) [ADC]", i, i%200);
    } else if (i < 800) {
      sprintf(spname, "%d; ch %d clean energy (e_ctc or trapmax) [ADC]", i, i%200);
    } else if (i < 1000) {
      sprintf(spname, "%d; ch %d dirty energy (e_ctc or trapmax) [0.5 keV]", i, i%200);
    } else if (i < 1200) {
      sprintf(spname, "%d; ch %d clean energy (e_ctc or trapmax) [0.5 keV]e", i, i%200);
    } else if (i < 1400) {
      sprintf(spname, "%d; ch %d A/E", i, i%200);
    } else if (i < 1600) {
      sprintf(spname, "%d; ch %d DCR and lamda", i, i%200);
    } else if (i < 1800) {
      sprintf(spname, "%d; ch %d baseline mean, RMS, slope (+4000 for e<50)", i, i%200);
    } else if (i < 2000) {
      sprintf(spname, "%d; ch %d trap_t_max and t_1 [10 ns] and E=0 [0.05 ADC]", i, i%200);
    } else {
      sprintf(spname, "%d; ch %d on-bd trap offset [0.01 & 0.2 ADC]", i, i%200);
    }
    write_his(his[i], 8192, i, spname, f_out);
  }
  fclose(f_out);

  /* presorted data subsets have few events, not enough to ectract trap offsets,
     thresholds, data cleaning limits, etc... */
  if (totevts < 2000) return 0; // skip remaining end-of-run processing

  /*  ======  extract on-board trap offsets and thresholds, and tmax  ======  */
  if (Dets[0].HGcalib[0] > 0) {
    printf("\n"
           "                               On-bd                    True               T_1"
           "               On-bd                     True              T_1\n"
           "    Detector    rise  HGChan   offset  (keV)         threshold  (keV)     time"
           "     LGChan    offset  (keV)         threshold  (keV)     time\n");
  } else {
    printf("\n"
           "                               On-bd            True        T_1"
           "              On-bd            True        T_1\n"
           "    Detector    rise  HGChan   offset       threshold      time"
           "     LGChan   offset       threshold      time\n");
  }
  for (chan=0; chan < 100+runInfo->nGe; chan++) {
    obt_offset[chan] = 0;
    if (chan >= runInfo->nGe && chan < 100) continue;
    /* look for peak around bin 1000 (0.01 ADC) */
    if ((e = peak_find(his[2000+chan], 1, 1999))) { // peak found
      if ((pos = autopeak(his[2000+chan], e, &area, &fwhm)) == 0) pos = e;
      pos = 0.01 * (pos - 1000.0);
    } else {  /* look for peak above bin 2000 (0.2 ADC)*/
      if (!(e = peak_find(his[2000+chan], 2001, 3999)))  continue; // no peak found
      if ((pos = autopeak(his[2000+chan], e, &area, &fwhm)) == 0) pos = e;
      pos = 0.2 * (pos - 3000.0);
    }
    obt_offset[chan] = pos;
    // tmax = peak_find(his[1800+chan], 1, 1000); // find position of trap tmax peak
    if ((t = peak_find(his[1800+chan], 2001, 4000))) // find position of t_1 peak
      t1[chan] = autopeak(his[1800+chan], t, &area, &fwhm) - 2000.0;
  }
  for (chan=0; chan < runInfo->nGe; chan++) {
    if (!Dets[chan].HGChEnabled || (t1[100+chan] == 0 && t1[chan] == 0)) continue;
    if (Dets[0].HGcalib[0] > 0) {
      printf("%8s %s %4d %6d", Dets[chan].DetName, Dets[chan].StrName, on_bd_rise[chan], chan);
      if (t1[chan] != 0) {
        printf("%9.2f (%5.2f) %7.2f ->%7.2f (%4.2f) %8.1f",
               obt_offset[chan], obt_offset[chan] * Dets[chan].HGcalib[0],
               Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan],
               Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan],
               (Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[chan]) * Dets[chan].HGcalib[0],
               t1[chan]);
      } else {
        printf("%51s", " ");
      }
      if (t1[100+chan] != 0) {
        printf("%10d %9.2f (%5.2f)  %7.2f ->%7.2f (%4.2f) %8.1f\n",
               chan+100, obt_offset[100+chan], obt_offset[100+chan] * Dets[chan].LGcalib[0],
               Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan],
               Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[100+chan],
               (Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[100+chan]) * Dets[chan].LGcalib[0],
               t1[100+chan]);
      } else {
        printf("\n");
      }
    } else { 
      printf("%8s %s %4d %6d", Dets[chan].DetName, Dets[chan%100].StrName, on_bd_rise[chan], chan);
      if (t1[chan] != 0) {
        printf("%9.2f %7.2f ->%7.2f %8.1f",
               obt_offset[chan],
               Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan],
               Dets[chan].HGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[100+chan],
               t1[chan]);
      } else {
        printf("%36s", " ");
      }
      if (t1[100+chan] != 0) {
        printf("%10d %9.2f %7.2f ->%7.2f %8.1f\n",
               chan+100, obt_offset[100+chan],
               Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan],
               Dets[chan].LGTrapThreshold/ (float) on_bd_rise[chan] - obt_offset[100+chan],
               t1[100+chan]);
      } else {
        printf("\n");
      }
    }
  }

  /*  ======  baseline (E=0) trapfixed position and FWHM ======  */
  for (i=0; i<200; i++) {
    e0pos[i] = e0fwhm[i] = bl[i] = t1[i] = 0;
  }

  printf("\n   NOTE that correct E=0 position and FWHM require 9us of flat signal baseline.\n");
  if (Dets[0].HGcalib[0] > 0) {
    printf("    Detector     HGChan   E=0  (keV)   FWHM   (keV)  area"
           "     LGChan   E=0  (keV)   FWHM   (keV)  area    diff\n");
  } else {
    printf("    Detector     HGChan   E=0  FWHM  area     LGChan   E=0  FWHM  area    diff\n");
  }
  for (chan=0; chan < runInfo->nGe; chan++) {
    if (!Dets[chan].HGChEnabled) continue;
    // high gain
    if ((pos = autopeak2(his[1800+chan], 4000, 8000, &area, &fwhm)) &&
        area > 100) {
      e0pos[chan] = 0.05 * (pos-6000);
      e0fwhm[chan] = 0.05 * fwhm;
    }
    // low gain
    if ((pos = autopeak2(his[1900+chan], 4000, 8000, &area2, &fwhm)) &&
        area2 > 100) {
      e0pos[100+chan] = 0.05 * (pos-6000);
      e0fwhm[100+chan] = 0.05 * fwhm;
    }
    if (area <= 100 && area2 <= 100) continue;
    printf("%8s %s %6d", Dets[chan].DetName, Dets[chan%100].StrName, chan);

    if (area > 100) {
      if (Dets[0].HGcalib[0] > 0) {
        printf("%7.2f (%5.2f) %5.2f (%5.3f) %5.0f", e0pos[chan], e0pos[chan] * Dets[chan].HGcalib[0],
               e0fwhm[chan], e0fwhm[chan] * Dets[chan].HGcalib[0], area);
      } else {
        printf("%7.2f %5.2f %5.0f", e0pos[chan], e0fwhm[chan], area);
      }
    } else {
      printf("%21s", " ");
      if (Dets[0].HGcalib[0] > 0) printf("%14s", " ");
    }
    // low gain
    if (area2 > 100) {
      if (Dets[0].HGcalib[0] > 0) {
        printf("%10d %6.2f (%5.2f) %5.2f (%5.3f) %5.0f %7.0f\n",
               chan+100, e0pos[100+chan], e0pos[100+chan] * Dets[chan].LGcalib[0],
               e0fwhm[100+chan], e0fwhm[100+chan] * Dets[chan].LGcalib[0], area2, area-area2);
      } else {
        printf("%10d %6.2f %5.2f %5.0f %7.0f\n",
               chan+100, e0pos[100+chan], e0fwhm[100+chan], area2, area-area2);
      }
    } else {
      printf("\n");
    }
  }
  printf("\n");

  /*  ======  data-cleaning limits (resting baseline mean, RMS, slope) ======  */
  /* initialize data-cleaning values to defaults */
  for (i=0; i<100; i++) {
    bl_lo[i]    = -99;   //
    bl_hi[i]    = 699;   // these seem like conservative (loose) limits
    blrms_hi[i] =  49;   //      for HG chs with no pulser signals
    blsl_lo[i]  = -99;   //
    blsl_hi[i]  =  99;   //
    //dcInfo->modified[i] = 0;
  }
  for (i=100; i<200; i++) {
    bl_lo[i]    = -299;  //
    bl_hi[i]    = 599;   // these seem like conservative (loose) limits
    blrms_hi[i] =  49;   //      for LG chs with no pulser signals
    blsl_lo[i]  = -99;   //
    blsl_hi[i]  =  99;   //
    //dcInfo->modified[i] = 0;
  }

  if (VERBOSE) {
    printf("\n  Data-cleaning limits from baselines:\n");
    printf("  Chan     Detector          baseline value      RMS      slope\n");
  }
  for (chan=0; chan < 100+runInfo->nGe; chan++) {
    if (chan >= runInfo->nGe && chan < 100) continue;
    if (chan == 100 && VERBOSE) 
      printf("  Chan     Detector          baseline value      RMS      slope\n");
    if ((pos = autopeak2(his[1600+chan], 0, 1999, &area, &fwhm)) && area > 100) {
      pos = autopeak2(his[1600+chan], pos-3.0*fwhm,  pos+3.0*fwhm, &area, &fwhm);
      dcInfo->modified[chan] = 1;
      // printf("dcl for chan %d modified...\n", chan);
      bl[chan] = 0.3 * (pos-1000.0);                  // mean baseline is in 0.3 ADC units
      if (chan > 99) bl[chan] = 0.1 * (pos-1000.0);   // mean baseline is in 0.1 ADC units
      find_limits(his[1600+chan], pos-2.0*fwhm, pos+2.0*fwhm, &bl_lo[chan], &bl_hi[chan]);
      bl_lo[chan] -= 1000;
      bl_hi[chan] -= 1000;
      if (chan > 99) find_limits(his[1600+chan], 2000, 2500, &j, &blrms_hi[chan]);

      blrms[chan] = s1 = s2 = 0;
      for (i=2000; i < 2100; i++) {
        s1 += (i-2000) * his[1600+chan][i];
        s2 += his[1600+chan][i];
      }
      if (s2 > 50) {
        blrms[chan] = 0.3 * (float) s1 / (float) s2;                // baseline RMS, 0.3 ADC units
        if (chan > 99) blrms[chan] = 0.1 * (float) s1 / (float) s2; // baseline RMS, 0.1 ADC units
        blrms_hi[chan] = 2.0 * (float) s1 / (float) s2;     // baseline RMS limit = twice mean value
      }

      find_limits(his[1600+chan], 2500, 3500, &blsl_lo[chan], &blsl_hi[chan]);
      if ((pos = autopeak2(his[1600+chan], 2500, 3500, &area, &fwhm)) && area > 100)
        find_limits(his[1600+chan], pos-1.7*fwhm, pos+1.7*fwhm, &blsl_lo[chan], &blsl_hi[chan]);
      blsl_lo[chan] -= 3000;
      blsl_hi[chan] -= 3000;
      if (VERBOSE)
        printf("%4d  %8s %s %s   %9.2f %4d %4d %7d %7d %4d\n",
               chan, Dets[chan%100].DetName, Dets[chan%100].StrName, cg[chan/100],
               bl[chan], bl_lo[chan], bl_hi[chan], blrms_hi[chan],
               blsl_lo[chan], blsl_hi[chan]);
    }
  }

  /*  ============  save data-cleaning-limit results in a file  ============  */
  // if (runInfo->analysisPass >= 0) {
  j = 0;
  for (i=0; i<200; i++) {
    if (dcInfo->modified[i]) {
      //printf("dcl for chan %d modified...\n", i);
      j = 1;
      dcInfo->bl[i]       = bl[i];
      dcInfo->bl_lo[i]    = bl_lo[i];
      dcInfo->bl_hi[i]    = bl_hi[i];
      dcInfo->blrms_hi[i] = blrms_hi[i];
      dcInfo->blsl_lo[i]  = blsl_lo[i];
      dcInfo->blsl_hi[i]  = blsl_hi[i];
    }
  }
  if (j) data_clean_info_write(runInfo, dcInfo);
  // }
  /* ====== extract live times from pulser counts ====== */

  /* ptag.nevts[chan][0]: energy-ungated HG pulser cts
     ptag.nevts[chan][1]: energy-ungated LG pulser cts
     ptag.nevts[chan][2]: energy-ungated HG&&LG pulser cts
     ptag.nevts[chan][3]: energy-gated HG pulser cts
     ptag.nevts[chan][4]: energy-gated LG pulser cts
     ptag.nevts[chan][5]: energy-gated HG&&LG pulser cts
     ptag.nevts[chan][6]: expected energy-ungated pulser cts (from finding max over the CC)
   */
  for (i=0; i < runInfo->nCC; i++) {
    // find max # pulser counts between all channels on a given controller card
    k = 0;
    for (j=0; j < runInfo->nGe; j++) {
      n = pt->nevts[j][0] + pt->nevts[j][1] - pt->nevts[j][2]; // count for either HG or LG
      if (j != 16 && Dets[j].CCnum == i && k < n) k = n;
    }
    //  and save in pt->nevts[][6] (expected # pulser counts)
    for (j=0; j < runInfo->nGe; j++)
      if (Dets[j].CCnum == i) pt->nevts[j][6] = k;
  }
  
#ifndef QUIET
  /* report results */
  printf("\n"
         "                           Pulser Counts                         LiveTime\n"
         "        Detector        HG    LG  Both  Expected     diff      Either   (Both)\n");
  for (chan=0; chan < runInfo->nGe; chan++) {
    if (!Dets[chan].HGChEnabled ||   // ch not enabled
        pt->nevts[chan][6] < 2 ||   // no pulser expected
        (pt->nevts[chan][0] == 0 && pt->nevts[chan][1] == 0)) continue;  // no pulser
    n = pt->nevts[chan][0] + pt->nevts[chan][1] - pt->nevts[chan][2]; // either
    printf("%3d %8s %s %6d %5d %5d %7d %6d %4d %10.4f  (%6.4f)\n",
           chan, Dets[chan].DetName, Dets[chan%100].StrName,
           pt->nevts[chan][0], pt->nevts[chan][1],
           pt->nevts[chan][2], pt->nevts[chan][6],
           pt->nevts[chan][6] - pt->nevts[chan][0],
           pt->nevts[chan][6] - pt->nevts[chan][1],
           (float) n / (float) pt->nevts[chan][6],
           (float) pt->nevts[chan][2] / (float) pt->nevts[chan][6]);
  }
#endif
  if (0) {     // change to 0 to 1 for summation of pulser counts over many data subsets, and...
    if ((f_out = fopen("ptagtot.dat", "r"))) {  // ... delete ptagtot.dat to reset sums to zero
      fread(&pt2, sizeof(pt2), 1, f_out);
      fclose(f_out);
      for (chan=0; chan < runInfo->nGe; chan++) {
        for (i=0; i<8; i++) pt2.nevts[chan][i] += pt->nevts[chan][i];
      }
      f_out = fopen("ptagtot.dat", "w");
      fwrite(&pt2, sizeof(pt2), 1, f_out);
      fclose(f_out);
      printf("\n"
             "                           Pulser Counts                         LiveTime\n"
             "        Detector        HG    LG  Both  Expected     diff      Either   (Both)\n");
      for (chan=0; chan < runInfo->nGe; chan++) {
        if (pt2.nevts[chan][6] < 2 ||   // no pulser expected
            (pt2.nevts[chan][0] == 0 && pt2.nevts[chan][1] == 0)) continue;  // no pulser
        n = pt2.nevts[chan][0] + pt2.nevts[chan][1] - pt2.nevts[chan][2]; // either
        printf("%3d %8s %s %6d %5d %5d %7d %6d %4d %10.4f  (%6.4f)\n",
               chan, Dets[chan].DetName, Dets[chan%100].StrName,
               pt2.nevts[chan][0], pt2.nevts[chan][1],
               pt2.nevts[chan][2], pt2.nevts[chan][6],
               pt2.nevts[chan][6] - pt2.nevts[chan][0],
               pt2.nevts[chan][6] - pt2.nevts[chan][1],
               (float) n / (float) pt2.nevts[chan][6],
               (float) pt2.nevts[chan][2] / (float) pt2.nevts[chan][6]);
      }
    } else {
      f_out = fopen("ptagtot.dat", "w");
      fwrite(pt, sizeof(pt2), 1, f_out);
      fclose(f_out);
    }
  }

  for (i=0; i < runInfo->nGe; i++) {
    if (Dets[i].HGcalib[0] < 0.1) Dets[i].HGcalib[0] = 1.0;
    if (Dets[i].LGcalib[0] < 0.1) Dets[i].LGcalib[0] = 1.0;
  }
  

  /*  ============  save some results in a file, for plotting   ============ */
  if (!(f_out = fopen("run_data_for_plots.dat", "w"))) return 0;
  for (i=0; i < runInfo->nGe; i++) {
    if (Dets[i].HGcalib[0] < 0.1) Dets[i].HGcalib[0] = 1.0;
    if (Dets[i].LGcalib[0] < 0.1) Dets[i].LGcalib[0] = 1.0;
  }

  c = runInfo->filename;
  while (strchr(c, '/')) c = strchr(c, '/')+1; // remove initial path from input file name 
  fprintf(f_out,
          "# Input file: %s\n"
          "#HGch  OBT-offset    Threshold     LGch    OBT-offset    Threshold     rise  detector\n"
          "#      ADC    keV    ADC    keV            ADC    keV    ADC    keV\n", c);
  for (i=0; i < runInfo->nGe; i++) {
    if (!Dets[i].HGChEnabled) continue;
    oh = obt_offset[i] * Dets[i].HGcalib[0];
    ol = obt_offset[i+100] * Dets[i].LGcalib[0];
    th = (Dets[i].HGTrapThreshold / (float) on_bd_rise[i] - obt_offset[i]) * Dets[i].HGcalib[0];
    tl = (Dets[i].LGTrapThreshold / (float) on_bd_rise[i] - obt_offset[i+100]) * Dets[i].LGcalib[0];
    if (oh < -2.0) oh = -2.0;
    if (oh >  2.0) oh =  2.0;
    if (ol < -2.0) ol = -2.0;
    if (ol >  2.0) ol =  2.0;
    if (th >  5.0) th =  5.0;
    if (tl >  5.0) tl =  5.0;
    fprintf(f_out, "%3d %6.2f %6.2f %6.2f %6.2f %7d %6.2f %6.2f %6.2f %6.2f %7d   %s\n",
            i, obt_offset[i], oh, //obt_offset[i] * Dets[i].HGcalib[0],
            Dets[i].HGTrapThreshold/ (float) on_bd_rise[i] - obt_offset[i], th,
            //(Dets[i].HGTrapThreshold/ (float) on_bd_rise[i] - obt_offset[i]) * Dets[i].HGcalib[0],
            i+100, obt_offset[i+100], ol, //obt_offset[i+100] * Dets[i].LGcalib[0],
            Dets[i].LGTrapThreshold/ (float) on_bd_rise[i] - obt_offset[i+100], tl,
            //(Dets[i].LGTrapThreshold/ (float) on_bd_rise[i] - obt_offset[i+100]) * Dets[i].LGcalib[0],
            on_bd_rise[i], Dets[i].StrName);
  }

  fprintf(f_out, "\n"
          "#HGch   E=0 energy         FWHM            RMS      LGch     E=0 energy         FWHM           RMS     detector\n"
          "#       ADC    keV      ADC   keV      ADC   keV             ADC    keV      ADC   keV      ADC   keV\n");
  for (i=0; i < runInfo->nGe; i++) {
    if (!Dets[i].HGChEnabled || e0fwhm[i] < 0.01) continue;
    fh = e0fwhm[i] * Dets[i].HGcalib[0];
    fl = e0fwhm[i+100] * Dets[i].LGcalib[0];
    rh = blrms[i] * Dets[i].HGcalib[0];
    rl = blrms[i+100] * Dets[i].LGcalib[0];
    if (fh > 1.0) fh = 1.0;
    if (fl > 1.0) fl = 1.0;
    if (rh > 3.0) rh = 2.0;
    if (rl > 3.0) rl = 2.0;
    fprintf(f_out, "%3d %7.2f %7.3f %7.2f %6.3f %7.2f %5.2f %7d %7.2f %7.3f %7.2f %6.3f %7.2f %5.2f   %s\n",
            i, e0pos[i], e0pos[i] * Dets[i].HGcalib[0],
            e0fwhm[i], fh, //e0fwhm[i] * Dets[i].HGcalib[0],
            blrms[i], blrms[i] * Dets[i].HGcalib[0],
            i+100, e0pos[i+100], e0pos[i+100] * Dets[i].LGcalib[0],
            e0fwhm[i+100], fl, //e0fwhm[i+100] * Dets[i].LGcalib[0],
            rh, rl, //blrms[i+100], blrms[i+100] * Dets[i].LGcalib[0],
            Dets[i].StrName);
  }
  fclose(f_out);

  /* ===== make a file listing deadtime vs. (neg_threshold/sigma) ===== */
  if (0) {  // change 0 to 1 to enable creation of DT_vs_thresh.dat file
    if (!(f_out = fopen("DT_vs_thresh.dat", "w"))) return 0;
    fprintf(f_out,
            "#                 HG                     LG\n"
            "# detId    dead%%  thresh/sgima   dead%%  thresh/sgima\n");
    for (i=0; i < runInfo->nGe; i++) {
      if (Dets[i].HGChEnabled && Dets[i].LGChEnabled &&
          pt->nevts[i][6] > 99 &&
          e0fwhm[i]> 0 && e0fwhm[i+100] > 0) {
        th = Dets[i].HGTrapThreshold / (float) on_bd_rise[i] - obt_offset[i];
        tl = Dets[i].LGTrapThreshold / (float) on_bd_rise[i] - obt_offset[i+100];
        th2 = Dets[i].HGTrapThreshold / (float) on_bd_rise[i] + obt_offset[i];
        tl2 = Dets[i].LGTrapThreshold / (float) on_bd_rise[i] + obt_offset[i+100];
        fh = e0fwhm[i]/2.355;
        fl = e0fwhm[i+100]/2.355;;
        rh = 100.0 * (float) (pt->nevts[i][6] - pt->nevts[i][0]) / (float) pt->nevts[i][6];
        rl = 100.0 * (float) (pt->nevts[i][6] - pt->nevts[i][1]) / (float) pt->nevts[i][6];
        fprintf(f_out, "%2d %6.1f %5.2f %6.3f %6.1f %5.2f %6.3f    %s\n",
                i, rh, th2/fh, th2/th, rl, tl2/fl, tl2/tl, Dets[i].StrName);
      }
    }
    fclose(f_out);
  }

  /* ===== make a file listing channels for which the threshold finder should be re-run ===== */
  if (0) {  // change 0 to 1 to enable creation of TFChannles.txt file
    if (!(f_out = fopen("TFChannels.txt", "w"))) return 0;
    for (i=0; i < runInfo->nGe; i++) {
      if (Dets[i].HGChEnabled &&
          obt_offset[i] < -0.3 * Dets[i].HGTrapThreshold / (float) on_bd_rise[i])
        fprintf(f_out, "%d,%d,%d    %s HG\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanHi, Dets[i].StrName);
      if (Dets[i].LGChEnabled &&
          obt_offset[i+100] < -0.3 * Dets[i].LGTrapThreshold / (float) on_bd_rise[i])
        fprintf(f_out, "%d,%d,%d    %s LG\n",
                Dets[i].crate, Dets[i].slot, Dets[i].chanLo, Dets[i].StrName);
    }
    fclose(f_out);
  }

  /* finally check for cases where, for the pulser, we count both > either channel by itself;
     or the energy gate misses a large fraction of the events.
     This can happen when there's lots of noise */
  /* ptag.nevts[chan][0]: energy-ungated HG pulser cts
     ptag.nevts[chan][1]: energy-ungated LG pulser cts
     ptag.nevts[chan][2]: energy-ungated HG&&LG pulser cts
     ptag.nevts[chan][3]: energy-gated HG pulser cts
     ptag.nevts[chan][4]: energy-gated LG pulser cts
     ptag.nevts[chan][5]: energy-gated HG&&LG pulser cts
     ptag.nevts[chan][6]: expected energy-ungated pulser cts (from finding max over the CC)
   */
  flag = 0;
  for (j=0; j < runInfo->nGe; j++) {
    if (pt->nevts[j][2] > pt->nevts[j][0] ||
        pt->nevts[j][2] > pt->nevts[j][1])
      printf("Error in pulser counts for detector %d; Both > HG or Both > LG\n", j);
    if ((j != 16 || !SPECIAL_DET_16) &&
        (pt->nevts[j][3] < pt->nevts[j][0]*4/5 ||
         pt->nevts[j][4] < pt->nevts[j][1]*4/5)) flag = j;
  }
 
  if (flag) {
    printf("\n >>>>>>  WARNING: The pulser tag is having problems with the energy gates.  <<<<<<\n"
           "   The energies may have changed, or the resolution may be a lot worse for\n"
           "       some detectors than what is given in the .pdt file.\n"
           "   Consider re-running pulser_tag_init, or discarding some detectors for this run.\n"
           "   Problem detectors are: ");
    for (j=0; j < runInfo->nGe; j++) {
      if (j != 16 &&
          (pt->nevts[j][3] < pt->nevts[j][0]*4/5 ||
           pt->nevts[j][4] < pt->nevts[j][1]*4/5))
        printf("%3d", j);
    }
    printf("\n");
  }

  /* find position of 2614.5-keV peak and use that to compute the calibration gains */
  for (chan=0; chan<100+runInfo->nGe; chan++) {
    if (chan >= runInfo->nGe && chan < 100) continue;
    if (chan%100 == 0) printf("\n");
    fwhm = 20;
    j = 5000;
    if (chan > 99) {
      fwhm = 8;
      j = 1700;
    }
    if ((pos = autopeak3(his[600+chan], j, 7000, &area, &fwhm))) {
      if (chan < 100) {
        Dets[chan].HGcalib[0] = 2614.5/pos;
        printf(" 2615: %3d P = %7.1f A = %7.0f; FWHM = %7.1f keV\n",
               chan, pos, area, fwhm*Dets[chan].HGcalib[0]);
      } else {
        Dets[chan-100].LGcalib[0] = 2614.5/pos;
        printf(" 2615: %3d P = %7.1f A = %7.0f; FWHM = %7.1f keV\n",
               chan, pos, area, fwhm*Dets[chan-100].LGcalib[0]);
      }
    }
  }
  if ((f_out = fopen("gains.output","w"))) {
    printf("\n Writing energy calibrations to gains.output\n");
    for (chan=0; chan<runInfo->nGe; chan++)
      fprintf(f_out,
              "%3d %10.8lf %10.8lf %s\n",
              chan, Dets[chan].HGcalib[0], Dets[chan].LGcalib[0], Dets[chan].StrName);
  }

  return 0;  // END of post-run processing

} /* ep_finalize() */

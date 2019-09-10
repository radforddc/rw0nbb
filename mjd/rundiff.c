/*
 *  rundiff.c -- read and decode the XML file headers from two MJD raw data files
 *               and look for differences in the important parameters
 *
 *    rundiff <file1> <file2> 
 *    returns: -1 on file open/read error, or number of differences found
 *              i.e. 0 for no significant difference
 *
 *  David Radford   Dec 2016
 *
 *  TODO: add veto voltages
 *        add detector names to output when detector values differ
 *        add verbose flag?
 *        add noise analysis of data
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"
#include "runBits.h"

#define VERBOSE 0

/* some convenience functions */
int diff_int(int value1, int value2, char *label, char *description) {
  if (value1 == value2) return 0;
  printf("-------- Number of %s is different!\n"
         "< %s %d\n"
         "> %s %d\n", description, label, value1, label, value2);
  return 1;
}

int diff_int_2(int value1, int value2, int i, char *label, char *description) {
  if (value1 == value2) return 0;
  printf("-------- %s number %d is different!\n"
         "< %s[%d] %d\n"
         "> %s[%d] %d\n", description, i, label, i, value1, label, i, value2);
  return 1;
}

int diff_int_det(MJDetInfo  *detInfo, int value1, int value2,
                 int i, char *label, char *description) {
  if (value1 == value2) return 0;
  printf("-------- %s for detector number %d (%s/%s) is different!\n"
         "< %s[%d] %d\n"
         "> %s[%d] %d\n", description, i,
         detInfo[i].DetName, detInfo[i].StrName,
         label, i, value1, label, i, value2);
  return 1;
}

int diff_char(char *value1, char *value2, int i, char *label, char *description) {
  if (!strncmp(value1, value2, 16)) return 0;
  printf("-------- %s number %d is different!\n"
         "< %s[%d] %s\n"
         "> %s[%d] %s\n", description, i, label, i, value1, label, i, value2);
  return 1;
}

/* ------------------------ main code ------------------------- */
/*
 *  rundiff.c -- read and decode the XML file headers from two MJD raw data files
 *               and look for differences in the important parameters
 *
 *    rundiff <file1> <file2> 
 *    returns: -1 on file open/read error, or number of differences found
 */

int main(int argc, char **argv) {

  FILE       *f_in;
  int        nDets1, nDets2;
  MJDetInfo  detInfo1[NMJDETS], detInfo2[NMJDETS];
  MJRunInfo  runInfo1, runInfo2;
  int        i, argn=1, ndiff=0;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [-v] runfile1 runfile2\n\n", argv[0]);
    return -1;
  }
  /* look for and skip over command line flags (ignored for now) */
  while (argn < argc-1 && argv[argn][0] == '-') argn ++;

  /* read in header from first file */
  if ((f_in = fopen(argv[argn],"r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", argv[argn]);
    return -1;
  }
  // printf("\n >>> Reading %s\n", argv[argn]);
  strncpy(runInfo1.filename, argv[argn], 256);
  runInfo1.argc = argc;
  runInfo1.argv = argv;
  nDets1 = decode_runfile_header(f_in, detInfo1, &runInfo1);
  fclose(f_in);
  if (nDets1 < 1) return -1;

  /* read in header from second file */
  argn++;
  if ((f_in = fopen(argv[argn],"r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", argv[argn]);
    return -1;
  }
  // printf(" >>> Reading %s\n", argv[argn]);
  strncpy(runInfo2.filename, argv[argn], 256);
  runInfo2.argc = argc;
  runInfo2.argv = argv;
  nDets2 = decode_runfile_header(f_in, detInfo2, &runInfo2);
  fclose(f_in);
  if (nDets2 < 1) return -1;

  printf("\n< file 1 from %s : %s\n"
         "> file 2 from %s : %s\n\n",
         runInfo1.date, runInfo1.filename,
         runInfo2.date, runInfo2.filename);

  /* compare results */

  /* compare general run info */
  // ndiff += diff_int(runInfo1.runType, runInfo2.runType, "runType", "run bits");
  if (runInfo1.runType != runInfo2.runType) {
    printf("-------- Run bits are different!\n"
           "< run bits 0x%8.8x\n", runInfo1.runType);
    for (i=0; i<32; i++) {
      if ((runInfo1.runType & 1<<i) && !(runInfo2.runType & 1<<i))
        printf("<          0x%8.8x - %s\n", 1<<i, runBitDesc[i]);
    }
    printf("> run bits 0x%8.8x\n", runInfo2.runType);
    for (i=0; i<32; i++) {
      if (!(runInfo1.runType & 1<<i) && (runInfo2.runType & 1<<i))
        printf(">          0x%8.8x - %s\n", 1<<i, runBitDesc[i]);
    }
    ndiff++;
  }

  ndiff += diff_int(nDets1, nDets2, "nDets", "detectors");
  ndiff += diff_int(runInfo1.nGe, runInfo2.nGe, "nGe", "Ge detectors");
  ndiff += diff_int(runInfo1.nPT, runInfo2.nPT, "nPT", "pulser tag chs");
  ndiff += diff_int(runInfo1.nV, runInfo2.nV, "nV", "veto VME modules");
  ndiff += diff_int(runInfo1.nCC, runInfo2.nCC, "nCC", "preamp controller cards");
  ndiff += diff_int(runInfo1.idNum, runInfo2.idNum, "idNum", "data-type IDs");

  /* compare pulser tag DAQ */
  if (runInfo1.nPT == runInfo2.nPT) {
    for (i=0; i<runInfo1.nPT; i++) {
      ndiff += diff_int_2(runInfo1.PTcrate[i], runInfo2.PTcrate[i], i,
                          "PTcrate", "Pulser tag crate");
      ndiff += diff_int_2(runInfo1.PTslot[i], runInfo2.PTslot[i], i,
                          "PTslot", "Pulser tag slot");
      ndiff += diff_int_2(runInfo1.PTchan[i], runInfo2.PTchan[i], i,
                          "PTchan", "Pulser tag channel");
    }
  }
  /* compare veto DAQ */
  if (runInfo1.nV == runInfo2.nV) {
    for (i=0; i<runInfo1.nV; i++) {
      ndiff += diff_int_2(runInfo1.Vcrate[i], runInfo2.Vcrate[i], i,
                          "Vcrate", "Veto module crate");
      ndiff += diff_int_2(runInfo1.Vslot[i], runInfo2.Vslot[i], i,
                          "Vslot", "Veto module slot");
    }
  }
  /* compare controller cards */
  if (runInfo1.nCC == runInfo2.nCC) {
    for (i=0; i<runInfo1.nCC; i++) {
      ndiff += diff_int_2(runInfo1.CCcrate[i], runInfo2.CCcrate[i], i,
                          "CCcrate", "Controller card crate");
      ndiff += diff_int_2(runInfo1.CCslot[i], runInfo2.CCslot[i], i,
                          "CCslot", "Controller card slot");
    }
  }
  /* compare DAQ dataID types */
  if (runInfo1.idNum == runInfo2.idNum) {
    for (i=0; i<runInfo1.idNum; i++) {
      ndiff += diff_int_2(runInfo1.dataId[i], runInfo2.dataId[i], i,
                          "dataId", "Data-type ID");
      ndiff += diff_char(runInfo1.decoder[i], runInfo2.decoder[i], i,
                         "decoder", "Data-type decoder/label");
    }
  }

  if (nDets1 == nDets2) {
  /* compare Ge digitizer DAQ */
    for (i=0; i < nDets1; i++) {
      ndiff += diff_int_det(detInfo1, detInfo1[i].OrcaDetID, detInfo2[i].OrcaDetID, i,
                          "OrcaDetID", "Orca detector ID number");
      ndiff += diff_int_det(detInfo1, detInfo1[i].crate, detInfo2[i].crate, i,
                          "crate", "Detector digitizer crate");
      ndiff += diff_int_det(detInfo1, detInfo1[i].slot, detInfo2[i].slot, i,
                          "slot", "Detector digitizer slot");
      ndiff += diff_int_det(detInfo1, detInfo1[i].chanLo, detInfo2[i].chanLo, i,
                          "chanLo", "Detector digitizer LG channel");
      ndiff += diff_int_det(detInfo1, detInfo1[i].chanHi, detInfo2[i].chanHi, i,
                          "chanHi", "Detector digitizer HG channel");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HVCrate, detInfo2[i].HVCrate, i,
                          "HVCrate", "Detector HV crate");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HVCard, detInfo2[i].HVCard, i,
                          "HVCard", "Detector HV card");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HVChan, detInfo2[i].HVChan, i,
                          "HVChan", "Detector HV channel");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HVtarget, detInfo2[i].HVtarget, i,
                          "HVtarget", "Detector HV target voltage");
      ndiff += diff_int_det(detInfo1, detInfo1[i].DetType, detInfo2[i].DetType, i,
                          "DetType", "Detector type (nat/enr)");

      ndiff += diff_int_det(detInfo1, detInfo1[i].HGChEnabled, detInfo2[i].HGChEnabled, i,
                          "HGChEnabled", "HG chan enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGPreSumEnabled, detInfo2[i].HGPreSumEnabled, i,
                          "HGPreSumEnabled", "HG chan enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGPostrecnt, detInfo2[i].HGPostrecnt, i,
                          "HGPostrecnt", "HG chan Postrecnt");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGPrerecnt, detInfo2[i].HGPrerecnt, i,
                          "HGPrerecnt", "HG chan Prerecnt");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGTrigPolarity, detInfo2[i].HGTrigPolarity, i,
                          "HGTrigPolarity", "HG chan trigger polarity");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGTrigMode, detInfo2[i].HGTrigMode, i,
                          "HGTrigMode", "HG chan trigger mode");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGLEDThreshold, detInfo2[i].HGLEDThreshold, i,
                          "HGLEDThreshold", "HG chan LED threshold");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGTrapThreshold, detInfo2[i].HGTrapThreshold, i,
                          "HGTrapThreshold", "HG chan trap threshold");
      ndiff += diff_int_det(detInfo1, detInfo1[i].HGTrapEnabled, detInfo2[i].HGTrapEnabled, i,
                          "HGTrapEnabled", "HG chan trap enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGChEnabled, detInfo2[i].LGChEnabled, i,
                          "LGChEnabled", "LG chan enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGPreSumEnabled, detInfo2[i].LGPreSumEnabled, i,
                          "LGPreSumEnabled", "LG chan enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGPostrecnt, detInfo2[i].LGPostrecnt, i,
                          "LGPostrecnt", "LG chan Postrecnt");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGPrerecnt, detInfo2[i].LGPrerecnt, i,
                          "LGPrerecnt", "LG chan Prerecnt");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGTrigPolarity, detInfo2[i].LGTrigPolarity, i,
                          "LGTrigPolarity", "LG chan trigger polarity");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGTrigMode, detInfo2[i].LGTrigMode, i,
                          "LGTrigMode", "LG chan trigger mode");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGLEDThreshold, detInfo2[i].LGLEDThreshold, i,
                          "LGLEDThreshold", "LG chan LED threshold");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGTrapThreshold, detInfo2[i].LGTrapThreshold, i,
                          "LGTrapThreshold", "LG chan trap threshold");
      ndiff += diff_int_det(detInfo1, detInfo1[i].LGTrapEnabled, detInfo2[i].LGTrapEnabled, i,
                          "LGTrapEnabled", "LG chan trap enabled");

      ndiff += diff_int_det(detInfo1, detInfo1[i].pulseHighTime, detInfo2[i].pulseHighTime, i,
                          "pulseHighTime", "Pulser high-time");
      ndiff += diff_int_det(detInfo1, detInfo1[i].pulseLowTime, detInfo2[i].pulseLowTime, i,
                          "pulseLowTime", "Pulser low-time");
      ndiff += diff_int_det(detInfo1, detInfo1[i].pulserEnabled, detInfo2[i].pulserEnabled, i,
                          "pulserEnabled", "Pulser enabled");
      ndiff += diff_int_det(detInfo1, detInfo1[i].amplitude, detInfo2[i].amplitude, i,
                          "amplitude", "Pulser amplitude");
      ndiff += diff_int_det(detInfo1, detInfo1[i].attenuated, detInfo2[i].attenuated, i,
                          "attenuated", "Pulser attenuation");
      ndiff += diff_int_det(detInfo1, detInfo1[i].finalAttenuated, detInfo2[i].finalAttenuated, i,
                          "finalAttenuated", "Pulser final attenuation");

      ndiff += diff_char(detInfo1[i].DetName, detInfo2[i].DetName, i,
                         "DetName", "Detector name");
      ndiff += diff_char(detInfo1[i].StrName, detInfo2[i].StrName, i,
                         "StrName", "String name");
    }
  }

  if (ndiff == 1) {
    printf("\n --------- %d difference  --------- \n", ndiff);
  } else {
    printf("\n --------- %d differences --------- \n", ndiff);
  }
  return ndiff;
}

/* 
   decode_runfile_header.c

   code to read through file header and extract detector info etc from the XML
   and populate the detector info data structures

   David Radford  Nov 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

#define NMJPTS   16  // max number of pulser-tag channels in MJ model
#define NMJSTRS  14  // max number of detector strings in MJ model
#define NMJVSEGS 32  // max number of veto segments in MJ model
#define NGEDIGS  20  // max number of GRETINA digitizers
#define NGEHV   100  // max number of Ge HV channels
#define NGeCCS    8  // max number of preamp controller cards

/*  ------------------------------------------------------------ */
/*  Use these definitions to adjust the function of this program */
#define VERBOSE   0
#define VERB_DIG  0
#define FIX_PT1_SLOT 1    // apply fix to pulser tag 1 (crate 1 slot 8 -> slot 10)
/*  ------------------------------------------------------------ */


/*  -- -- -- first some convenience functions -- -- -- */

#define CHECK_FOR_STRING(A) if (!strstr(line, A)) {printf("ERROR: Missing %s\n %s\n", A, line); return 1;}

/* read_int(): function to read a single integer
   into dest
   from file f_in
   using reads into string line.   (NOTE: line must be at least 256 chars!) */
int read_int(FILE *f_in, int *dest, char *line) {
  char   *c;

  fgets(line, 256, f_in);
  if (strstr(line, "<false/>")) {
    *dest = 0;
  } else if (strstr(line, "<true/>")) {
    *dest = 1;
  } else if ((c=strstr(line, "<real>")) &&
             (1 == sscanf(c+6, "%d", dest))) {
  } else if (!(c=strstr(line, "<integer>")) ||
             (1 != sscanf(c+9, "%d", dest))) {
    fprintf(stderr, "\n ERROR decoding read_int:\n %s\n", line);
    return -1;
  }
  fgets(line, 256, f_in);  // read next line of input file, for use by calling function
  return 0;
}  /* read_int() */


/* read_int_array(): function to read array of integers or booleans
   of length dest_size
   into array dest[]
   from file f_in
   using reads into string line.   (NOTE: line must be at least 256 chars!) */
int read_int_array(FILE *f_in, int *dest, int dest_size, char *line) {
  int    i;
  float  f;
  char   *c;

  fgets(line, 256, f_in);
  if (!strstr(line, "<array>")) {
    fprintf(stderr, "\n ERROR in read_int_array! Missing <array>...\n %s\n", line);
    return -1;
  }
  fgets(line, 256, f_in);
  for (i = 0; i<dest_size; i++) {
    if (strstr(line, "<false/>")) {
      dest[i] = 0;
    } else if (strstr(line, "<true/>")) {
      dest[i] = 1;
    } else if ((c=strstr(line, "<real>"))) {
      if (1 != sscanf(c+6, "%f", &f)) {
        fprintf(stderr, "\n ERROR decoding read_int_array:\n %s\n", line);
        return -1;
      }
      dest[i] = lrintf(f);
    } else if (!(c=strstr(line, "<integer>")) ||
               (1 != sscanf(c+9, "%d", &dest[i]))) {
      fprintf(stderr, "\n ERROR decoding read_int_array:\n %s\n", line);
      return -1;
    }
    fgets(line, 256, f_in);
  }
  if (!strstr(line, "</array>")) {
    fprintf(stderr, "\n ERROR in read_int_array! Missing </array>... dest_size too small?\n %s\n", line);
    return -1;
  }

  fgets(line, 256, f_in);  // read next line of input file, for use by calling function
  // printf("read_int_array() returning %d values\n", i);
  return i;
}  /* read_int_array() */


/* read_float_array(): function to read array of floats
   of length dest_size
   into array dest[]
   from file f_in
   using reads into string line.   (NOTE: line must be at least 256 chars!) */
int read_float_array(FILE *f_in, float *dest, int dest_size, char *line) {
  int    i, j;
  char   *c;

  fgets(line, 256, f_in);
  if (!strstr(line, "<array>")) {
    fprintf(stderr, "\n ERROR in read_float_array! Missing <array>...\n %s\n", line);
    return -1;
  }
  fgets(line, 256, f_in);
  for (i = 0; i<dest_size; i++) {
    j=0;
    if ((c=strstr(line, "<integer>"))) j = sscanf(c+9, "%f", &dest[i]);
    if ((c=strstr(line, "<real>")))    j = sscanf(c+6, "%f", &dest[i]);
    if (j!=1) {
      fprintf(stderr, "\n ERROR decoding read_float_array [%d]:\n %s\n", i, line);
      return -1;
    }
    fgets(line, 256, f_in);
  }
  if (!strstr(line, "</array>")) {
    fprintf(stderr, "\n ERROR in read_float_array! Missing </array>... dest_size too small?\n %s\n", line);
    return -1;
  }

  fgets(line, 256, f_in);  // read next line of input file, for use by calling function
  // printf("read_float_array() returning %d values\n", i);
  return i;
}  /* read_float_array() */


/* check_false(): function to check array of booleans are all false
   and that string check is in first line.
   array length = size
   uses reads into string line   (NOTE: line must be at least 256 chars!)
   from file f_in */
int check_false(FILE *f_in, int size, char *line, char *check) {
  int    i;

  CHECK_FOR_STRING(check);
  fgets(line, 256, f_in);
  if (!strstr(line, "<array>")) {
    fprintf(stderr, "\n ERROR in check_false %s! Missing <array>...\n %s\n", check, line);
    return -1;
  }
  for (i = 0; i<size; i++) {
    fgets(line, 256, f_in);
    if (!strstr(line, "<false/>")) {
      fprintf(stderr, "\n check_false %s: %d NOT false\n %s\n", check, i, line);
      return -1;
    }
  }
  fgets(line, 256, f_in);
  if (!strstr(line, "</array>")) {
    fprintf(stderr, "\n ERROR in check_false %s! Missing </array>...\n %s\n", check, line);
    return -1;
  }
  fgets(line, 256, f_in);  // read next line of input file, for use by calling function
  return 0;
}  /* check_false() */


/* discard(): function to discard num lines from file f_in
   and make sure that string check is in first line.
   uses reads into string line   (NOTE: line must be at least 256 chars!) */
int discard(FILE *f_in, int num, char *line, char *check) {
  int    i;

  CHECK_FOR_STRING(check);
  for (i=0; i<num; i++) fgets(line, 256, f_in);
  return 0;
}  /* discard() */


/* ---------------------------------------------------------------------------
   decode_runfile_header():
   read through file header and extract detector info etc from the XML
   and populate the data structure array DetsReturn.
   input:   f_in (opened file pointer)
   output:  populated data structures DetsReturn[NMJDETS]
   returns: -1 on error
   .        otherwise the actual number of detectors found in the header
   --------------------------------------------------------------------------- */

int decode_runfile_header(FILE *f_in, MJDetInfo *DetsReturn, MJRunInfo *runInfo) {

  int    i, j, k, l, jj, kk, crateNum, slotNum, chNum, buf[25], reclen, reclen2;
  int    dataId[32], idNum = 0;
  char   *c, *c2, line[256];

  typedef struct{    // pulser-tag-channel info from MJ orca model
    char Description[32];
    int  VME, Slot, Chan, PADig, PAChan;
    // char Cable[32], Type[32];
  } MJModelPT;

  typedef struct{    // Detector string info from MJ orca model
    int  Index, Det[5];
    char Name[32];
  } MJModelStr;

  typedef struct{    // Veto segment info from MJ orca model
    int  kSegmentNumber, crate, slot, chan, HVCrate, HVCard, HVChan;
  } MJModelVeto;

  typedef struct{    // Detector HV info from OREHS8260pModel
    int  crate, slot, ch, target;
  } GeHVinfo;

  typedef struct{    // Ge WF digitizer info from ORGretina4MModel or ORGretina4AModel
    int  crate, slot, SerialNumber;
    int  type;       // 0 for ORGretina4M, 1 for ORGretina4A
    int  ChEnabled[10], LEDThreshold[10];
    //----------------- for type 0 = 4M
    int  ClockSource, ClockPhase, CollectionTime, IntegrationTime, DownSample;
    int  PreSumEnabled[10], Postrecnt[10], Prerecnt[10];
    int  TrigPolarity[10], TrigMode[10], TrapThreshold[10], TrapEnabled[10];
    //----------------- for type 1 = 4A
    int clkSelect0, clkSelect1, downSampleHoldOffTime, extDiscriminatorSrc, holdOffTime;
    int p2Window, rawDataLength, rawDataWindow, triggerConfig;
    int d3Window[10], dWindow[10], decimationFactor[10], discCountMode[10], discWidth[10];
    int kWindow[10], mWindow[10], p1Window[10], pileupMode[10];
  } GeDigInfo;

  typedef struct{    // ControlerCard info from ORMJDPreAmpModel
    int   crate, slot;     // serial control/readout link goes to this VME crate and slot
    int   adcEnabledMask, enabled[2], attenuated[2], finalAttenuated[2];
    int   preampID, pulseHighTime, pulseLowTime, pulserMask;
    int   amplitudes[16];
    float baselineVoltages[16];
    char  detectorNames[16][16];
    int   GeSlot[2];       // slots where the Ge _signals_ are digitized
  } CCInfo;

  MJDetInfo   MJMDets[NMJDETS], tmpDet;
  MJModelPT   MJMPTs[NMJPTS];
  MJModelStr  MJMStrs[NMJSTRS];
  MJModelVeto MJMVSegs[NMJVSEGS];
  GeHVinfo    GeHV[NGEHV];
  GeDigInfo   GeDig[NGEDIGS];
  CCInfo      GeCC[NGeCCS];
  int         nMJDets=0, nMJPTs=0, nMJStrs=0, nMJVSegs=0;
  int         nGeHV=0, nGeDig=0, nGeCC=0;

  /* initialize a few things */
  runInfo->dataIdGM = runInfo->dataIdGA = runInfo->idNum = runInfo->nV = 0;
  runInfo->flashcam = 0;

  if (strstr(runInfo->filename, ".fciops")) { // presorted flashcam data, no header
    runInfo->flashcam = 2;
    reclen = reclen2 = 0;
    runInfo->fileHeaderLen = 0;
  } else {
    /* read unformatted data and plist at start of orca file */
    fread(buf, sizeof(buf), 1, f_in);
    reclen = (unsigned int) buf[0];   // length of header in long words
    reclen2 = (unsigned int) buf[1];  // length of header in bytes
  }
  /* ----------------FlashCam data setup------------------ */
  if (runInfo->flashcam || reclen == -1000000001 || strstr(((char *)buf) + 4, "FlashCam")) {
    if (!runInfo->flashcam) runInfo->flashcam = 1;
    sprintf(MJMDets[0].DetName, "Det00");  // detector name
    nMJDets = 1;
    MJMDets[0].type            = 2;
    MJMDets[0].HGChEnabled     = 1;
    MJMDets[0].LGChEnabled     = 0;
    MJMDets[0].HGPreSumEnabled = MJMDets[0].LGPreSumEnabled = 0;
    MJMDets[0].HGPostrecnt     = MJMDets[0].LGPostrecnt     = 0;
    MJMDets[0].HGPrerecnt      = MJMDets[0].LGPrerecnt      = 0;
    /* copy results to returned data structures */ 
    for (i=0; i<nMJDets; i++)
      memcpy(&DetsReturn[i], &MJMDets[i], sizeof(MJDetInfo));
    runInfo->nGe = nMJDets;
    runInfo->nGD = 1;
    runInfo->nPT = 0;
    runInfo->nCC = 0;
    if (runInfo->flashcam < 2) {
      runInfo->fileHeaderLen = 128/4;
      fseek(f_in, 128, SEEK_SET);
    }
    return nMJDets;
  }

  /* loop through the lines of the XML data until we find the end of the plist */
  while (fgets(line, sizeof(line), f_in) && strncmp(line, "</plist>", 8)) {

    /* decode some run information */
    /* =========================== */
    if (strstr(line, "<key>date</key>")) {
      fgets(line, sizeof(line), f_in);
      CHECK_FOR_STRING("</string>");
      *(strstr(line, "</string>")) = '\0';
      strncpy(runInfo->date, strstr(line, "<string>")+8, 32);      // run date and time
    }
    if (strstr(line, "<key>documentName</key>")) {
      fgets(line, sizeof(line), f_in);
      CHECK_FOR_STRING("</string>");
      *(strstr(line, "</string>")) = '\0';
      strncpy(runInfo->orcaname, strstr(line, "<string>")+8, 256); // ORCA setup file name
    }
    if (strstr(line, "<key>RunNumber</key>"))
      if (read_int(f_in, &runInfo->runNumber, line)) return -1;    // run number
    if (strstr(line, "<key>quickStart</key>"))
      if (read_int(f_in, &runInfo->quickStart, line)) return -1;   // quick start flag
    if (strstr(line, "<key>refTime</key>"))
      if (read_int(f_in, &runInfo->refTime, line)) return -1;      // reference time? CHECK ME; enough bits?
    if (strstr(line, "<key>runType</key>"))
      if (read_int(f_in, &runInfo->runType, line)) return -1;      // run bits?
    if (strstr(line, "<key>startTime</key>"))
      if (read_int(f_in, &runInfo->startTime, line)) return -1;    // start time;     CHECK ME; enough bits?
      
    /* decode readout dataId's */
    /* ======================= */
    if (strstr(line, "<key>dataId</key>")) {
      if (idNum >= 32) {
        printf("ERROR: More than 32 dataIDs found in file header! Extra ones ingnored.\n\n");
      } else {
        while (!(c = strstr(line, "<integer>"))) fgets(line, sizeof(line), f_in);
        sscanf(c+9, "%d", &j);
        runInfo->dataId[idNum] = dataId[idNum] = j >> 18;
        while (!strstr(line, "<key>decoder</key>")) fgets(line, sizeof(line), f_in);
        while (!(c = strstr(line, "</string>"))) fgets(line, sizeof(line), f_in);
        *c = '\0';
        c = strstr(line, "<string>")+8;
        strncpy(runInfo->decoder[idNum], c, 32);
        if (VERBOSE) printf(" dataId = %2d  -> %s\n", j>>18, c);
        if (strstr(line, "Gretina4M")) runInfo->dataIdGM = dataId[idNum];
        if (strstr(line, "Gretina4A")) runInfo->dataIdGA = dataId[idNum];
        runInfo->idNum = ++idNum;
      }
    }

    /* decode DetectorGeometry */
    /* ======================= */
    if (strstr(line, "<key>DetectorGeometry</key>")) {
      fgets(line, sizeof(line), f_in);
      // encoded parameters given in next line:
      if (!strstr(line, "<string>kSegmentNumber,kVME,kCardSlot,kChanLo,kChanHi,kPreAmpChan,kHVCrate,"
                  "kHVCard,kHVChan,kMaxVoltage,kDetectorName,kDetectorType,kPreAmpDigitizer")) {
        fprintf(stderr, "\n ERROR in DetectorGeometry format!\n %s\n", line);
        return -1;
      }
      fgets(line, sizeof(line), f_in);

      while (!strstr(line, "</string>")) {
	//        if (!strstr(line, "--,--,--")) {  // empty data, so no detector present; skip this line
        if (!strstr(line, "--,--,--") &&
	    !strstr(line, ",,,,,,,,,,,")) {  // RLV Mostly empty data, so no detector present; skip this line
          if (nMJDets >= NMJDETS) {
            printf("ERROR: Too many detectors (>%d) found in ORCA MJD model!\n\n", NMJDETS);
            return -1;
          }
          c = strchr(line, ',');
          for (j=1; j<10; j++) c = strchr(c+1, ',');
          strncpy(MJMDets[nMJDets].DetName, c+1, sizeof(MJMDets[nMJDets].DetName));  // detector name
          c = strchr(c+1, ',');
        //if (10 != sscanf(line, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
          if (7 > sscanf(line, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",                    // various integers
                           &MJMDets[nMJDets].OrcaDetID, &MJMDets[nMJDets].crate,
                           &MJMDets[nMJDets].slot,      &MJMDets[nMJDets].chanLo,
                           &MJMDets[nMJDets].chanHi,    &j,                        // &MJMDets[nMJDets].PreAmpChan,
                           &MJMDets[nMJDets].HVCrate,   &MJMDets[nMJDets].HVCard,
                           &MJMDets[nMJDets].HVChan,    &MJMDets[nMJDets].HVMax) ||
              2  != sscanf(c+1, "%d,%d",
                           &MJMDets[nMJDets].DetType, &MJMDets[nMJDets].PreAmpDigitizer)) {
            fprintf(stderr, "\n ERROR decoding DetectorGeometry:\n %s\n", line);
            return -1;
          }
          if ((c = strchr(MJMDets[nMJDets].DetName, ','))) *c = 0;

          if (VERBOSE)  // report results
            printf(" Detector #%3d %2d %9s:  crate %d  slot %2d  chs %d %d"
                   "   PAch %2d   HV crate %d  slot %2d  ch %d   MaxV %d\n",
                   nMJDets, MJMDets[nMJDets].OrcaDetID,
                   MJMDets[nMJDets].DetName, MJMDets[nMJDets].crate,
                   MJMDets[nMJDets].slot,    MJMDets[nMJDets].chanLo,
                   MJMDets[nMJDets].chanHi,  j,                                   // MJMDets[nMJDets].PreAmpChan,
                   MJMDets[nMJDets].HVCrate, MJMDets[nMJDets].HVCard,
                   MJMDets[nMJDets].HVChan,  MJMDets[nMJDets].HVMax);
          MJMDets[nMJDets].PTcrate = MJMDets[nMJDets].PTslot = MJMDets[nMJDets].PTchan = -1;
          nMJDets++;
        }
        fgets(line, sizeof(line), f_in);
      }
      if (VERBOSE) printf("\n==> %d detectors found in model\n\n", nMJDets);
      if (nMJDets < NMJDETS) MJMDets[nMJDets].OrcaDetID = -1;
    }
    
    /* decode StringGeometry */
    /* ===================== */
    if (strstr(line, "<key>StringGeometry</key>")) {
      fgets(line, sizeof(line), f_in);
      // encoded parameters given in next line:
      if (!strstr(line, "<string>Index,Det1,Det2,Det3,Det4,Det5,Name")) {
        fprintf(stderr, "\n ERROR in StringGeometry format!\n %s\n", line);
        return -1;
      }
      fgets(line, sizeof(line), f_in);
      while (!strstr(line, "</string>")) {
        if (!strstr(line, "-,-,-,-,-,-")) {  // empty data, so no string present; skip this line
          if (nMJStrs >= NMJSTRS) {
            printf("ERROR: Too many detector strings (>%d) found in ORCA MJD model!\n\n", NMJSTRS);
            return -1;
          }
          for (j=0; j<5; j++) MJMStrs[nMJStrs].Det[j] = -1;
          c = strchr(line, ',');
          for (j=1; j<6; j++) c = strchr(c+1, ',');
          strncpy(MJMStrs[nMJStrs].Name, c+1, 4);                       // string name
          MJMStrs[nMJStrs].Name[4] = 0;
          if (2 >= sscanf(line, "%d,%d,%d,%d,%d,%d",                    // various integers
                          &MJMStrs[nMJStrs].Index,  &MJMStrs[nMJStrs].Det[0],
                          &MJMStrs[nMJStrs].Det[1], &MJMStrs[nMJStrs].Det[2],
                          &MJMStrs[nMJStrs].Det[3], &MJMStrs[nMJStrs].Det[4])) {
            fprintf(stderr, "\n ERROR decoding StringGeometry:\n %s\n", line);
            return -1;
          }

          if (VERBOSE) {  // report results
            printf(" String #%3d %2d %5s:  detectors",
                   nMJStrs, MJMStrs[nMJStrs].Index, MJMStrs[nMJStrs].Name);
            for (j=0; j<5 && MJMStrs[nMJStrs].Det[j] >= 0; j++)
              printf(" %2d", MJMStrs[nMJStrs].Det[j]);
            printf("\n");
          }
          nMJStrs++;
        }
        fgets(line, sizeof(line), f_in);
      }
      if (VERBOSE) printf("==> %d strings found in model\n\n", nMJStrs);
      if (nMJStrs < NMJSTRS) MJMStrs[nMJStrs].Index = -1;
    }
    
    /* decode VetoGeometry */
    /* =================== */
    if (strstr(line, "<key>VetoGeometry</key>")) {
      fgets(line, sizeof(line), f_in);
      // encoded parameters given in next line:
      if (!strstr(line, "<string>kSegmentNumber,kVME,kCardSlot,kChannel,kHVCrate,kHVCard,kHVChan")) {
        fprintf(stderr, "\n ERROR in VetoGeometry format!\n %s\n", line);
        return -1;
      }
      fgets(line, sizeof(line), f_in);
      while (!strstr(line, "</string>")) {
        if (!strstr(line, "--,--,--,--,--,--")) {  // empty data, so no segment present; skip this line
          if (nMJVSegs >= NMJVSEGS) {
            printf("ERROR: Too many veto segments (>%d) found in ORCA MJD model!\n\n", NMJVSEGS);
            return -1;
          }
          if (7 != sscanf(line, "%d,%d,%d,%d,%d,%d,%d",
                          &MJMVSegs[nMJVSegs].kSegmentNumber, &MJMVSegs[nMJVSegs].crate,
                          &MJMVSegs[nMJVSegs].slot,           &MJMVSegs[nMJVSegs].chan,
                          &MJMVSegs[nMJVSegs].HVCrate,        &MJMVSegs[nMJVSegs].HVCard,
                          &MJMVSegs[nMJVSegs].HVChan)) {
            fprintf(stderr, "\n ERROR decoding VetoGeometry:\n %s\n", line);
            return -1;
          }

          if (VERBOSE)  // report results
            printf(" Veto segment %2d %2d:  crate %d  slot %2d  ch %2d     HV crate %d  slot %2d  ch %2d\n",
                   nMJVSegs,
                   MJMVSegs[nMJVSegs].kSegmentNumber, MJMVSegs[nMJVSegs].crate,
                   MJMVSegs[nMJVSegs].slot,           MJMVSegs[nMJVSegs].chan,
                   MJMVSegs[nMJVSegs].HVCrate,        MJMVSegs[nMJVSegs].HVCard,
                   MJMVSegs[nMJVSegs].HVChan);
          nMJVSegs++;
        }
        fgets(line, sizeof(line), f_in);
      }
      if (VERBOSE) printf("==> %d veto segments found in model\n\n", nMJVSegs);
      if (nMJVSegs < NMJVSEGS) MJMVSegs[nMJVSegs].kSegmentNumber = -1;
    }

    /* decode pulser tagging channels */
    /* ============================== */
    if (strstr(line, "<key>SpecialStrings</key>")) {
      fgets(line, sizeof(line), f_in);
      // encoded parameters given in next line:
      if (!strstr(line, "<string>Index,Description,VME,Slot,Chan,PADig,PAChan")) {
        fprintf(stderr, "\n ERROR in decode_runfile_header; pulser tag SpecialStrings format!\n"
                " %s\n", line);
        return -1;
      }
      fgets(line, sizeof(line), f_in);
      while (!strstr(line, "</string>")) {
        if (!strstr(line, "-,-,-")) {  // empty data, so no PT present; skip this line
          if (nMJPTs >= NMJPTS) {
            printf("ERROR: Too many pulser tags (>%d) found in ORCA MJD model!\n\n", NMJPTS);
            return -1;
          }
          c = strchr(line, ',');
          c2 = strchr(c+1, ',');
          *c2 = '\0';
          strncpy(MJMPTs[nMJPTs].Description, c+1, sizeof(MJMPTs[nMJPTs].Description)); // text

          /* --------- special hack for missing VME number
             (example line: 0,M1-1 alpha3,-,8,4,6,5,red-blue,Internal,
             where the third value should be a VME crate) */
          if (*(c2+1) == '-') *(c2+1) = *(strstr(line,",M")+2);
          /* ---------- */
          if (5 != sscanf(c2+1, "%d,%d,%d,%d,%d",                                      // various integers
                          &MJMPTs[nMJPTs].VME,   &MJMPTs[nMJPTs].Slot, &MJMPTs[nMJPTs].Chan,
                          &MJMPTs[nMJPTs].PADig, &MJMPTs[nMJPTs].PAChan)) {
            fprintf(stderr, "\n ERROR decoding MotherBoard SpecialStrings:\n %s,%s\n", line, c2+1);
            // return -1;
            fgets(line, sizeof(line), f_in);
            continue;
          }
          if (VERBOSE)  // report results
            printf(" Pulser tag %2d %12s:  crate %2d  slot %2d  ch %2d     PADig %2d  ch %2d\n",
                   nMJPTs, MJMPTs[nMJPTs].Description,
                   MJMPTs[nMJPTs].VME,   MJMPTs[nMJPTs].Slot, MJMPTs[nMJPTs].Chan,
                   MJMPTs[nMJPTs].PADig, MJMPTs[nMJPTs].PAChan);
          nMJPTs++;
        }
        fgets(line, sizeof(line), f_in);
      }
      if (VERBOSE) printf("==> %d pulser-tag channels found in model\n\n", nMJPTs);
    }
    
    /* decode preamp controller cards */
    /* ============================== */
    if (strstr(line, "<string>ORMJDPreAmpModel</string>")) {
      if (nGeCC >= NGeCCS) {
        fprintf(stderr, "\n ERROR: Maximum number of Controller cards (%d) reached!\n\n", nGeCC);
        return -1;
      }
      fgets(line, sizeof(line), f_in);

      //This code is fixed to deal with an occasional missing "ConnectedTo" key is good data
      //The solution is not to return if missing.
      //if (discard(f_in, 2, line, "<key>ConnectedTo</key>")) return -1;
      if (strstr(line, "<key>ConnectedTo</key>")) {
        fgets(line, sizeof(line), f_in);
        fgets(line, sizeof(line), f_in);
	}
      CHECK_FOR_STRING("<key>adcEnabledMask</key>");
      if (read_int(f_in, &GeCC[nGeCC].adcEnabledMask, line)) return -1;
      CHECK_FOR_STRING("<key>amplitudes</key>");
      if (16 != read_int_array(f_in, &GeCC[nGeCC].amplitudes[0], 16, line)) return -1;
      CHECK_FOR_STRING("<key>attenuated0</key>");
      if (read_int(f_in, &GeCC[nGeCC].attenuated[0], line)) return -1;
      GeCC[nGeCC].attenuated[0] = 1 - GeCC[nGeCC].attenuated[0];  // seems to be encoded as inverse?
      CHECK_FOR_STRING("<key>attenuated1</key>");
      if (read_int(f_in, &GeCC[nGeCC].attenuated[1], line)) return -1;
      GeCC[nGeCC].attenuated[1] = 1 - GeCC[nGeCC].attenuated[1];  // seems to be encoded as inverse?
      CHECK_FOR_STRING("<key>baselineVoltages</key>");
      if (16 != read_float_array(f_in, &GeCC[nGeCC].baselineVoltages[0], 16, line)) return -1;
      if (discard(f_in, 2, line, "<key>boardRev</key>")) return -1;
      CHECK_FOR_STRING("<key>crate</key>");
      if (read_int(f_in, &GeCC[nGeCC].crate, line)) return -1;
      if (discard(f_in, 19, line, "<key>dacs</key>")) return -1;
      for (j=0; j<16; j++) {                           // ---- extract 16 detector names
        if (!(c=strstr(line, "<key>detectorName")) ||
            1 != sscanf(c+17, "%d", &i) || i < 0 || i > 15 ||
            !fgets(line, 256, f_in) ||
            !(c=strstr(line, "<string>")) || !(c2=strstr(line, "</string>"))) {
          fprintf(stderr, "\n ERROR decoding read_string:\n %s\n", line);
          return -1;
        }
        *c2 = 0;
        strncpy(GeCC[nGeCC].detectorNames[i], c+8, 16);
        if (VERBOSE && strlen(GeCC[nGeCC].detectorNames[i]) > 2)
          printf("CCDetName %d %2d %9s -> Baseline %6.2f\n",
                 nGeCC, i, GeCC[nGeCC].detectorNames[i], GeCC[nGeCC].baselineVoltages[i]);
        fgets(line, 256, f_in);
      }                                                // ----------------
      CHECK_FOR_STRING("<key>enabled0</key>");
      if (read_int(f_in, &GeCC[nGeCC].enabled[0], line)) return -1;
      GeCC[nGeCC].enabled[0] = 1 - GeCC[nGeCC].enabled[0];        // seems to be encoded as inverse?
      CHECK_FOR_STRING("<key>enabled1</key>");
      if (read_int(f_in, &GeCC[nGeCC].enabled[1], line)) return -1;
      GeCC[nGeCC].enabled[1] = 1 - GeCC[nGeCC].enabled[1];        // seems to be encoded as inverse?
      if (discard(f_in, 19, line, "<key>feedBackResistors</key>")) return -1;
      CHECK_FOR_STRING("<key>finalAttenuated0</key>");
      if (read_int(f_in, &GeCC[nGeCC].finalAttenuated[0], line)) return -1;
      GeCC[nGeCC].finalAttenuated[0] = 1 - GeCC[nGeCC].finalAttenuated[0];  // seems to be encoded as inverse?
      CHECK_FOR_STRING("<key>finalAttenuated1</key>");
      if (read_int(f_in, &GeCC[nGeCC].finalAttenuated[1], line)) return -1;
      GeCC[nGeCC].finalAttenuated[1] = 1 - GeCC[nGeCC].finalAttenuated[1];  // seems to be encoded as inverse?
      if (discard(f_in, 10, line, "<key>firmwareRev</key>")) return -1;
      CHECK_FOR_STRING("<key>pulseHighTime</key>");
      if (read_int(f_in, &GeCC[nGeCC].pulseHighTime, line)) return -1;
      CHECK_FOR_STRING("<key>pulseLowTime</key>");
      if (read_int(f_in, &GeCC[nGeCC].pulseLowTime, line)) return -1;
      CHECK_FOR_STRING("<key>pulserMask</key>");
      if (read_int(f_in, &GeCC[nGeCC].pulserMask, line)) return -1;
      CHECK_FOR_STRING("<key>slot</key>");
      if (read_int(f_in, &GeCC[nGeCC].slot, line)) return -1;
      if (VERBOSE)
        printf("    GeCC %d link is crate %d, slot %d\n", nGeCC, GeCC[nGeCC].crate, GeCC[nGeCC].slot);
      nGeCC++;
    } // end of preamp controller card info

    /* decode card info             */
    /* Ge HV and GRETINA digitizers */
    /* ============================ */
    /* first read GRETINA digitizer BLREnabled array
       - seems to be out of order? */
    if (strstr(line, "<key>Baseline Restore Enabled</key>")) {
      // just check that these are all false, don't bother to keep the values
      if (check_false(f_in, 10, line, "<key>Baseline Restore Enabled</key>")) return -1;
    }

    // read Card slot number
    if (strstr(line, "<key>Card</key>")) {
      fgets(line, sizeof(line), f_in);
      if (!(c = strstr(line, "<integer>"))) {
        fprintf(stderr, "\n ERROR in Card format! Missing <integer>...\n %s\n", line);
        return -1;
      }
      if (1 != sscanf(c+9, "%d", &slotNum)) {
        fprintf(stderr, "\n ERROR decoding Card slot number:\n %s\n", line);
        return -1;
      }
      fgets(line, sizeof(line), f_in);
      if (strstr(line, "<key>Chpsdv</key>")) {
        //GRETINA 4M v1.07 card (why does Chpdsv come before Class Name??)
        /* ------------- GRETINA 4M (LBL firmware) digitizer card info ------------ */
        if (VERBOSE) printf(" GRETINA 4M card found! Slot %d\n", slotNum);
        if (nGeDig >= NGEDIGS) {
          fprintf(stderr, "\n ERROR: Maximum number of GRETINA digitizers (%d) reached!\n\n", nGeDig);
          return -1;
        }

        GeDig[nGeDig].crate = -99;     // crate number not yet known; -99 is a space holder to be replaced later
        GeDig[nGeDig].slot = slotNum;
        if (discard(f_in, 13, line, "<key>Chpsdv</key>")) return -1;
        if (discard(f_in, 13, line, "<key>Chpsrt</key>")) return -1;
        CHECK_FOR_STRING("<key>Class Name</key>");
        fgets(line, sizeof(line), f_in);
        CHECK_FOR_STRING("<string>ORGretina4MModel</string>");
        fgets(line, sizeof(line), f_in);
        CHECK_FOR_STRING("<key>Clock Phase</key>");
        if (read_int(f_in, &GeDig[nGeDig].ClockPhase, line)) return -1;
        CHECK_FOR_STRING("<key>Clock Source</key>");
        if (read_int(f_in, &GeDig[nGeDig].ClockSource, line)) return -1;
        CHECK_FOR_STRING("<key>Collection Time</key>");
        if (read_int(f_in, &GeDig[nGeDig].CollectionTime, line)) return -1;
        if (check_false(f_in, 10, line, "<key>Debug Mode</key>")) return -1;
        CHECK_FOR_STRING("<key>Down Sample</key>");
        if (read_int(f_in, &GeDig[nGeDig].DownSample, line)) return -1;
        CHECK_FOR_STRING("<key>Enabled</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].ChEnabled[0], 10, line)) return -1;
        if (discard(f_in, 2, line, "<key>Ext Trig Length</key>")) return -1;
        if (discard(f_in, 2, line, "<key>External Window</key>")) return -1;
        if (discard(f_in, 13, line, "<key>FtCnt</key>")) return -1;
        if (discard(f_in, 2, line, "<key>Hist E Multiplier</key>")) return -1;
        CHECK_FOR_STRING("<key>Integration Time</key>");
        if (read_int(f_in, &GeDig[nGeDig].IntegrationTime, line)) return -1;
        CHECK_FOR_STRING("<key>LED Threshold</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].LEDThreshold[0], 10, line)) return -1;
        if (discard(f_in, 13, line, "<key>Mrpsdv</key>")) return -1;
        if (discard(f_in, 13, line, "<key>Mrpsrt</key>")) return -1;
        if (discard(f_in, 2, line, "<key>Noise Window</key>")) return -1;
        // just check that these are all false, don't bother to keep the values
        if (check_false(f_in, 10, line, "<key>PZ Trace Enabled</key>")) return -1;
        if (check_false(f_in, 10, line, "<key>Pile Up</key>")) return -1;
        if (discard(f_in, 2, line, "<key>Pile Up Window</key>")) return -1;
        if (check_false(f_in, 10, line, "<key>Pole Zero Enabled</key>")) return -1;
        if (discard(f_in, 13, line, "<key>Pole Zero Multiplier</key>")) return -1;
        CHECK_FOR_STRING("<key>Postrecnt</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].Postrecnt[0], 10, line)) return -1;
        CHECK_FOR_STRING("<key>PreSum Enabled</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].PreSumEnabled[0], 10, line)) return -1;
        CHECK_FOR_STRING("<key>Prerecnt</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].Prerecnt[0], 10, line)) return -1;
        if (strstr(line, "<key>Serial Number</key>")) {
          if (read_int(f_in, &GeDig[nGeDig].SerialNumber, line)) return -1;
          // printf("Digitizer %d has Serial Number %d\n", nGeDig, GeDig[nGeDig].SerialNumber);
          //discard(f_in, 2, line,  "<key>Serial Number</key>");
        }
        CHECK_FOR_STRING("<key>TPol</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].TrigPolarity[0], 10, line)) return -1;
        CHECK_FOR_STRING("<key>TRAP Threshold</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].TrapThreshold[0], 10, line)) return -1;
        CHECK_FOR_STRING("<key>Trap Enabled</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].TrapEnabled[0], 10, line)) return -1;
        CHECK_FOR_STRING("<key>Trigger Mode</key>");
        if (10 != read_int_array(f_in, &GeDig[nGeDig].TrigMode[0], 10, line)) return -1;
        if (discard(f_in, 2, line, "<key>baseAddress</key>")) return -1;
        if (discard(f_in, 13, line, "<key>forceFullInit</key>")) return -1;
        if (discard(f_in, 2, line, "<key>forceFullInitCard</key>")) return -1;
        GeDig[nGeDig].type = 0;
        nGeDig++;
      } // end of 4M digitizer card info

      if (strstr(line, "<key>Class Name</key>")) {
        fgets(line, sizeof(line), f_in);
        if (strstr(line, "<string>ORGretina4AModel</string>")) {
          /* ----------- GRETINA 4A card (ANL firmware) ------------- */
          if (VERBOSE) printf(" GRETINA 4A card found! Slot %d\n", slotNum);
          if (nGeDig >= NGEDIGS) {
            fprintf(stderr, "\n ERROR: Maximum number of GRETINA digitizers (%d) reached!\n\n", nGeDig);
            return -1;
          }
          GeDig[nGeDig].crate = -99;     // crate number not yet known; -99 is a space holder to be replaced later
          GeDig[nGeDig].slot = slotNum;
          fgets(line, sizeof(line), f_in);
          CHECK_FOR_STRING("<key>Clock Source</key>");
          if (read_int(f_in, &GeDig[nGeDig].ClockSource, line)) return -1;
          if (discard(f_in, 13, line, "<key>aHitCountMode</key>")) return -1;
          if (discard(f_in,  2, line, "<key>autoMode</key>")) return -1;
          if (discard(f_in,  2, line, "<key>auxIoConfig</key>")) return -1;
          if (discard(f_in,  2, line, "<key>auxIoRead</key>")) return -1;
          if (discard(f_in,  2, line, "<key>auxIoWrite</key>")) return -1;
          if (discard(f_in,  2, line, "<key>baseAddress</key>")) return -1;
          if (discard(f_in,  2, line, "<key>baselineDelay</key>")) return -1;
          if (discard(f_in, 13, line, "<key>baselineStart</key>")) return -1;
          if (discard(f_in,  2, line, "<key>channelPulsedControl</key>")) return -1;
          CHECK_FOR_STRING("<key>clkSelect0</key>");
          if (read_int(f_in, &GeDig[nGeDig].clkSelect0, line)) return -1;
          CHECK_FOR_STRING("<key>clkSelect1</key>");
          if (read_int(f_in, &GeDig[nGeDig].clkSelect1, line)) return -1;
          CHECK_FOR_STRING("<key>d3Window</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].d3Window[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>dWindow</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].dWindow[0], 10, line)) return -1;
          if (discard(f_in,  2, line, "<key>dacAttenuation</key>")) return -1;
          if (discard(f_in,  2, line, "<key>dacChannelSelect</key>")) return -1;
          CHECK_FOR_STRING("<key>decimationFactor</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].decimationFactor[0], 10, line)) return -1;
          if (discard(f_in,  2, line, "<key>diagChannelEventSel</key>")) return -1;
          if (discard(f_in,  2, line, "<key>diagInput</key>")) return -1;
          if (discard(f_in,  2, line, "<key>diagMuxControl</key>")) return -1;
          CHECK_FOR_STRING("<key>discCountMode</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].discCountMode[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>discWidth</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].discWidth[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>downSampleHoldOffTime</key>");
          if (read_int(f_in, &GeDig[nGeDig].downSampleHoldOffTime, line)) return -1;
          if (discard(f_in, 13, line, "<key>droppedEventCountMode</key>")) return -1;
          CHECK_FOR_STRING("<key>enabled</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].ChEnabled[0], 10, line)) return -1;
          if (discard(f_in, 13, line, "<key>eventCountMode</key>")) return -1;
          CHECK_FOR_STRING("<key>extDiscriminatorSrc</key>");
          if (read_int(f_in, &GeDig[nGeDig].extDiscriminatorSrc, line)) return -1;
          if (discard(f_in,  2, line, "<key>forceFullCardInit</key>")) return -1;
          if (discard(f_in, 13, line, "<key>forceFullInit</key>")) return -1;
          CHECK_FOR_STRING("<key>holdOffTime</key>");
          if (read_int(f_in, &GeDig[nGeDig].holdOffTime, line)) return -1;
          CHECK_FOR_STRING("<key>kWindow</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].kWindow[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>ledThreshold</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].LEDThreshold[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>mWindow</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].mWindow[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>p1Window</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].p1Window[0], 10, line)) return -1;
          CHECK_FOR_STRING("<key>p2Window</key>");
          if (read_int(f_in, &GeDig[nGeDig].p2Window, line)) return -1;
          if (discard(f_in,  2, line, "<key>peakSensitivity</key>")) return -1;
          if (discard(f_in, 13, line, "<key>pileExtensionMode</key>")) return -1;
          CHECK_FOR_STRING("<key>pileupMode</key>");
          if (10 != read_int_array(f_in, &GeDig[nGeDig].pileupMode[0], 10, line)) return -1;
          if (discard(f_in, 13, line, "<key>pileupWaveformOnlyMode</key>")) return -1;
          CHECK_FOR_STRING("<key>rawDataLength</key>");
          if (read_int(f_in, &GeDig[nGeDig].rawDataLength, line)) return -1;
          CHECK_FOR_STRING("<key>rawDataWindow</key>");
          if (read_int(f_in, &GeDig[nGeDig].rawDataWindow, line)) return -1;
          if (discard(f_in, 4*2, line, "<key>rj45SpareIoDir</key>")) return -1;
          CHECK_FOR_STRING("<key>triggerConfig</key>");
          if (read_int(f_in, &GeDig[nGeDig].triggerConfig, line)) return -1;
          if (discard(f_in, 4*2, line, "<key>userPackageData</key>")) return -1;
          GeDig[nGeDig].type = 1;
          nGeDig++;
        }  else if (strstr(line, "<string>OREHS8260pModel</string>")) { // Ge HV card
          /* decode Ge HV card info */
          /* ====================== */
          if (nGeHV >= NGEHV) {
            fprintf(stderr, "\n ERROR: Maximum number of Ge HV channels (%d) reached!\n\n", nGeHV);
            return -1;
          }

          while (!strstr(line, "<key>targets</key>")) fgets(line, sizeof(line), f_in);
          fgets(line, sizeof(line), f_in);
          if (!strstr(line, "<array>")) {
            fprintf(stderr, "\n ERROR in OREHS8260pModel format! Missing <array>...\n %s\n", line);
            return -1;
          }
          fgets(line, sizeof(line), f_in);
          for (chNum = 0; (c=strstr(line, "<integer>")); chNum++) {
            if (1 != sscanf(c+9, "%d", &GeHV[nGeHV].target)) {
              fprintf(stderr, "\n ERROR decoding HV target:\n %s\n", line);
              return -1;
            }
            GeHV[nGeHV].crate = -99;    // crate number not yet known; -99 is a space holder to be replaced later
            GeHV[nGeHV].slot = slotNum;
            GeHV[nGeHV].ch = chNum;

            nGeHV++;
            fgets(line, sizeof(line), f_in);
          }

        } else if (strstr(line, "<string>ORCaen792Model</string>") ||  // Veto modules
                   strstr(line, "<string>ORCV977Model</string>") ||
                   strstr(line, "<string>ORCV830Model</string>")) {
          /* decode veto DAQ card info */
          /* ====================== */
          runInfo->Vcrate[runInfo->nV] = -99;    // crate number not yet known
          runInfo->Vslot[runInfo->nV]  = slotNum;
          *strstr(line, "</string>") = 0;
          strncpy(runInfo->Vtype[runInfo->nV], strstr(line, ">OR")+1, 32);
          runInfo->nV++;
        }
      } // end of HV card info
    }

    /* decode crate info */
    /* ================ */
    if (strstr(line, "<key>CrateNumber</key>")) {
      fgets(line, sizeof(line), f_in);
      if (!(c = strstr(line, "<integer>"))) {
        fprintf(stderr, "\n ERROR in CrateNumber format! Missing <integer>...\n %s\n", line);
        return -1;
      }
      if (1 != sscanf(c+9, "%d", &crateNum)) {
        fprintf(stderr, "\n ERROR decoding CrateNumber:\n %s\n", line);
        return -1;
      }
      for (j=0; j<nGeHV; j++) {
        if (GeHV[j].crate == -99) {  // replace space holders now that we know the crate number
          GeHV[j].crate = crateNum;
          if (VERBOSE)  // report results
            printf(" GeHV %2d crate %d, slot %2d, ch: %d -> target %dV\n",
                   j, GeHV[j].crate, GeHV[j].slot, GeHV[j].ch, GeHV[j].target);
        }
      }
      for (j=0; j<nGeDig; j++) {
        if (GeDig[j].crate == -99) {  // replace space holders now that we know the crate number
          GeDig[j].crate = crateNum;
          if (VERBOSE)
            printf(" >>> GeDig # %2d crate %d, slot %2d\n", j, GeDig[j].crate, GeDig[j].slot);
        }
      }
      for (j=0; j<runInfo->nV; j++) {
        if (runInfo->Vcrate[j] == -99) {  // replace space holders now that we know the crate number
          runInfo->Vcrate[j] = crateNum;
          if (VERBOSE)
            printf(" >>> Veto module # %2d crate %d, slot %2d    %s\n",
                   j, runInfo->Vcrate[j], runInfo->Vslot[j], runInfo->Vtype[j]);
        }
      }
    }
  }  /* ============== end of loop reading XML ================= */

  /* report results summary */
  if (runInfo->dataIdGM == 0 && runInfo->dataIdGA == 0)
    printf("\n No data ID found for Gretina4M or 4A data!\n\n");
  if (VERBOSE && runInfo->dataIdGM)
    printf("\n Data ID %d found for Gretina4M data\n\n", runInfo->dataIdGM);
  if (VERBOSE && runInfo->dataIdGA)
    printf("\n Data ID %d found for Gretina4A data\n\n", runInfo->dataIdGA);
  if (!strstr(runInfo->argv[0], "rundiff")) {
    printf("%3d veto segments\n", nMJVSegs);
    printf("%3d Ge HV channels\n", nGeHV);
    printf("%3d Controller cards and %d Pulser-tag chs\n", nGeCC, nMJPTs);
    printf("%3d GRETINA WF digitizers\n", nGeDig);
    printf("%3d Detectors in %d strings\n\n", nMJDets, nMJStrs);
  }

  /* add CnPnDn names and re-order detector IDs in CnPnDn order */
  k = 0;
  for (i=0; i<nMJStrs; i++) {
    for (j=0; j<5 && MJMStrs[i].Det[j] >= 0; j++) {
      for (l = 0; l < nMJDets && MJMDets[l].OrcaDetID != MJMStrs[i].Det[j]; l++);
      if (l >= nMJDets) {
        printf("\n ERROR; Could not find detector number %d!\n", MJMStrs[i].Det[j]);
        return -1;
      }
      sprintf(MJMDets[l].StrName, "%sD%d", MJMStrs[i].Name, j+1);
      memcpy(&tmpDet, &MJMDets[k], sizeof(tmpDet));
      memcpy(&MJMDets[k], &MJMDets[l], sizeof(tmpDet));
      memcpy(&MJMDets[l], &tmpDet, sizeof(tmpDet));
      
      if (VERBOSE)  // report results
        printf(" Detector %2d  was %2d %8s %9s   HG ch %d,%2.2d,%d"
               "   HV %d,%2.2d,%d   MaxV %d\n",
               k, MJMDets[k].OrcaDetID,  MJMDets[k].StrName,
               MJMDets[k].DetName, MJMDets[k].crate,  MJMDets[k].slot,   MJMDets[k].chanHi,
               MJMDets[k].HVCrate, MJMDets[k].HVCard, MJMDets[k].HVChan, MJMDets[k].HVMax);
      k++;
    }
  }
  /* from DS7 onward, the detector table can have holes in it.  The sorting just done
     will remove the holes, because the empty entries have no name and no bias 
     entries.  Test replacing nMJDet by k. 
  */
  nMJDets = k;  //fixup empty entries in the raw detector table

  /*  --------------- give details of 4M digizer settings ------------------ */
  if (VERB_DIG) {
    for (jj=0; jj<nGeDig; jj++) {
      if (GeDig[jj].type != 0) continue;
      for (kk=0; kk<10; kk++) {
        if (GeDig[jj].ChEnabled[kk]) break;  // find first channel that is enabled
      }
      if (kk < 10) break;
    }
    printf("Common Gretina4M digitizer values:\n"
           "  ClockSource       %4d\n"
           "  CollectionTime    %4d        IntegrationTime   %4d\n"
           "  DownSample        %4d        PreSumEnabled[]   %4d\n"
           "  Postrecnt[]       %4d        Prerecnt[]        %4d\n"
           "  TrigPolarity[]    %4d        TrigMode[]        %4d\n"
           "  TrapEnabled[]     %4d\n",
           GeDig[jj].ClockPhase, GeDig[jj].CollectionTime, GeDig[jj].IntegrationTime,
           GeDig[jj].DownSample, GeDig[jj].PreSumEnabled[kk], GeDig[jj].Postrecnt[kk],
           GeDig[jj].Prerecnt[kk], GeDig[jj].TrigPolarity[kk], GeDig[jj].TrigMode[kk],
           GeDig[jj].TrapEnabled[kk]);
    printf("Variable digitizer value:   ClockPhase        %4d\n\n",
           GeDig[jj].ClockSource);
 
    for (j=0; j<nGeDig; j++) {  // slot
      if (GeDig[j].type != 0) continue;
      for (k=0; k<10; k++) {    // ch
        if (GeDig[j].ChEnabled[k]) {
          if (GeDig[j].ClockSource         != GeDig[jj].ClockSource        ||
              GeDig[j].CollectionTime      != GeDig[jj].CollectionTime     ||
              GeDig[j].IntegrationTime     != GeDig[jj].IntegrationTime    ||
              GeDig[j].DownSample          != GeDig[jj].DownSample         ||
              GeDig[j].PreSumEnabled[k]    != GeDig[jj].PreSumEnabled[kk]  ||
              GeDig[j].Postrecnt[k]        != GeDig[jj].Postrecnt[kk]      ||
              GeDig[j].Prerecnt[k]         != GeDig[jj].Prerecnt[kk]       ||
              GeDig[j].TrigPolarity[k]     != GeDig[jj].TrigPolarity[kk]   ||
              GeDig[j].TrigMode[k]         != GeDig[jj].TrigMode[kk]       ||
              GeDig[j].TrapEnabled[k]      != GeDig[jj].TrapEnabled[kk]) {
            printf("Digitizer in crate %d slot %2d ch %d is enabled,"
                   " but has different values than above!\n",
                   GeDig[j].crate, GeDig[j].slot, k);
            printf("\n Gretina4M Digitizer values:\n"
                   " ClockSource       %4d\n"
                   " CollectionTime    %4d        IntegrationTime   %4d\n"
                   " DownSample        %4d        PreSumEnabled[]   %4d\n"
                   " Postrecnt[]       %4d        Prerecnt[]        %4d\n"
                   " TrigPolarity[]    %4d        TrigMode[]        %4d\n"
                   " TrapEnabled[]     %4d\n\n",
                   GeDig[j].ClockSource, GeDig[j].CollectionTime, GeDig[j].IntegrationTime,
                   GeDig[j].DownSample, GeDig[j].PreSumEnabled[k], GeDig[j].Postrecnt[k],
                   GeDig[j].Prerecnt[k], GeDig[j].TrigPolarity[k], GeDig[j].TrigMode[k],
                   GeDig[j].TrapEnabled[k]);
          }
        }
      }
    }
  }

  /*  --------------- give details of 4A digizer settings ------------------ */
  if (VERB_DIG) {
    for (jj=0; jj<nGeDig; jj++) {
      if (GeDig[jj].type != 1) continue;
      for (kk=0; kk<10; kk++) {
        if (GeDig[jj].ChEnabled[kk]) break;  // find first channel that is enabled
      }
      if (kk < 10) break;
    }
    printf("Common Gretina4A digitizer values:\n"
           "  clkSelect0             %4d         clkSelect1           %4d\n"
           "  downSampleHoldOffTime  %4d         extDiscriminatorSrc  %4d\n"
           "  holdOffTime            %4d         triggerConfig        %4d\n"
           "  rawDataLength          %4d         rawDataWindow        %4d\n"
           "  decimationFactor       %4d\n"
           "  discCountMode          %4d         discWidth            %4d\n"
           "  LEDThreshold           %4d\n"
           "  mWindow                %4d         kWindow              %4d\n"
           "  d3Window               %4d         dWindow              %4d\n"
           "  p1Window               %4d         p2Window             %4d\n"
           "  pileupMode             %4d\n",
           GeDig[jj].clkSelect0, GeDig[jj].clkSelect1, GeDig[jj].downSampleHoldOffTime,
           GeDig[jj].extDiscriminatorSrc, GeDig[jj].holdOffTime, GeDig[jj].triggerConfig,
           GeDig[jj].rawDataLength, GeDig[jj].rawDataWindow, GeDig[jj].decimationFactor[kk],
           GeDig[jj].discCountMode[kk], GeDig[jj].discWidth[kk], GeDig[jj].LEDThreshold[kk],
           GeDig[jj].mWindow[kk], GeDig[jj].kWindow[kk], GeDig[jj].d3Window[kk],
           GeDig[jj].dWindow[kk], GeDig[jj].p1Window[kk], GeDig[jj].p2Window,
           GeDig[jj].pileupMode[kk]);
 
    for (j=0; j<nGeDig; j++) {  // slot
      if (GeDig[j].type != 1) continue;
      for (k=0; k<10; k++) {    // ch
        if (GeDig[j].ChEnabled[k]) {
          if (GeDig[j].clkSelect0            != GeDig[jj].clkSelect0            ||
              GeDig[j].clkSelect1            != GeDig[jj].clkSelect1            ||
              GeDig[j].downSampleHoldOffTime != GeDig[jj].downSampleHoldOffTime ||
              GeDig[j].extDiscriminatorSrc   != GeDig[jj].extDiscriminatorSrc   ||
              GeDig[j].holdOffTime           != GeDig[jj].holdOffTime           ||
              GeDig[j].p2Window              != GeDig[jj].p2Window              ||
              GeDig[j].rawDataLength         != GeDig[jj].rawDataLength         ||
              GeDig[j].rawDataWindow         != GeDig[jj].rawDataWindow         ||
              GeDig[j].triggerConfig         != GeDig[jj].triggerConfig         ||
              GeDig[j].d3Window[k]           != GeDig[jj].d3Window[kk]          ||
              GeDig[j].dWindow[k]            != GeDig[jj].dWindow[kk]           ||
              GeDig[j].decimationFactor[k]   != GeDig[jj].decimationFactor[kk]  ||
              GeDig[j].discCountMode[k]      != GeDig[jj].discCountMode[kk]     ||
              GeDig[j].discWidth[k]          != GeDig[jj].discWidth[kk]         ||
              GeDig[j].kWindow[k]            != GeDig[jj].kWindow[kk]           ||
              GeDig[j].LEDThreshold[k]       != GeDig[jj].LEDThreshold[kk]      ||
              GeDig[j].mWindow[k]            != GeDig[jj].mWindow[kk]           ||
              GeDig[j].p1Window[k]           != GeDig[jj].p1Window[kk]          ||
              GeDig[j].pileupMode[k]         != GeDig[jj].pileupMode[kk]){
            printf("Digitizer in crate %d slot %2d ch %d is enabled,"
                   " but has different values than above!\n",
                   GeDig[j].crate, GeDig[j].slot, k);
            printf("\n Gretina4A digitizer values:\n"
                   "  clkSelect0             %4d         clkSelect1           %4d\n"
                   "  downSampleHoldOffTime  %4d         extDiscriminatorSrc  %4d\n"
                   "  holdOffTime            %4d         triggerConfig        %4d\n"
                   "  rawDataLength          %4d         rawDataWindow        %4d\n"
                   "  decimationFactor       %4d\n"
                   "  discCountMode          %4d         discWidth            %4d\n"
                   "  LEDThreshold           %4d\n"
                   "  mWindow                %4d         kWindow              %4d\n"
                   "  d3Window               %4d         dWindow              %4d\n"
                   "  p1Window               %4d         p2Window             %4d\n"
                   "  pileupMode             %4d\n",
                   GeDig[j].clkSelect0, GeDig[j].clkSelect1, GeDig[j].downSampleHoldOffTime,
                   GeDig[j].extDiscriminatorSrc, GeDig[j].holdOffTime, GeDig[j].triggerConfig,
                   GeDig[j].rawDataLength, GeDig[j].rawDataWindow, GeDig[j].decimationFactor[k],
                   GeDig[j].discCountMode[k], GeDig[j].discWidth[k], GeDig[j].LEDThreshold[k],
                   GeDig[j].mWindow[k], GeDig[j].kWindow[k], GeDig[j].d3Window[k],
                   GeDig[j].dWindow[k], GeDig[j].p1Window[k], GeDig[j].p2Window,
                   GeDig[j].pileupMode[k]);
          }
        }
      }
    }

    printf("Digitizer channels that are enabled, but not assigned to Ge:\n"
           "  Crate  slot   ch\n");
    for (j=0; j<nGeDig; j++) {  // slot
      for (k=0; k<10; k++) {    // ch
        l = 0;
        if (GeDig[j].ChEnabled[k]) {
          l = 1;
          for (i=0; i<nMJDets; i++) {
            if (MJMDets[i].crate == GeDig[j].crate &&
                MJMDets[i].slot  == GeDig[j].slot) {
              if (k == MJMDets[i].chanHi || k == MJMDets[i].chanLo) {
                l = 0;
                break;
              }
            }
          }
        }
        if (l) {
          printf("%7d %5d %4d", GeDig[j].crate, GeDig[j].slot, k);
          for (i=0; i<nMJPTs; i++) {
            if (GeDig[j].crate == MJMPTs[i].VME &&
                GeDig[j].slot  == MJMPTs[i].Slot &&
                k == MJMPTs[i].Chan) printf("    -> pulser-tag number %d", i);
          }
          printf("\n");
        }
      }
    }
  } // --------------- end of digizer settings report --------------------

  /* this is needed for detector characterization data sets */
  if (nMJDets == 0) { // no detectors defined yet; look for enabled dig. chs to add as a new detector
    for (j=0; j<nGeDig; j++) {  // slot
      for (k=0; k<10; k++) {    // ch
        l = 0;
        if (GeDig[j].ChEnabled[k]) {
          l = 1;
          for (i=0; i<nMJDets; i++) {
            if (MJMDets[i].crate == GeDig[j].crate &&
                MJMDets[i].slot  == GeDig[j].slot) {
              if (k == MJMDets[i].chanHi || k == MJMDets[i].chanLo) {
                l = 0;
                break;
              }
            }
          }
        }
        if (l) {
          // add this channel as a new detector
          sprintf(MJMDets[nMJDets].DetName, "Det%2.2d", nMJDets);  // detector name
          sprintf(MJMDets[nMJDets].StrName, "Det%2.2d", nMJDets);  // detector name
          MJMDets[nMJDets].crate= GeDig[j].crate;
          MJMDets[nMJDets].slot = GeDig[j].slot;
          MJMDets[nMJDets].chanHi = k;
          MJMDets[nMJDets].chanLo = k+1;
          nMJDets++;
        }
      }
    }
  }

  /* Figure out which digitizer slots digitize the WFs
     from the preamps controlled by each of the controller cards  */
  for (j=0; j<nGeCC; j++) {
    GeCC[j].GeSlot[0] = GeCC[j].GeSlot[1] = -1;
  }

  //Set the GeSlots, checking for errors first
  for (j=0; j<nGeCC; j++) {

    for (k=0; k<16; k++) {

      for (i=0; i<nMJDets; i++) {
	
	if (strstr(GeCC[j].detectorNames[k], MJMDets[i].DetName)) {
          if (GeCC[j].crate != MJMDets[i].crate) {
	    //Check for crate==0, which implies that the ADC was removed from the Acq list
	    //Disable these channels which should not be in the data stream
	    if (GeCC[j].crate == 0) {
	      for (int kkk=0; kkk< nGeDig; kkk++) {
		if ((GeDig[kkk].crate == MJMDets[i].crate) &&
		    (GeDig[kkk].slot == MJMDets[i].slot)) {
		  GeDig[kkk].ChEnabled[MJMDets[i].chanHi]=0;
		  if (VERBOSE)
		    printf("NOTICE: Disable detector %s in crate %d slot %d because card removed from ACQ\n",
			   MJMDets[i].DetName, MJMDets[i].crate, MJMDets[i].slot);
		}
	      }
	    } else {
	      printf("\nERROR: CC %d crate %d != Detector %d crate %d!\n\n",
		     j, GeCC[j].crate, i,MJMDets[i].crate);
	      return -1;
	    }
          }
          if (GeCC[j].GeSlot[0] < 0) {
            GeCC[j].GeSlot[0] = MJMDets[i].slot;
          } else if (GeCC[j].GeSlot[0] != MJMDets[i].slot) {
            if (GeCC[j].GeSlot[1] < 0) {
              GeCC[j].GeSlot[1] = MJMDets[i].slot;
              if (VERBOSE)
                printf("CC %d waveforms -> Ge crate %d slots %2d and %2d\n",
                       j, GeCC[j].crate, GeCC[j].GeSlot[0], GeCC[j].GeSlot[1]);
            } else if (GeCC[j].GeSlot[1] != MJMDets[i].slot) {
              printf("\nERROR: More than two Ge WF digitizer slots for controller card %d!\n\n", j);
              return -1;
            }
          }
          break;
        }
      }
    }
    for (k=0; k<nMJPTs; k++) {
      if (MJMPTs[k].VME  == GeCC[j].crate &&
          MJMPTs[k].PADig == GeCC[j].slot) {
        if (VERBOSE)
          printf("CC %d has pulser tag number %2d for digitizer slots %2d, %2d in crate %d\n",
                 j, k, GeCC[j].GeSlot[0], GeCC[j].GeSlot[1], GeCC[j].crate);
        for (i=0; i<nMJDets; i++) {
          if (GeCC[j].GeSlot[0] == MJMDets[i].slot ||
              GeCC[j].GeSlot[1] == MJMDets[i].slot) {
            MJMDets[i].PTcrate = MJMPTs[i].VME;
            MJMDets[i].PTslot  = MJMPTs[i].Slot;
            MJMDets[i].PTchan  = MJMPTs[i].Chan;
          }
        }
        break;
      }
    }
    if (k == nMJPTs) {
      /* --------- special hack for PT # 1 slot number
         (pulser tag #1 (for M1) had PAChan entered as 8 instead of 10) */
      if (FIX_PT1_SLOT && j == 2 &&  // FIXME - add run number requirement?
          MJMPTs[1].PAChan == 5&&  
          MJMPTs[1].VME  == GeCC[j].crate &&
          MJMPTs[1].PADig == 8 &&
          GeCC[j].slot == 10) {
        MJMPTs[1].PADig = 10;
        if (VERBOSE) printf("  Fixing PAChan number for pulser tag 1 (8 -> 10)\n\n");
        j--;
        continue;
      }
      /* ---------- */

      printf(" Warning: Controller card %d has no matching pulser tag!\n", j);
    }
  }
    
  /* store actual HV target values, digitizer channel info, and CC card info
     ----------------------------------------------------------------------- */
  k = 0;
  for (i=0; i<nMJDets; i++) {
    // HV target values
    for (j=0; j<nGeHV; j++) {
      if (MJMDets[i].HVCrate == GeHV[j].crate &&
          MJMDets[i].HVCard  == GeHV[j].slot &&
          MJMDets[i].HVChan  == GeHV[j].ch) {
        MJMDets[i].HVtarget = GeHV[j].target;
        break;
      }
    }
   // GRETINA digitizer values
    for (j=0; j<nGeDig; j++) {
      if (MJMDets[i].crate == GeDig[j].crate &&
          MJMDets[i].slot  == GeDig[j].slot) {
        k = MJMDets[i].chanHi;
        MJMDets[i].DigSerialNum    = GeDig[j].SerialNumber;
        MJMDets[i].type            = GeDig[j].type;
        MJMDets[i].HGChEnabled     = GeDig[j].ChEnabled[k];
        MJMDets[i].HGLEDThreshold  = GeDig[j].LEDThreshold[k];
        if (GeDig[j].type == 0) {
          MJMDets[i].HGPreSumEnabled = GeDig[j].PreSumEnabled[k];
          MJMDets[i].HGPostrecnt     = GeDig[j].Postrecnt[k];
          MJMDets[i].HGPrerecnt      = GeDig[j].Prerecnt[k];
          MJMDets[i].HGTrigPolarity  = GeDig[j].TrigPolarity[k];
          MJMDets[i].HGTrigMode      = GeDig[j].TrigMode[k];
          MJMDets[i].HGLEDThreshold  = GeDig[j].LEDThreshold[k];
          MJMDets[i].HGTrapThreshold = GeDig[j].TrapThreshold[k];
          MJMDets[i].HGTrapEnabled   = GeDig[j].TrapEnabled[k];
        } else {
          MJMDets[i].HGTrigMode         = GeDig[j].triggerConfig;
          MJMDets[i].HGdecimationFactor = GeDig[j].decimationFactor[k];
          if (MJMDets[i].HGdecimationFactor > 1) MJMDets[i].HGPreSumEnabled = 1;
        }
        k = MJMDets[i].chanLo;
        MJMDets[i].LGChEnabled     = GeDig[j].ChEnabled[k];
        MJMDets[i].LGLEDThreshold  = GeDig[j].LEDThreshold[k];
        if (GeDig[j].type == 0) {
          MJMDets[i].LGPreSumEnabled = GeDig[j].PreSumEnabled[k];
          MJMDets[i].LGPostrecnt     = GeDig[j].Postrecnt[k];
          MJMDets[i].LGPrerecnt      = GeDig[j].Prerecnt[k];
          MJMDets[i].LGTrigPolarity  = GeDig[j].TrigPolarity[k];
          MJMDets[i].LGTrigMode      = GeDig[j].TrigMode[k];
          MJMDets[i].LGTrapThreshold = GeDig[j].TrapThreshold[k];
          MJMDets[i].LGTrapEnabled   = GeDig[j].TrapEnabled[k];
        } else {
          MJMDets[i].LGTrigMode         = GeDig[j].triggerConfig;
          MJMDets[i].LGdecimationFactor = GeDig[j].decimationFactor[k];
          if (MJMDets[i].LGdecimationFactor > 1) MJMDets[i].LGPreSumEnabled = 1;
        }
        break;
      }
    }
    if (j >= nGeDig) {
      printf(" Error: Detector %d has no associated digitizer!"
             "   Crate, slot = %d %d\n", i, MJMDets[i].crate, MJMDets[i].slot);
      MJMDets[i].DigSerialNum    = 0;
      MJMDets[i].type            = 0;
      MJMDets[i].HGChEnabled     = MJMDets[i].LGChEnabled     = 0;
      MJMDets[i].HGLEDThreshold  = MJMDets[i].LGLEDThreshold  = 0;
      MJMDets[i].HGPreSumEnabled = MJMDets[i].LGPreSumEnabled = 0;
      MJMDets[i].HGPostrecnt     = MJMDets[i].LGPostrecnt     = 0;
      MJMDets[i].HGPrerecnt      = MJMDets[i].LGPrerecnt      = 0;
      MJMDets[i].HGTrigPolarity  = MJMDets[i].LGTrigPolarity  = 0;
      MJMDets[i].HGTrigMode      = MJMDets[i].LGTrigMode      = 0;
      MJMDets[i].HGLEDThreshold  = MJMDets[i].LGLEDThreshold  = 0;
      MJMDets[i].HGTrapThreshold = MJMDets[i].LGTrapThreshold = 0;
      MJMDets[i].HGTrapEnabled   = MJMDets[i].LGTrapEnabled   = 0;
    }
    // Controller card values
    for (j=0; j<nGeCC; j++) {
      for (k=0; k<16; k++) {
        if (strstr(GeCC[j].detectorNames[k], MJMDets[i].DetName)) {
          MJMDets[i].CCnum           = j;
          MJMDets[i].pulseHighTime   = GeCC[j].pulseHighTime;
          MJMDets[i].pulseLowTime    = GeCC[j].pulseLowTime;
          MJMDets[i].pulserEnabled   = 0;
          if (GeCC[j].enabled[k/8] && (GeCC[j].pulserMask & 1<<k)) MJMDets[i].pulserEnabled = 1;
          MJMDets[i].amplitude       = GeCC[j].amplitudes[k];
          MJMDets[i].attenuated      = GeCC[j].attenuated[k/8];
          MJMDets[i].finalAttenuated = GeCC[j].finalAttenuated[k/8];
          MJMDets[i].baselineVoltage = GeCC[j].baselineVoltages[k];
        }
      }
    }

    /* write out final results */
    if (VERBOSE) {
      if (i%20 == 0)
        printf("#    DetID      pos      name     HiGain GAT Enab Thresh     HVch  MaxV Target"
               "     Pulser times enab  ampl atten\n");
      printf(" %3d  was %2d %8s %9s   %d,%2.2d,%d %4d %3d %6d %3d,%2.2d,%d"
             " %5d %6d %9d %7d %3d %5d %3d %d\n", i,
             MJMDets[i].OrcaDetID,       MJMDets[i].StrName,
             MJMDets[i].DetName,         MJMDets[i].crate,
             MJMDets[i].slot,            MJMDets[i].chanHi,
             MJMDets[i].crate*512 + MJMDets[i].slot*16 + MJMDets[i].chanHi, 
             MJMDets[i].HGChEnabled,     MJMDets[i].HGTrapThreshold,
             MJMDets[i].HVCrate,         MJMDets[i].HVCard,
             MJMDets[i].HVChan,          MJMDets[i].HVMax,
             MJMDets[i].HVtarget,
             MJMDets[i].pulseHighTime,   MJMDets[i].pulseLowTime,
             MJMDets[i].pulserEnabled,   MJMDets[i].amplitude,
             MJMDets[i].attenuated,      MJMDets[i].finalAttenuated);
    }
  }

  /* copy results to returned data structures */ 
  for (i=0; i<nMJDets; i++)
    memcpy(&DetsReturn[i], &MJMDets[i], sizeof(MJDetInfo));
  runInfo->nGe = nMJDets;

  for (i=0; i<nGeDig; i++) {
    runInfo->GDcrate[i] = GeDig[i].crate;
    runInfo->GDslot[i]  = GeDig[i].slot;
  }
  runInfo->nGD = nGeDig;

  for (i=0; i<nMJPTs; i++) {
    runInfo->PTcrate[i] = MJMPTs[i].VME;
    runInfo->PTslot[i]  = MJMPTs[i].Slot;
    runInfo->PTchan[i]  = MJMPTs[i].Chan;
  }
  runInfo->nPT = nMJPTs;

  for (i=0; i<nGeCC; i++) {
    runInfo->CCcrate[i] = GeCC[i].crate;
    runInfo->CCslot[i]  = GeCC[i].slot;
  }
  runInfo->nCC = nGeCC;

  if (VERBOSE)
    printf("File is at %d %ld; moving to reclen*4 = %d\n\n", reclen2+8, ftell(f_in), reclen*4);
  //fseek(f_in, reclen*4, SEEK_SET); // position file at start of events data
  fread(line, reclen*4-reclen2-8, 1, f_in);
  if (VERBOSE) printf("\n All Done.\n\n");
  runInfo->fileHeaderLen = reclen;

  /*
  if (nMJDets < 10) {
    for (i=0; i<NMJDETS; i++) {
      if (MJMDets[i].StrName[0] != 'C' || MJMDets[i].StrName[2] != 'P' || MJMDets[i].StrName[4] != 'D')
        sprintf(MJMDets[i].DetName, "Det%2.2d", i);
        sprintf(MJMDets[i].StrName, "Det%2.2d", i);
    }
  }
  */

  return nMJDets;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

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

#define ERR {printf("ERROR reading rms file at %ld!\n", ftell(file)); return 1;}
/* ======================================================================= */
int read_his(int *his, int idimsp, int spid, char *namesp, FILE *file)
{

  int  i, iy, nsp, dptr, spdir[2];
  int  id, spmode, xlen, expand;
  int  x0, nx, numch;
  char txt[64];
  char *modes[5] = {"shorts", "ints", "floats", "doubles", "bytes"};


  numch = 0;
  iy = spid;

  /* read spectrum ID directory */
  /* directory(ies):   2*(n+1) ints
     int n = number_of_IDs_in_header; int pointer_to_next_directory
     int sp_ID_1    int pointer_to_sp_data_1
     int sp_ID_2    int pointer_to_sp_data_2
     int sp_ID_3    int pointer_to_sp_data_3
     . . .
     int sp_ID_n    int pointer_to_sp_data_n
     [-1            0]  if no more spectra
  */
  //printf(" >> Looking for ID %d\n", iy); fflush(stdout);
  rewind(file);
  while (1) {
    if (1 != fread(&nsp,  sizeof(int), 1, file) ||
        1 != fread(&dptr, sizeof(int), 1, file)) ERR;
    for (i=0; i<nsp; i++) {  // read through directory entries
      if (2 != fread(spdir, sizeof(int), 2, file)) ERR;
      // printf("... found ID %d at %d\n", spdir[0], spdir[1]); fflush(stdout);
      if (spdir[0] == iy) break;  // found it!
    }
    if (spdir[0] == iy || dptr <= 0) break;
    fseek(file, dptr, SEEK_SET);
  }
  if (spdir[0] != iy) {
    printf("Spectrum ID %d not found!\n", iy);
    return 1;
  }

  /* if we are here, then we have found the spectrum we want */
  /* initialize output spectrum to zero */
  memset(his, 0, sizeof(int) * idimsp);
  if (spdir[1] < 1) return 1; // spectrum has zero length / does not exist

  fseek(file, spdir[1], SEEK_SET);
  /* now read the spectrum header */
  if (1 != fread(&id,      sizeof(int), 1, file) || // sp ID again
      1 != fread(&spmode,  sizeof(int), 1, file) || // storage mode (short/int/float...)
      1 != fread(&xlen,    sizeof(int), 1, file) || // x dimension
      1 != fread(&expand,  sizeof(int), 1, file) || // extra space for expansion
      1 != fread(txt,      sizeof(txt), 1, file)) ERR;  // description/name text
  if (id != iy) {
    printf("Spectrum ID mismatch (%d, %d)!\n", iy, id);
    return 1;
  }
  if (spmode != 2) {
    printf("spmode = %d\n", spmode);
    ERR;  // spectrum is not ints
  }

  if (0) printf("Sp ID %d; %d %s; %s\n", iy, xlen, modes[spmode-1], txt);
  if ((numch = xlen) > idimsp) {
    printf("First %d (of %d) chs only taken.", idimsp, numch);
    numch = idimsp;
  }

  i = 0;
  /* now read the actual spectrum data */
  while (i < numch) {
    if (1 != fread(&x0, sizeof(int), 1, file) ||    // starting ch for this chunk of spectrum
        1 != fread(&nx, sizeof(int), 1, file)) ERR; // number of chs in this chunk of spectrum
    if (x0 < 0 || nx < 1) break;   // no more non-zero bins in the spectrum
    if (nx > numch - x0) nx = numch - x0;

    if (nx != fread(his+x0, sizeof(int), nx, file)) ERR;
    i = x0 + nx;
  }
  strncpy(namesp, txt, 8);
  return 0;
#undef ERR
} /* rmsread */

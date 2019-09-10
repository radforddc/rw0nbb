#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>


/* ======================================================================= */
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
    nout = 0;
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

int rmsread(char *fn, int sp_id, int *sp, char *namesp, int *numch, int idimsp)
{
  /* subroutine to read spectrum of ints from radware multi-spectrum file
     into array sp, of dimension idimsp
     fn = file name
     sp_id = spectrum id to read
     numch = number of channels read
     namesp = title of spectrum (char*64)
     file extension must be .rms */


  int  i, iy, nsp, dptr, spdir[2];
  int  id, spmode, xlen, expand;
  int  x0, nx;
  char filnam[80];
  FILE *file;


  *numch = 0;
  strncpy(filnam, fn, sizeof(filnam));
  if (!(file = fopen(filnam, "r"))) return 1;

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
  iy = sp_id;
  while (1) {
    if (1 != fread(&nsp,  sizeof(int), 1, file) ||
        1 != fread(&dptr, sizeof(int), 1, file)) {
      printf("Read error0 for spectrum ID %d in file %s\n", iy, filnam);
      return 1;
    }
    // printf("file %s nsp = %d  dptr = %d\n", filnam, nsp, dptr);
    for (i=0; i<nsp; i++) {  // read through directory entries
      if (2 != fread(spdir, sizeof(int), 2, file)) {
        printf("Read error0 for spectrum ID %d in file %s\n", iy, filnam);
        return 1;
      }
      if (spdir[0] == iy) break;  // found it!
    }
    if (spdir[0] == iy || dptr <= 0) break;
    fseek(file, dptr, SEEK_SET);
  }
  if (spdir[0] != iy) {
    printf("Spectrum ID %d not found in file %s\n", iy, filnam);
    return 1;
  }

  /* if we are here, then we have found the spectrum we want */
  if (spdir[1] < 1) return 1; // spectrum has zero length / does not exist

  fseek(file, spdir[1], SEEK_SET);
  /* now read the spectrum header */
  if (1 != fread(&id,      sizeof(int), 1, file) || // sp ID again
      1 != fread(&spmode,  sizeof(int), 1, file) || // storage mode (short/int/float...)
      1 != fread(&xlen,    sizeof(int), 1, file) || // x dimension
      1 != fread(&expand,  sizeof(int), 1, file) || // extra space for expansion
      1 != fread(namesp,   64,          1, file)) {
    printf("Read error for spectrum ID %d in file %s\n", iy, filnam);
    return 1;
  }
  if (spmode != 2) {
    printf("Spectrum ID %d in file %s is not ints\n", iy, filnam);
    return 1;
  }
  if ((*numch = xlen) > idimsp) {
    printf("Read error for spectrum ID %d in file %s; too long(%d)\n", iy, filnam, xlen);
    *numch = idimsp;
    return 1;
  }

  /* initialize input spectrum to zero */
  memset(sp, 0, 4*xlen);

  /* now read the actual spectrum data */
  while (i < *numch) {
    if (1 != fread(&x0, sizeof(int), 1, file) ||    // starting ch for this chunk of spectrum
        1 != fread(&nx, sizeof(int), 1, file)) { // number of chs in this chunk of spectrum
      printf("Read error2 for spectrum ID %d in file %s\n", iy, filnam);
      return 1;
    }
    if (x0 < 0 || nx < 1) break;   // no more non-zero bins in the spectrum
    if (nx > *numch - x0) nx = *numch - x0;
    if (nx != fread(sp+x0, sizeof(int), nx, file)) {
      printf("Read error3 for spectrum ID %d in file %s\n", iy, filnam);
      return 1;
    }
    i = x0 + nx;
  }
    
  fclose(file);
  return 0;
} /* rmsread */


int main(int argc, char **argv)
{

  FILE *f;
  int   i, j, k, n = 0, n2 = 0;
  int   his[16384], his2[16384];
  char  txt[64];
  
  if (argc < 2 || strstr(argv[1], "-h")) {
    printf("\n  Usage: %s <input_file_1.rms> <input_file_2.rms> ...\n"
           "   Sums up both/all the his.rms files and\n"
           " creates a new file called sum.rms\n", argv[0]);
    return 0;
  }

  if (!(f=fopen(argv[1], "r"))) {
    printf("\n  File1 %s does not exist.\n", argv[1]);
    return 0;
  }
  fclose(f);
  if (!(f=fopen(argv[2], "r"))) {
    printf("\n  File2 %s does not exist.\n", argv[2]);
    return 0;
  }
  fclose(f);
  if (!(f=fopen("sum.rms", "w"))) {
    printf("\n  Cannot open file sum.rms\n");
    return 0;
  }

  for (i=0; i < 4096; i++) {
    if (rmsread(argv[1], i, his, txt, &n, 16384) || n < 10) break;
    for (j = 2; j < argc; j++) {
      if (rmsread(argv[j], i, his2, txt, &n2, 16384) || n2 != n) break;
      for (k=0; k<n; k++) his[k] += his2[k];
    }
    write_his(his, n, i, txt, f);
  }

  printf("\n Wrote sums of %d spectra in %d files to sum.rms\n\n", i, argc-1);

  return 0;
}

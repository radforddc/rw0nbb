#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>


/* ======================================================================= */
/* ----- write histograms to file, to look at with radware (gf3) ----- */
int write_his(int *his, int nch, int sp_id, char *sp_name, FILE *file) {

  static int dptr1 = 0, sptr1 = 0;
  static int dptr2 = 0, sptr2 = 0;
  int        dptr, sptr, id;
  static FILE *fsave1 = 0, *fsave2 = 0;
  typedef struct {
    int num;
    int nextptr;
    int dir[100][2];
  } rmsdir;
  static rmsdir dir1, dir2, *dir;
  static struct {
    int id;
    int mode;
    int xlen;
    int expand;
    char txt[64];
  } shead;
  int i, j=0, k, dd[2];


  if (fsave1 == 0) {  // new file; initialize
    dir1.num = 0;                      // up to 100 sp IDs per directory
    dir1.nextptr = -1;                 // pointer to next directory
    for (i=0; i<100; i++) dir1.dir[i][0] = dir1.dir[i][1] = -1;
    fwrite(&dir1, sizeof(dir1), 1, file);
    sptr1 = sizeof(dir1);               // pointer to where to write next spectrum
    fsave1 = file;
    printf("New file 1\n");
  } else if (file != fsave1 && fsave2 == 0) {  // new file; initialize
    dir2.num = 0;                      // up to 100 sp IDs per directory
    dir2.nextptr = -1;                 // pointer to next directory
    for (i=0; i<100; i++) dir2.dir[i][0] = dir2.dir[i][1] = -1;
    fwrite(&dir2, sizeof(dir2), 1, file);
    sptr2 = sizeof(dir2);               // pointer to where to write next spectrum
    fsave2 = file;
    printf("New file 2\n");
  } else if (file != fsave1 && file != fsave2) {
    printf("ERROR: file is not equal to fsave1 or fsave2!\n");
    return 1;
  }

  id = sp_id;
  if (file == fsave1) {
    dptr = dptr1;
    sptr = sptr1;
    dir = &dir1;
    printf("Writing ID %d to file 1; dptr, sptr, num = %d, %d, %d\n", id, dptr,sptr, dir->num);
  } else {
    dptr = dptr2;
    sptr = sptr2;
    dir = &dir2;
    printf("Writing ID %d to file 2; dptr, sptr, num = %d, %d, %d\n", id, dptr,sptr, dir->num);
  }

  /* find first non-zero channel */
  for (k=0; k<nch; k++) if (his[k] != 0) break;
  printf("First non-zero ch = %d\n", k);

  dir->dir[dir->num][0] = id;
  dir->dir[dir->num][1] = sptr;

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
  dir->num++;                         // up to 100 sp IDs per directory
  dir->nextptr = -1;                  // pointer to next directory
  if (dir->num == 100) {
    dir->nextptr = dptr = sptr;
  }
  fwrite(dir, sizeof(rmsdir), 1, file); // do re-write

  if (dir->num == 100) {
    fseek(file, dptr, SEEK_SET);      // go to the start of NEW directory
    dir->num = 0;
    dir->nextptr = -1;                // pointer to next directory
    for (i=0; i<100; i++) dir->dir[i][0] = dir->dir[i][1] = -1;
    fwrite(dir, sizeof(rmsdir), 1, file);
    sptr += sizeof(rmsdir);              // pointer to where to write next spectrum
  }

  if (file == fsave1) {
    sptr1 = sptr;                 // pointer to where to write next spectrum
    dptr1 = dptr;
    
  } else {
    sptr2 = sptr;                 // pointer to where to write next spectrum
    dptr2 = dptr;
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

  FILE  *f, *f1, *f2;
  int   i, n = 0;
  int   his[16384], z[16384] = {0};
  char  txt[64], fn[256];
  int   nat[100] = {0}, veto[100] = {0};
  int   nat_id[100] = {4,8,12,13,14,15,16,17,25,29,33,34,35,36,37,41,42,43,44,45,46,54,57, -1};
  // int   veto_id[100] = {8, 28, 30, 31, 49, 50, -1};  // DS5a
  // int   veto_id[100] = {8, 30, 32, 50, -1};          // DS5b
  int   veto_id[100] = {8, 30, 31, 32, 50, -1};  // DS5c
  
  if (argc < 2 || strstr(argv[1], "-h")) {
    printf("\n  Usage: %s <input_file.rms> \n", argv[0]);
    return 0;
  }

  if (!(f=fopen(argv[1], "r"))) {
    printf("\n  File %s does not exist.\n", argv[1]);
    return 0;
  }
  fclose(f);
  for (i=0; i < 100 && nat_id[i] >= 0; i++) nat[nat_id[i]] = 1;
  for (i=0; i < 100 && veto_id[i] >= 0; i++) veto[veto_id[i]] = 1;

  sprintf(fn, "enr_%s", argv[1]);
  if (!(f1=fopen(fn, "w"))) {
    printf("\n  Cannot open file %s\n", fn);
    return 0;
  }
  sprintf(fn, "nat_%s", argv[1]);
  if (!(f2=fopen(fn, "w"))) {
    printf("\n  Cannot open file %s\n", fn);
    return 0;
  }

  for (i=0; i < 4096; i++) {
    if (rmsread(argv[1], i, his, txt, &n, 16384) || n < 10) break;
    if (veto[i%100]) {
      write_his(z, n, i, txt, f1);
      write_his(z, n, i, txt, f2);
    } else if (nat[i%100]) {
      write_his(z, n, i, txt, f1);
      write_his(his, n, i, txt, f2);
    } else {
      write_his(his, n, i, txt, f1);
      write_his(z, n, i, txt, f2);
    }
  }

  printf("\n Done; %d spectra\n\n", i);

  return 0;
}

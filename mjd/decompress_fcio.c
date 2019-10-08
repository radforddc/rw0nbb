#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define VERBOSE 0

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


int main(int argc, char **argv) {

  FILE     *f_in, *f_out;
  int      n = 0, len;
  unsigned short    buf[16384];
  short    buf2[16384];
  char     fname[256], fname2[256], line[256];

  if (argc < 2) {
    fprintf(stderr,
            "\nusage: %s <input_file_name>\n\n", argv[0]);
    return -1;
  }
  strncpy(fname, argv[1], sizeof(fname));
  strncpy(fname2, argv[1], sizeof(fname2));
  strncat(fname2, ".decompress", sizeof(fname2) - strlen(fname) - 1);

  /* open raw data file as input */
  if ((f_in = fopen(fname,"r")) == NULL) {
    fprintf(stderr, "\n Failed to open input file %s\n", fname);
    return 0;
  }
  if ((f_out = fopen(fname2,"w")) == NULL) {
    fprintf(stderr, "\n Failed to open output file %s\n", fname2);
    return 0;
  }
  printf("\n >>> Reading %s, writing %s\n\n", fname, fname2);

  /* read file header (128 bytes) */
  fread(line, 128, 1, f_in);
  if (!(strstr(line+4, "FlashCamV1"))) {
    fprintf(stderr,
            "This file does not appear to be a FlashCamV1 data file!\n"
            "This record should say FlashCamV1: %s\n\n", line+4);
    return 0;
  }
  fwrite(line, 128, 1, f_out);

  while (1) {
    /* read record header */
    if (fread(&len,  4, 1, f_in) != 1) {
      printf("\n No more data, %d records read\n\n", n);
      fclose(f_in);
      fclose(f_out);
      return 0;
    }
    if (len < 1) {
      fwrite(&len,  4, 1, f_out);
      continue;
    }
    if (len > 32768) {
      printf("\n Error after %d records read: len = %d\n\n", n, len);
      fclose(f_in);
      fclose(f_out);
      return 0;
    }
    if (VERBOSE) printf("record %d; len = %d\n", n, len);
    /* read waveform */
    if (fread(buf, len, 1, f_in) != 1) {
      printf("\n No more data, %d records read\n\n", n);
      fclose(f_in);
      fclose(f_out);
      return 0;
    }
    if (len > 40) {
      len = 2 * decompress_signal(buf, buf2, len/2);
      if (++n % 10000 == 0) {printf("\r %6d records...", n); fflush(stdout);}
      fwrite(&len,  4, 1, f_out);
      fwrite(buf2, len, 1, f_out);
    } else {
      fwrite(&len,  4, 1, f_out);
      fwrite(buf, len, 1, f_out);
    }
  }
  return -1;
}
    

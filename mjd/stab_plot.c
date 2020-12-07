/* code to make a combined data file, which will in turn be used to create a series of plots
    to check stability of PZ.input, ctc.input, gains.input, and psa.input parameters 

    example:
    cd ~/analysis/mjd/ds6b/cal
    ls -1d ds* > j
    stab_plot j 05      (makes stab05.dat, for detector ID 5))

    cp stab05.dat stab_det.dat; doplot1.sh stab_det; sed -e "s/XYZ/ 5   C1P2D2/g" stab_det.ps > stab05.ps
    cat stab[0-5][0-9].ps > stab_all.ps
    ps2pdf stab_all.ps; open stab_all.pdf
    rm stab_det.dat stab_det.ps stab[0-5][0-9].ps stab[0-5][0-9].dat stab_all.ps
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {

  FILE   *f_in, *f_in2, *f_out;
  char   fn_in[256], fn_out[256], line[256];
  char   *fn2[5] = {"gains.input", "ctc.input", "PZ.input", "fwhm_ctc.txt", "psa.input"};
  int    i, j, k, det_id=1;
  double gain[1000]={1.0}, fwhm;
  int    data_good[1000] = {1};

  if (argc < 3) {
    fprintf(stderr,
            "\nusage: %s <dir_list_fname_in> <det_num>\n\n",
            argv[0]);
    return -1;
  }
  if (!(f_in = fopen(argv[1], "r"))) {
    printf("File %s does not exist??\n", argv[1]);
    return -1;
  }
  det_id = atoi(argv[2]);

  sprintf(fn_out, "stab%2.2d.dat", det_id);
  f_out = fopen(fn_out, "w");

  for (j=0; j<5; j++) {
    i = 1;
    while (fgets(fn_in, 256, f_in)) {
      sprintf(fn_in + strlen(fn_in) - 1, "/%s", fn2[j]);
      if (!(f_in2 = fopen(fn_in, "r"))) {
        printf("Error: File %s does not exist?", fn_in);
        return 0;
      }

      fgets(line, 256, f_in2);
      if (i == 1 && line[0] == '#') fputs(line, f_out);  // # header line
      if (i == 1 && j == 3) fprintf(f_out, "#  Chan     gain    sigma     fwhm\n");
      k = -1;
      while (k < det_id && fgets(line, 256, f_in2)) {
        if (line[0] == '#') continue;
        k = atoi(line);
      }
      if (k > det_id && (j == 0 || data_good[i])) {
        printf("Data missing for detector ID %2d, %3d %s\n", det_id, i, fn_in);
      } else {
        if (j == 0) {                                // reading gains.input
          gain[i] = 1.0;
          data_good[i] = 1;
          sscanf(line, "%d %lf", &k, &gain[i]);
          // printf(" >> %lf %lf %2d %s | %s\n", gain[i-1], gain[i], i, fn_in, line);
          if (gain[i] > 0.99 ||                         // if gain == 1 or gain is exactly the same as previous calib
              (0 && i > 1 && gain[i] == gain[i-1])) {   //    then this detector is missing from this data subset
            printf("Data missing for detector ID %2d %3d, %s\n", det_id, i, fn_in);
            data_good[i] = 0;
          }
        } else if (data_good[i]) {    // gain data is valid for this data subset and detector
          if (j == 3) {
            sscanf(line, "%d %lf", &k, &fwhm);
            fprintf(f_out, " %2d %4d %9.5lf %8.5lf %6.2lf\n",
                    i, det_id, gain[i], gain[i]*fwhm/2.355/2614.5, fwhm);
          } else if ((j == 1 && !strstr(line, " 1.00   1.00000000 ")) ||    // ctc.input (very reliable for identifying missing data)
                     (j == 2 && !strstr(line, "72.5000")) ||                // PZ.input
                     (j == 4 && !strstr(line, " 0.00  500.00  500.00 "))) { // psa.input
            fprintf(f_out, " %2d %s", i, line);
          } else {
            printf("Data missing for detector ID %2d, %3d %s\n", det_id, i, fn_in);
            data_good[i] = 0;
          }
        }
      }
      i++;
      fclose(f_in2);
    }
    if (j > 0) fprintf(f_out, "\n");
    printf("%d %s files read\n", --i, fn2[j]);
    rewind(f_in);
  }
  fclose(f_out);
  fclose(f_in);
  return 0;
}
   

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

  FILE *f_in, *f_in2, *f_out;
  char fn_in[256], fn_out[256], line[256], *fn2[4] = {"PZ", "ctc", "gains", "psa"};
  int  i, j, k, det_id=1;

  if (argc < 3) {
    fprintf(stderr,
            "\nusage: %s <fname_in> <det_num>\n\n",
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

  for (j=0; j<4; j++) {
    i = 1;
    while (fgets(fn_in, 256, f_in)) {
      sprintf(fn_in + strlen(fn_in) - 1, "/%s.input", fn2[j]);
      if (!(f_in2 = fopen(fn_in, "r"))) {
        printf("Error: File %s does not exist?", fn_in);
        return 0;
      }

      fgets(line, 256, f_in2);
      if (i == 1 && line[0] == '#') fputs(line, f_out);  // # header line
      while (fgets(line, 256, f_in2)) {
        if (line[0] == '#') continue;
        k = atoi(line);
        if (k == det_id) {
          fprintf(f_out, " %2d %s", i++, line);
          break;
        }
      }
      fclose(f_in2);
    }
    fprintf(f_out, "\n");
    printf("%d %s.input files read\n", --i, fn2[j]);
    rewind(f_in);
  }
  fclose(f_out);

  sprintf(fn_out, "stab%2.2d.pdc", det_id);
  //f_out = fopen(fn_out, "w");

  //fclose(f_out);
  fclose(f_in);
  return 0;
}
   

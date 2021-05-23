#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>

#include "MJDSort.h"

/*
 * Looks at veto_only_runs.txt and ds.../a[12].txt files to estimate fractional exposure
 *  To run, do: check_veto_only <veto_only_runs.txt> <evl.txt> [additional_det_IDs]
 */

#define VERBOSE 0

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  int           i, j, k;
  char          line[1024], subset_dir[256];
  FILE          *fp1, *fp2, *f_out;
  struct dirent *de;  // Pointer for directory entry
  DIR           *dp;
  float         a1[100], a2[100], a3[100];

  /* get physics run directories */
  if ((dp = opendir(".")) == NULL) {
    printf("Could not open current directory" );
    return 0;
  }
  while ((de = readdir(dp)) != NULL) {
    if (strstr(de->d_name, "ds") == de->d_name)
      strncpy(subset_dir, de->d_name, sizeof(subset_dir));
    sprintf(line, "%s/a1.txt",subset_dir);
    if (!(fp1 = fopen(line, "r"))) {
      // printf("%s does not exist?\n", line);
      continue;
    }
    printf(" %s + ", line);
    sprintf(line, "%s/a2.txt",subset_dir);
    if (!(fp2 = fopen(line, "r"))) {
      printf("\n%s does not exist?\n", line);
      return 0;
    }
    printf("%s", line);
    sprintf(line, "%s/a3.txt",subset_dir);
    f_out = fopen(line, "w");
    fprintf(f_out, "# ch ratio sum\n");
    printf(" -> %s\n", line);

    for (i=0; i < 100; i++)  a1[i] = a2[i] = a3[i] = 0;
    i = j = 0;
    while (fgets(line, sizeof(line), fp1)) {
      if (strstr(line, "8192 chs read") &&
          sscanf(line, " Sp. %d", &j) == 1 &&
          j > 800 && j < 900) {
        i = j-800;
        k = 0;
      } else if (strstr(line, "Chs 500 to 6000    Area = 0")) {
        j = 0;
        if (k == 0) {
          a1[i] = j;
          k++;
        } else {
          a3[i] = a1[i] + j;
        }
      } else if (strstr(line, "Chs 500 to 6000,  Area: "))  {
        sscanf(line+24, "%d", &j);
        if (k == 0) {
          a1[i] = j;
          k++;
        } else {
          a3[i] = a1[i] + j;
        }
      }
    }

    while (fgets(line, sizeof(line), fp2)) {
      if (strstr(line, "8192 chs read") &&
          sscanf(line, " Sp. %d", &j) == 1 &&
          j > 800 && j < 900) {
        i = j-800;
        k = 0;
      } else if (strstr(line, "Chs 500 to 6000    Area = 0")) {
        j = 0;
        if (k == 0) {
          a2[i] = j;
          k++;
        } else if (a3[i] != a2[i] + j) {
          printf("ERROR: i = %d; a1, a2, a3, j: %f %f %f %d\n", i, a1[i], a2[i], a3[i], j);
          return(0);
        }
      } else if (strstr(line, "Chs 500 to 6000,  Area: "))  {
        sscanf(line+24, "%d", &j);
        if (k == 0) {
          a2[i] = j;
          k++;
        } else if (a3[i] != a2[i] + j) {
          printf("ERROR: i = %d; a1, a2, a3, j: %f %f %f %d\n", i, a1[i], a2[i], a3[i], j);
          return(0);
        }
      }
    }
    fclose(fp1);
    fclose(fp2);

    /* calculate ratio of events that are cut by extra data cleaning */
    for (i=0; i < 100; i++) {
      if (a3[i] > 0.5) {
        a1[i] = (a1[i] - a2[i]) / a3[i];
        if (a1[i] >0.0001) fprintf(f_out, "%3d %7.3f  %.0f\n", i, a1[i], a3[i]);
      }
    }
    fclose(f_out);
  }


  closedir(dp);
  return 0;
}

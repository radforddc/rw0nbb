#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

/*
 * Looks at veto_only_runs.txt and evl.txt files to see if any events should be treated as veto-only
 *  To run, do: check_veto_only <veto_only_runs.txt> <evl.txt> [additional_det_IDs]
 */

#define VERBOSE 0

/*  ------------------------------------------------------------ */

int main(int argc, char **argv) {

  int        i, nveto[200] = {0}, runveto[200][100][2] = {{{0}}}, r1, r2, id;
  char       line[1024], *c, veto_str[200][100][10];
  FILE       *f_evl, *f_out, *f_veto;


  if (argc <3 ) {
    printf("\n Usage: %s  <veto_only_runs.txt> <evl.txt> [additional_det_IDs]\n", argv[0]);
    return 0;
  }

  /* open file for veto-only * detector list */
  if (!(f_veto = fopen(argv[1], "r"))) {
    printf("%s does not exist?\n", argv[1]);
    return 0;
  }
  /* open file for event list */
  if (!(f_evl = fopen(argv[2], "r"))) {
    printf("%s does not exist?\n", argv[2]);
    return 0;
  }
  f_out = fopen("evl_cvo.txt", "w");

  /* --------- read file for veto-only * detector list -------- */
  while (fgets(line, sizeof(line), f_veto)) {
    // 50   ds36710  ds37004L    26  C1P7D2 P42574B   t
    if (line[0] == '#') continue;
    c = strstr(line, "ds") + 2;
    sscanf(c, "%5d", &r1);
    c = strstr(c, "ds") + 2;
    sscanf(c, "%5d", &r2);
    c = strstr(c, " ") + 2;
    sscanf(c, "%d", &id);
    c = strstr(c, "C") + 15;

    runveto[id][nveto[id]][0] = r1;
    runveto[id][nveto[id]][1] = r2;
    strncpy(veto_str[id][nveto[id]], c, 9);
    nveto[id]++;
  }
  fclose(f_veto);
    
  /* --------- read file for event list -------- */
  while (fgets(line, sizeof(line), f_evl)) {
    if (line[0] == '#' || !(c = strstr(line, "Enr HG"))) {
      fprintf(f_out, "%s", line);
      continue;
    }
    sscanf(line, "%d", &id);
    sscanf(c+34, "%d", &r1);
    for (i=0; i<nveto[id]; i++) {
      if (runveto[id][i][0] < r1 && runveto[id][i][1] > r1) {
        sprintf(c, "Veto  %s", c+6);
        c += strlen(c);
        sprintf(c-1, "  %s\n", veto_str[id][i]);
        break;
      }
    }
    fprintf(f_out, "%s", line);
  }

  fclose(f_evl);
  fclose(f_out);
  return 0;
}

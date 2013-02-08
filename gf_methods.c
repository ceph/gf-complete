/*
 * gf_mult.c
 *
 * Multiplies two numbers in gf_2^w
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "gf_complete.h"
#include "gf_method.h"

#define NMULTS (14)
static char *mults[NMULTS] = { "SHIFT", "GROUP44", "GROUP48", "BYTWO_p", "BYTWO_b",
                               "TABLE", "LOG", "LOG_ZERO", "SPLIT2", "SPLIT4", "SPLIT8", "SPLIT88", "COMPOSITE-0", "COMPOSITE-1" };

#define NREGIONS (96) 
static char *regions[NREGIONS] = { "-", "SINGLE", "DOUBLE", "QUAD",
"LAZY", "SINGLE,LAZY", "DOUBLE,LAZY", "QUAD,LAZY", "SSE",
"SINGLE,SSE", "DOUBLE,SSE", "QUAD,SSE", "LAZY,SSE",
"SINGLE,LAZY,SSE", "DOUBLE,LAZY,SSE", "QUAD,LAZY,SSE", "NOSSE",
"SINGLE,NOSSE", "DOUBLE,NOSSE", "QUAD,NOSSE", "LAZY,NOSSE",
"SINGLE,LAZY,NOSSE", "DOUBLE,LAZY,NOSSE", "QUAD,LAZY,NOSSE",
"STDMAP", "SINGLE,STDMAP", "DOUBLE,STDMAP", "QUAD,STDMAP",
"LAZY,STDMAP", "SINGLE,LAZY,STDMAP", "DOUBLE,LAZY,STDMAP",
"QUAD,LAZY,STDMAP", "SSE,STDMAP", "SINGLE,SSE,STDMAP",
"DOUBLE,SSE,STDMAP", "QUAD,SSE,STDMAP", "LAZY,SSE,STDMAP",
"SINGLE,LAZY,SSE,STDMAP", "DOUBLE,LAZY,SSE,STDMAP",
"QUAD,LAZY,SSE,STDMAP", "NOSSE,STDMAP", "SINGLE,NOSSE,STDMAP",
"DOUBLE,NOSSE,STDMAP", "QUAD,NOSSE,STDMAP", "LAZY,NOSSE,STDMAP",
"SINGLE,LAZY,NOSSE,STDMAP", "DOUBLE,LAZY,NOSSE,STDMAP",
"QUAD,LAZY,NOSSE,STDMAP", "ALTMAP", "SINGLE,ALTMAP", "DOUBLE,ALTMAP",
"QUAD,ALTMAP", "LAZY,ALTMAP", "SINGLE,LAZY,ALTMAP",
"DOUBLE,LAZY,ALTMAP", "QUAD,LAZY,ALTMAP", "SSE,ALTMAP",
"SINGLE,SSE,ALTMAP", "DOUBLE,SSE,ALTMAP", "QUAD,SSE,ALTMAP",
"LAZY,SSE,ALTMAP", "SINGLE,LAZY,SSE,ALTMAP",
"DOUBLE,LAZY,SSE,ALTMAP", "QUAD,LAZY,SSE,ALTMAP", "NOSSE,ALTMAP",
"SINGLE,NOSSE,ALTMAP", "DOUBLE,NOSSE,ALTMAP", "QUAD,NOSSE,ALTMAP",
"LAZY,NOSSE,ALTMAP", "SINGLE,LAZY,NOSSE,ALTMAP",
"DOUBLE,LAZY,NOSSE,ALTMAP", "QUAD,LAZY,NOSSE,ALTMAP", "CAUCHY",
"SINGLE,CAUCHY", "DOUBLE,CAUCHY", "QUAD,CAUCHY", "LAZY,CAUCHY",
"SINGLE,LAZY,CAUCHY", "DOUBLE,LAZY,CAUCHY", "QUAD,LAZY,CAUCHY",
"SSE,CAUCHY", "SINGLE,SSE,CAUCHY", "DOUBLE,SSE,CAUCHY",
"QUAD,SSE,CAUCHY", "LAZY,SSE,CAUCHY", "SINGLE,LAZY,SSE,CAUCHY",
"DOUBLE,LAZY,SSE,CAUCHY", "QUAD,LAZY,SSE,CAUCHY", "NOSSE,CAUCHY",
"SINGLE,NOSSE,CAUCHY", "DOUBLE,NOSSE,CAUCHY", "QUAD,NOSSE,CAUCHY",
"LAZY,NOSSE,CAUCHY", "SINGLE,LAZY,NOSSE,CAUCHY",
"DOUBLE,LAZY,NOSSE,CAUCHY", "QUAD,LAZY,NOSSE,CAUCHY" };

#define NDIVS (3)
static char *divides[NDIVS] = { "-", "MATRIX", "EUCLID" }; 

int main()
{
  int m, r, d, w, i, sa, j;
  char *argv[20];
  gf_t gf;
  char divs[200], ks[10], ls[10];

  methods_to_stderr();

  printf("\n");
  printf("Implemented Methods: \n\n");
  
  for (i = 2; i < 8; i++) {
    w = (1 << i);
    argv[0] = "-";
    if (create_gf_from_argv(&gf, w, 1, argv, 0) > 0) {
      printf("w=%d: -\n", w);
      gf_free(&gf, 1);
    }
    for (m = 0; m < NMULTS; m++) {
      sa = 0;
      if (strcmp(mults[m], "GROUP44") == 0) {
        argv[sa++] = "GROUP";
        argv[sa++] = "4";
        argv[sa++] = "4";
      } else if (strcmp(mults[m], "GROUP48") == 0) {
        argv[sa++] = "GROUP";
        argv[sa++] = "4";
        argv[sa++] = "8";
      } else if (strcmp(mults[m], "SPLIT2") == 0) {
        argv[sa++] = "SPLIT";
        sprintf(ls, "%d", w);
        argv[sa++] = ls;
        argv[sa++] = "2";
      } else if (strcmp(mults[m], "SPLIT4") == 0) {
        argv[sa++] = "SPLIT";
        sprintf(ls, "%d", w);
        argv[sa++] = ls;
        argv[sa++] = "4";
      } else if (strcmp(mults[m], "SPLIT8") == 0) {
        argv[sa++] = "SPLIT";
        sprintf(ls, "%d", w);
        argv[sa++] = ls;
        argv[sa++] = "8";
      } else if (strcmp(mults[m], "SPLIT88") == 0) {
        argv[sa++] = "SPLIT";
        argv[sa++] = "8";
        argv[sa++] = "8";
      } else if (strcmp(mults[m], "COMPOSITE-0") == 0) {
        argv[sa++] = "COMPOSITE";
        argv[sa++] = "2";
        argv[sa++] = "0";
        argv[sa++] = "-";
      } else if (strcmp(mults[m], "COMPOSITE-1") == 0) {
        argv[sa++] = "COMPOSITE";
        argv[sa++] = "2";
        argv[sa++] = "1";
        argv[sa++] = "-";
      } else {
        argv[sa++] = mults[m];
      }
      for (r = 0; r < NREGIONS; r++) {
        argv[sa++] = regions[r]; 
        strcpy(divs, "");
        for (d = 0; d < NDIVS; d++) {
          argv[sa++] = divides[d];
/*          printf("w=%d:", w);
          for (j = 0; j < sa; j++) printf(" %s", argv[j]);
          printf("\n"); */
          if (create_gf_from_argv(&gf, w, sa, argv, 0) > 0) {
            strcat(divs, "|");
            strcat(divs, divides[d]);
            gf_free(&gf, 1);
          } 
          sa--;
        }
        if (strlen(divs) > 0) {
          printf("w=%d:", w);
          for (j = 0; j < sa; j++) printf(" %s", argv[j]);
          printf(" %s\n", divs+1);
        }
        sa--;
      }
      sa--;
    }
  }
}

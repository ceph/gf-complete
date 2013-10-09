/*
 * gf_methods.c
 *
 * Lists supported methods (incomplete w.r.t. GROUP and COMPOSITE 
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "gf_complete.h"
#include "gf_method.h"
#include "gf_int.h"

#define NMULTS (16)
static char *mults[NMULTS] = { "SHIFT", "CARRY_FREE", "GROUP44", "GROUP48", "BYTWO_p", "BYTWO_b",
                               "TABLE", "LOG", "LOG_ZERO", "LOG_ZERO_EXT", "SPLIT2",
                               "SPLIT4", "SPLIT8", "SPLIT16", "SPLIT88", "COMPOSITE" };

#define NREGIONS (7) 
static char *regions[NREGIONS] = { "DOUBLE", "QUAD", "LAZY", "SSE", "NOSSE", 
                                   "ALTMAP", "CAUCHY" };

#define NDIVS (2)
static char *divides[NDIVS] = { "MATRIX", "EUCLID" }; 

int main() 
{
  int m, r, d, w, i, sa, j, k, reset;
  char *argv[50];
  gf_t gf;
  char divs[200], ks[10], ls[10];
  
  for (i = 2; i < 8; i++) {
    w = (1 << i);
    argv[0] = "-";
    if (create_gf_from_argv(&gf, w, 1, argv, 0) > 0) {
      printf("w=%d: -\n", w);
      gf_free(&gf, 1);
    } else if (_gf_errno == GF_E_DEFAULT) {
      fprintf(stderr, "Unlabeled failed method: w=%d: -\n", 2);
      exit(1);
    }

    for (m = 0; m < NMULTS; m++) {
      sa = 0;
      argv[sa++] = "-m";
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
      } else if (strcmp(mults[m], "SPLIT16") == 0) {
        argv[sa++] = "SPLIT";
        sprintf(ls, "%d", w);
        argv[sa++] = ls;
        argv[sa++] = "16";
      } else if (strcmp(mults[m], "SPLIT88") == 0) {
        argv[sa++] = "SPLIT";
        argv[sa++] = "8";
        argv[sa++] = "8";
      } else if (strcmp(mults[m], "COMPOSITE") == 0) {
        argv[sa++] = "COMPOSITE";
        argv[sa++] = "2";
        argv[sa++] = "-";
      } else {
        argv[sa++] = mults[m];
      }
      reset = sa;
      for (r = 0; r < (1 << NREGIONS); r++) {
        sa = reset;
        for (k = 0; k < NREGIONS; k++) {
          if (r & 1 << k) {
            argv[sa++] = "-r";
            argv[sa++] = regions[k];
          }
        }
        argv[sa++] = "-";
        if (create_gf_from_argv(&gf, w, sa, argv, 0) > 0) {
          printf("w=%d:", w);
          for (j = 0; j < sa; j++) printf(" %s", argv[j]);
          printf("\n");
          gf_free(&gf, 1);
        } else if (_gf_errno == GF_E_DEFAULT) {
          fprintf(stderr, "Unlabeled failed method: w=%d:", w);
          for (j = 0; j < sa; j++) fprintf(stderr, " %s", argv[j]);
          fprintf(stderr, "\n");
          exit(1);
        }
        sa--;
        for (d = 0; d < NDIVS; d++) {
          argv[sa++] = "-d";
          argv[sa++] = divides[d];
          /*          printf("w=%d:", w);
                      for (j = 0; j < sa; j++) printf(" %s", argv[j]);
                      printf("\n"); */
          argv[sa++] = "-";
          if (create_gf_from_argv(&gf, w, sa, argv, 0) > 0) {
            printf("w=%d:", w);
            for (j = 0; j < sa; j++) printf(" %s", argv[j]);
            printf("\n");
            gf_free(&gf, 1);
          } else if (_gf_errno == GF_E_DEFAULT) {
            fprintf(stderr, "Unlabeled failed method: w=%d:", w);
            for (j = 0; j < sa; j++) fprintf(stderr, " %s", argv[j]);
            fprintf(stderr, "\n");
            exit(1);
          } 
          sa-=3;
        }
      }
    }
  }
  return 0;
}

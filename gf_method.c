/*
 * gf_method.c
 *
 * Parses argv to figure out the mult_type and arguments.  Returns the gf.
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "gf_complete.h"
#include "gf_method.h"

void methods_to_stderr()
{
  fprintf(stderr, "To specify the methods, do one of the following: \n");
  fprintf(stderr, "       - leave empty to use defaults\n");
  fprintf(stderr, "       - use a single dash to use defaults\n");
  fprintf(stderr, "       - specify MULTIPLY REGION DIVIDE\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Legal values of MULTIPLY:\n");
  fprintf(stderr, "       SHIFT: shift\n");
  fprintf(stderr, "       GROUP g_mult g_reduce: the Group technique - see the paper\n");
  fprintf(stderr, "       BYTWO_p: BYTWO doubling the product.\n");
  fprintf(stderr, "       BYTWO_b: BYTWO doubling b (more efficient thatn BYTWO_p)\n");
  fprintf(stderr, "       TABLE: Full multiplication table\n");
  fprintf(stderr, "       LOG:   Discrete logs\n");
  fprintf(stderr, "       LOG_ZERO: Discrete logs with a large table for zeros\n");
  fprintf(stderr, "       SPLIT g_a g_b: Split tables defined by g_a and g_b\n");
  fprintf(stderr, "       COMPOSITE k rec METHOD: Composite field.  GF((2^l)^k), l=w/k.\n");
  fprintf(stderr, "                               rec = 0 means inline single multiplication\n");
  fprintf(stderr, "                               rec = 1 means recursive single multiplication\n");
  fprintf(stderr, "                               METHOD is the method of the base field in GF(2^l)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Legal values of REGION: Specify multiples with commas e.g. 'DOUBLE,LAZY'\n");
  fprintf(stderr, "       -: Use defaults\n");
  fprintf(stderr, "       SINGLE/DOUBLE/QUAD: Expand tables\n");
  fprintf(stderr, "       LAZY: Lazily create table (only applies to TABLE and SPLIT)\n");
  fprintf(stderr, "       SSE/NOSSE: Use 128-bit SSE instructions if you can\n");
  fprintf(stderr, "       CAUCHY/ALTMAP/STDMAP: Use different memory mappings\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Legal values of DIVIDE:\n");
  fprintf(stderr, "       -: Use defaults\n");
  fprintf(stderr, "       MATRIX: Use matrix inversion\n");
  fprintf(stderr, "       EUCLID: Use the extended Euclidian algorithm.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "See the user's manual for more information.\n");
  fprintf(stderr, "There are many restrictions, so it is better to simply use defaults in most cases.\n");
}

int create_gf_from_argv(gf_t *gf, int w, int argc, char **argv, int starting)
{
  int mult_type, divide_type, region_type;
  uint32_t prim_poly = 0;
  int arg1, arg2, subrg_size;
  gf_t *base;
  char *crt, *x, *y;

  if (argc <= starting || strcmp(argv[starting], "-") == 0) {
    mult_type = GF_MULT_DEFAULT;
    if (!gf_init_easy(gf, w, mult_type)) return 0;
    return (argc <= starting) ? starting : starting+1;
  }

  region_type = GF_REGION_DEFAULT;
  divide_type = GF_DIVIDE_DEFAULT;

  arg1 = 0;
  arg2 = 0;
  prim_poly = 0;
  base = NULL;
  subrg_size = 0;
  
  if (argc < starting+3) return 0;

  if (strcmp(argv[starting], "SHIFT") == 0) {
    mult_type = GF_MULT_SHIFT;
    starting++;
  } else if (strcmp(argv[starting], "GROUP") == 0) {
    mult_type = GF_MULT_GROUP;
    if (argc < starting+5) return 0;
    if (sscanf(argv[starting+1], "%d", &arg1) == 0 ||
        sscanf(argv[starting+2], "%d", &arg2) == 0 ||
        arg1 <= 0 || arg2 <= 0 || arg1 >= w || arg2 >= w) return 0;
    starting += 3;
  } else if (strcmp(argv[starting], "BYTWO_p") == 0) {
    mult_type = GF_MULT_BYTWO_p;
    starting++;
  } else if (strcmp(argv[starting], "BYTWO_b") == 0) {
    mult_type = GF_MULT_BYTWO_b;
    starting++;
  } else if (strcmp(argv[starting], "TABLE") == 0) {
    mult_type = GF_MULT_TABLE;
    starting++;
  } else if (strcmp(argv[starting], "LOG") == 0) {
    mult_type = GF_MULT_LOG_TABLE;
    starting++;
  } else if (strcmp(argv[starting], "LOG_ZERO") == 0) {
    mult_type = GF_MULT_LOG_TABLE;
    arg1 = 1;
    starting++;
  } else if (strcmp(argv[starting], "SPLIT") == 0) {
    mult_type = GF_MULT_SPLIT_TABLE;
    if (argc < starting+5) return 0;
    if (sscanf(argv[starting+1], "%d", &arg1) == 0 ||
        sscanf(argv[starting+2], "%d", &arg2) == 0 ||
        arg1 <= 0 || arg2 <= 0 || w % arg1 != 0 || w % arg2 != 0) return 0;
    starting += 3;
  } else if (strcmp(argv[starting], "COMPOSITE") == 0) {
    mult_type = GF_MULT_COMPOSITE;
    if (argc < starting+6) return 0;
    if (sscanf(argv[starting+1], "%d", &arg1) == 0 ||
        sscanf(argv[starting+2], "%d", &arg2) == 0 ||
        arg1 <= 1 || w %arg1 != 0 || ((arg2 | 1) != 1)) return 0;
    base = (gf_t *) malloc(sizeof(gf_t));
    starting = create_gf_from_argv(base, w/arg1, argc, argv, starting+3);
    if (starting == 0) { free(base); return 0; }
  } else {
    return 0;
  }

  if (argc < starting+2) {
    if (base != NULL) gf_free(base, 1);
    return 0;
  }

  if (strcmp(argv[starting], "-") == 0) {
    region_type = GF_REGION_DEFAULT;
  } else {
    crt = strdup(argv[starting]);
    region_type = 0;
    x = crt;
    do { 
      y = strchr(x, ','); 
      if (y != NULL) *y = '\0';
      if (strcmp(x, "DOUBLE") == 0) {
        region_type |= GF_REGION_DOUBLE_TABLE;
      } else if (strcmp(x, "QUAD") == 0) {
        region_type |= GF_REGION_QUAD_TABLE;
      } else if (strcmp(x, "SINGLE") == 0) {
        region_type |= GF_REGION_SINGLE_TABLE;
      } else if (strcmp(x, "LAZY") == 0) {
        region_type |= GF_REGION_LAZY;
      } else if (strcmp(x, "SSE") == 0) {
        region_type |= GF_REGION_SSE;
      } else if (strcmp(x, "NOSSE") == 0) {
        region_type |= GF_REGION_NOSSE;
      } else if (strcmp(x, "CAUCHY") == 0) {
        region_type |= GF_REGION_CAUCHY;
      } else if (strcmp(x, "ALTMAP") == 0) {
        region_type |= GF_REGION_ALTMAP;
      } else if (strcmp(x, "STDMAP") == 0) {
        region_type |= GF_REGION_STDMAP;
      } else {
        if (base != NULL) gf_free(base, 1);
        free(crt);
        return 0;
      }
      if (y != NULL) x = y+1;
    } while (y != NULL);
    free(crt);
  }

  starting++;

  if (strcmp(argv[starting], "-") == 0) {
    divide_type = GF_DIVIDE_DEFAULT;
  } else if (strcmp(argv[starting], "MATRIX") == 0) {
    divide_type = GF_DIVIDE_MATRIX;
  } else if (strcmp(argv[starting], "EUCLID") == 0) {
    divide_type = GF_DIVIDE_EUCLID;
  } else {
    if (base != NULL) gf_free(base, 1);
    return 0;
  }
  starting++;

  if (!gf_init_hard(gf, w, mult_type, region_type, divide_type, prim_poly, arg1, arg2, base, NULL)) {
    if (base != NULL) gf_free(base, 1);
    return 0;
  }
  return starting;
}

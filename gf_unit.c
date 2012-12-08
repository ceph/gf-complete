/*
 * gf_unit.c
 *
 * Performs unit testing for gf arithmetic
 */

#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "gf.h"
#include "gf_int.h"
#include "gf_method.h"
#include "gf_rand.h"
#include "gf_general.h"

#define REGION_SIZE (16384) 

void problem(char *s)
{
  fprintf(stderr, "Unit test failed.\n");
  fprintf(stderr, "%s\n", s);
  exit(1);
}

void usage(char *s)
{
  fprintf(stderr, "usage: gf_unit w tests seed [method] - does unit testing in GF(2^w)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Legal w are: 1 - 32, 64 and 128\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Tests may be any combination of:\n");
  fprintf(stderr, "       A: All\n");
  fprintf(stderr, "       S: Single operations (multiplication/division)\n");
  fprintf(stderr, "       R: Region operations\n");
  fprintf(stderr, "       V: Verbose Output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Use -1 for time(0) as a seed.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "For method specification, type gf_methods\n");
  fprintf(stderr, "\n");
  if (s != NULL) fprintf(stderr, "%s\n", s);
  exit(1);
}

int main(int argc, char **argv)
{
  int w, i, verbose, single, region, tested, top;
  int start, end, xor;
  gf_t   gf, gf_def;
  time_t t0;
  gf_internal_t *h;
  gf_general_t *a, *b, *c, *d, *ai, *bi;
  char as[50], bs[50], cs[50], ds[50], ais[50], bis[50];
  uint32_t mask;
  char *ra, *rb, *rc, *rd, *target;
  int align;

  if (argc < 4) usage(NULL);
  if (sscanf(argv[1], "%d", &w) == 0) usage("Bad w\n");
  if (sscanf(argv[3], "%ld", &t0) == 0) usage("Bad seed\n");
  if (t0 == -1) t0 = time(0);
  MOA_Seed(t0);

  if (w > 32 && w != 64 && w != 128) usage("Bad w");
  
  if (create_gf_from_argv(&gf, w, argc, argv, 4) == 0) usage("Bad Method");

  for (i = 0; i < strlen(argv[2]); i++) {
    if (strchr("ASRV", argv[2][i]) == NULL) usage("Bad test\n");
  }

  h = (gf_internal_t *) gf.scratch;
  a = (gf_general_t *) malloc(sizeof(gf_general_t));
  b = (gf_general_t *) malloc(sizeof(gf_general_t));
  c = (gf_general_t *) malloc(sizeof(gf_general_t));
  d = (gf_general_t *) malloc(sizeof(gf_general_t));
  ai = (gf_general_t *) malloc(sizeof(gf_general_t));
  bi = (gf_general_t *) malloc(sizeof(gf_general_t));

  ra = (char *) malloc(sizeof(char)*REGION_SIZE);
  rb = (char *) malloc(sizeof(char)*REGION_SIZE);
  rc = (char *) malloc(sizeof(char)*REGION_SIZE);
  rd = (char *) malloc(sizeof(char)*REGION_SIZE);

  if (w <= 32) {
    mask = 0;
    for (i = 0; i < w; i++) mask |= (1 << i);
  }

  verbose = (strchr(argv[2], 'V') != NULL);
  single = (strchr(argv[2], 'S') != NULL || strchr(argv[2], 'A') != NULL);
  region = (strchr(argv[2], 'R') != NULL || strchr(argv[2], 'A') != NULL);

  if (!gf_init_easy(&gf_def, w, GF_MULT_DEFAULT)) problem("No default for this value of w");
  
  if (verbose) printf("Seed: %ld\n", t0);

  if (single) {
    
    if (gf.multiply.w32 == NULL) problem("No multiplication operation defined.");
    if (verbose) { printf("Testing single multiplications/divisions.\n"); fflush(stdout); }
    if (w <= 10) {
      top = (1 << w)*(1 << w);
    } else {
      top = 1024*1024;
    }
    for (i = 0; i < top; i++) {
      if (w <= 10) {
        a->w32 = i % (1 << w);
        b->w32 = (i >> w);
      } else if (i < 10) {
        gf_general_set_zero(a, w);
        gf_general_set_random(b, w, 1);
      } else if (i < 20) {
        gf_general_set_random(a, w, 1);
        gf_general_set_zero(b, w);
      } else if (i < 30) {
        gf_general_set_one(a, w);
        gf_general_set_random(b, w, 1);
      } else if (i < 40) {
        gf_general_set_random(a, w, 1);
        gf_general_set_one(b, w);
      } else {
        gf_general_set_random(a, w, 1);
        gf_general_set_random(b, w, 1);
      }

      tested = 0;
      gf_general_multiply(&gf, a, b, c);
      
      /* If this is not composite, then first test against the default: */

      if (h->mult_type != GF_MULT_COMPOSITE) {
        tested = 1;
        gf_general_multiply(&gf_def, a, b, d);

        if (!gf_general_are_equal(c, d, w)) {
          gf_general_val_to_s(a, w, as);
          gf_general_val_to_s(b, w, bs);
          gf_general_val_to_s(c, w, cs);
          gf_general_val_to_s(d, w, ds);
          printf("Error in single multiplication (all numbers in hex):\n\n");
          printf("  gf.multiply(gf, %s, %s) = %s\n", as, bs, cs);
          printf("  The default gf multiplier returned %s\n", ds);
          exit(1);
        }
      }

      /* Now, we also need to double-check by other means, in case the default is wanky, 
         and when we're performing composite operations. Start with 0 and 1, where we know
         what the result should be. */

      if (gf_general_is_zero(a, w) || gf_general_is_zero(b, w) || 
          gf_general_is_one(a, w)  || gf_general_is_one(b, w)) {
        tested = 1;
        if (((gf_general_is_zero(a, w) || gf_general_is_zero(b, w)) && !gf_general_is_zero(c, w)) ||
            (gf_general_is_one(a, w) && !gf_general_are_equal(b, c, w)) ||
            (gf_general_is_one(b, w) && !gf_general_are_equal(a, c, w))) {
          gf_general_val_to_s(a, w, as);
          gf_general_val_to_s(b, w, bs);
          gf_general_val_to_s(c, w, cs);
          printf("Error in single multiplication (all numbers in hex):\n\n");
          printf("  gf.multiply(gf, %s, %s) = %s, which is clearly wrong.\n", as, bs, cs);
;
          exit(1);
        }
      }

      /* Dumb check to make sure that it's not returning numbers that are too big: */

      if (w < 32 && (c->w32 & mask) != c->w32) {
        gf_general_val_to_s(a, w, as);
        gf_general_val_to_s(b, w, bs);
        gf_general_val_to_s(c, w, cs);
        printf("Error in single multiplication (all numbers in hex):\n\n");
        printf("  gf.multiply.w32(gf, %s, %s) = %s, which is too big.\n", as, bs, cs);
        exit(1);
      }
    }
  }

  if (region) {
    if (verbose) { printf("Testing region multiplications\n"); fflush(stdout); }
    for (i = 0; i < 1000; i++) {
      if (i < 20) {
        gf_general_set_zero(a, w);
      } else if (i < 40) {
        gf_general_set_one(a, w);
      } else {
        gf_general_set_random(a, w, 1);
      }
      MOA_Fill_Random_Region(ra, REGION_SIZE);
      MOA_Fill_Random_Region(rb, REGION_SIZE);
      xor = i%2;
      align = w/8;
      if (align == 0) align = 1;
      if (align > 16) align = 16;
      if ((h->region_type & GF_REGION_CAUCHY) || (w < 32 && w != 4 && w != 8 && w != 16)) {
        start = MOA_Random_W(5, 1);
        end = REGION_SIZE - MOA_Random_W(5, 1);
        target = rb;
        while ((end-start)%w != 0) end--;
      } else {
        start = MOA_Random_W(5, 1) * align;
        end = REGION_SIZE - (MOA_Random_W(5, 1) * align);
        if (h->mult_type == GF_MULT_COMPOSITE && (h->region_type & GF_REGION_ALTMAP)) {
          target = rb ;
        } else {
          target = ((i%4)/2) ? rb : ra;
        }
      }
      memcpy(rc, ra, REGION_SIZE);
      memcpy(rd, target, REGION_SIZE);
      gf_general_do_region_multiply(&gf, a, ra+start, target+start, end-start, xor);
      gf_general_do_region_check(&gf, a, rc+start, rd+start, target+start, end-start, xor);
    }
  }
}

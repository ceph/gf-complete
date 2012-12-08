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

#define REGION_SIZE (65536) 

static
uint8_t get_alt_map_2w8(int offset, uint8_t *buf, int region_size)
{
  uint8_t symbol = 0;
  int bit_off = offset % 2;

  if (bit_off == 0) {
    symbol = buf[offset / 2] & 0x0f | ((buf[(offset / 2)+region_size] & 0x0f) << 4);
  } else {
    symbol = ((buf[offset / 2] & 0xf0) >> 4) | (buf[(offset / 2)+region_size] & 0xf0);
  }
  
  return symbol;
}

static
uint16_t get_alt_map_2w16(int offset, uint8_t *buf, int region_size)
{
  uint16_t symbol = 0;

  symbol = buf[offset] | (buf[offset+region_size] << 8);

  return symbol;
}

static
uint32_t get_alt_map_2w32(int offset, uint8_t *buf, int region_size)
{
  uint32_t symbol = 0;
  uint16_t buf_a = buf[offset] | (buf[offset + 1] << 8);
  uint16_t buf_b = buf[offset + region_size] | (buf[offset + region_size + 1] << 8);
  
  symbol = buf_a | (buf_b << 16);

  return symbol;
}

static
void test_alt_map()
{
  uint8_t* buf = (uint8_t*)malloc(sizeof(uint8_t)*REGION_SIZE);  
  int i=0;
  uint8_t c=1, next_c;

  for (i=0; i < REGION_SIZE/2;i++) {
    if (c == 255) c = 1;
    buf[i] = c;
    buf[i+(REGION_SIZE/2)] = c; 
    c++;
  }


  c = 1;
  for (i=0; i < REGION_SIZE;i++) {
    uint8_t sym_w8 = get_alt_map_2w8(i, buf, REGION_SIZE/2); 
    uint8_t c_val = ((i % 2) == 0) ? (c & 0x0f) : ((c & 0xf0) >> 4);
    uint8_t exp_sym_w8 = c_val | c_val << 4;
    
    if (exp_sym_w8 != sym_w8) {
      fprintf(stderr, "Alt mapping failure (w=8,c=%d,i=%d): %u != %u\n", c, i, exp_sym_w8, sym_w8);
      exit(1);
    }
    
    if ((i % 2) == 1) {
      c++;
    }
    if (c == 255) {
      c = 1;
    }
  }
  
  c = 1;

  for (i=0; i < REGION_SIZE/2;i++) {
    uint16_t sym_w16 = get_alt_map_2w16(i, buf, REGION_SIZE/2); 
    uint16_t exp_sym_w16 = c | c << 8;

    if (exp_sym_w16 != sym_w16) {
      fprintf(stderr, "Alt mapping failure (w=16,c=%d,i=%d): %u != %u\n", c, i, exp_sym_w16, sym_w16);
      exit(1);
    }
    
    c++;
    if (c == 255) {
      c = 1;
    }
  }

  c = 1;
  next_c = 2;

  for (i=0; i < REGION_SIZE/4;i++) {
    uint32_t sym_w32 = get_alt_map_2w32(i, buf, REGION_SIZE/2); 
    uint32_t exp_sym_w32 = c | (next_c << 8) | c << 16 | (next_c << 24);

    if (exp_sym_w32 != sym_w32) {
      fprintf(stderr, "Alt mapping failure (w=32,c=%d,i=%d): %u != %u\n", c, i, exp_sym_w32, sym_w32);
      exit(1);
    }
    c++;
    next_c++;
    if (c == 255) {
      c = 1;
      next_c = 2;
    } else if (c == 254) {
      next_c = 1;
    } 
  }

}

void fill_random_region(void *reg, int size)
{
  uint32_t *r;
  int i;

  r = (uint32_t *) reg;
  for (i = 0; i < size/sizeof(uint32_t); i++) {
    r[i] = MOA_Random_32();
  }
}

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
  fprintf(stderr, "Legal w are: 4, 8, 16, 32, 64 and 128\n");
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
  int w, i, j, verbose, single, region, xor, off, size, sindex, eindex, tested, top;
  uint32_t a, b, c, d, ai, da, bi, mask;
  uint64_t a64, b64, c64, d64;
  uint64_t a128[2], b128[2], c128[2], d128[2], e128[2];
  gf_t                gf, gf_def;
  uint8_t *r8b, *r8c, *r8d;
  uint16_t *r16b, *r16c, *r16d;
  uint32_t *r32b, *r32c, *r32d;
  uint64_t *r64b, *r64c, *r64d;
  uint64_t *r128b, *r128c, *r128d;
  time_t t0;
  gf_internal_t *h;

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
  if (w <= 32) {
    mask = 0;
    for (i = 0; i < w; i++) mask |= (1 << i);
  }

  verbose = (strchr(argv[2], 'V') != NULL);
  single = (strchr(argv[2], 'S') != NULL || strchr(argv[2], 'A') != NULL);
  region = (strchr(argv[2], 'R') != NULL || strchr(argv[2], 'A') != NULL);

  if (((h->region_type & GF_REGION_ALTMAP) != 0) && (h->mult_type == GF_MULT_COMPOSITE)) {
    test_alt_map();
  }

  if (!gf_init_easy(&gf_def, w, GF_MULT_DEFAULT)) problem("No default for this value of w");
  
  if (verbose) printf("Seed: %ld\n", t0);

  if (single) {
    
    if (w <= 32) {
      if (gf.multiply.w32 == NULL) problem("No multiplication operation defined.");
      if (verbose) { printf("Testing single multiplications/divisions.\n"); fflush(stdout); }
      if (w <= 10) {
        top = (1 << w)*(1 << w);
      } else {
        top = 1000000;
      }
      for (i = 0; i < top; i++) {
        if (w <= 10) {
          a = i % (1 << w);
          b = i >> w;
        } else if (i < 10) {
          a = 0;
          b = MOA_Random_W(w, 1);
        } else if (i < 20) {
          b = 0;
          a = MOA_Random_W(w, 1);
        } else if (i < 30) {
          a = 1;
          b = MOA_Random_W(w, 1);
        } else if (i < 40) {
          b = 1;
          a = MOA_Random_W(w, 1);
        } else {
          a = MOA_Random_W(w, 1);
          b = MOA_Random_W(w, 1);
        }

        c = gf.multiply.w32(&gf, a, b);
        tested = 0;

        /* If this is not composite, then first test against the default: */

        if (h->mult_type != GF_MULT_COMPOSITE) {
          tested = 1;
          d = gf_def.multiply.w32(&gf_def, a, b);

          if (c != d) {
            printf("Error in single multiplication (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
            printf("  The default returned %x\n", d);
            exit(1);
          }
        }

        /* Now, we also need to double-check, in case the default is wanky, and when
           we're performing composite operations. Start with 0 and 1: */

        if (a == 0 || b == 0 || a == 1 || b == 1) {
          tested = 1;
          if (((a == 0 || b == 0) && c != 0) ||
              (a == 1 && c != b) ||
              (b == 1 && c != a)) {
            printf("Error in single multiplication (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x, which is clearly wrong.\n", a, b, c);
            exit(1);
          }
        
        /* If division or inverses are defined, let's test all combinations to make sure
           that the operations are consistent with each other. */
            
        } else {
          if ((c & mask) != c) {
            printf("Error in single multiplication (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x, which is too big.\n", a, b, c);
            exit(1);
          }
            
        }
        if (gf.inverse.w32 != NULL && (a != 0 || b != 0)) {
          tested = 1;
          if (a != 0) {
            ai = gf.inverse.w32(&gf, a);

            if (gf.multiply.w32(&gf, c, ai) != b) {
              printf("Error in single multiplication (all numbers in hex):\n\n");
              printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
              printf("  The inverse of %x is %x, and gf_multiply.w32() of %x and %x equals %x\n",
                        a, ai, c, ai, gf.multiply.w32(&gf, c, ai));
              exit(1);
            }
          }
          if (b != 0) {
            bi = gf.inverse.w32(&gf, b);
            if (gf.multiply.w32(&gf, c, bi) != a) {
              printf("Error in single multiplication (all numbers in hex):\n\n");
              printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
              printf("  The inverse of %x is %x, and gf_multiply.w32() of %x and %x equals %x\n",
                        b, bi, c, bi, gf.multiply.w32(&gf, c, bi));
              exit(1);
            }
          }
        }
        if (gf.divide.w32 != NULL && (a != 0 || b != 0)) {
          tested = 1;
        
          if (a != 0) {
            ai = gf.divide.w32(&gf, c, a);

            if (ai != b) {
              printf("Error in single multiplication/division (all numbers in hex):\n\n");
              printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
              printf("  gf.divide.w32() of %x and %x returned %x\n", c, a, ai);
              exit(1);
            }
          }
          if (b != 0) {
            bi = gf.divide.w32(&gf, c, b);

            if (bi != a) {
              printf("Error in single multiplication/division (all numbers in hex):\n\n");
              printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
              printf("  gf.divide.w32() of %x and %x returned %x\n", c, b, bi);
              exit(1);
            }
          }
        }

        if (!tested) problem("There is no way to test multiplication.\n");
      }

    } else if (w == 64) {
      if (verbose) { printf("Testing single multiplications/divisions.\n"); fflush(stdout); }
      if (gf.multiply.w64 == NULL) problem("No multiplication operation defined.");
      for (i = 0; i < 1000; i++) { 
        for (j = 0; j < 1000; j++) { 
          a64 = MOA_Random_64();
          b64 = MOA_Random_64();
          c64 = gf.multiply.w64(&gf, a64, b64);
          if ((a64 == 0 || b64 == 0) && c64 != 0) problem("Single Multiplication by zero Failed");
          if (a64 != 0 && b64 != 0) {
            d64 = (gf.divide.w64 == NULL) ? gf_def.divide.w64(&gf_def, c64, b64) : gf.divide.w64(&gf, c64, b64);
            if (d64 != a64) {
              printf("0x%llx * 0x%llx =? 0x%llx (check-a: 0x%llx)\n", a64, b64, c64, d64);
              problem("Single multiplication/division failed");
            }
          }
        }
      }
      if (gf.inverse.w64 == NULL) {
        printf("No inverse defined for this method.\n");
      } else {
        if (verbose) { printf("Testing Inversions.\n"); fflush(stdout); }
        for (i = 0; i < 1000; i++) {
          do { a64 = MOA_Random_64(); } while (a64 == 0);
          b64 = gf.inverse.w64(&gf, a64);
          if (gf.multiply.w64(&gf, a64, b64) != 1) problem("Inversion failed.\n"); 
        }
      }
    } else if (w == 128) {
      if (verbose) { printf("Testing single multiplications/divisions.\n"); fflush(stdout); }
      if (gf.multiply.w128 == NULL) problem("No multiplication operation defined.");
      for (i = 0; i < 500; i++) { 
        for (j = 0; j < 500; j++) { 
          MOA_Random_128(a128);
          MOA_Random_128(b128);
          gf.multiply.w128(&gf, a128, b128, c128);
          if ((GF_W128_IS_ZERO(a128) && GF_W128_IS_ZERO(b128)) && !(GF_W128_IS_ZERO(c128))) problem("Single Multiplication by zero Failed");
          if (!GF_W128_IS_ZERO(a128) && !GF_W128_IS_ZERO(b128)) {
            gf.divide.w128 == NULL ? gf_def.divide.w128(&gf_def, c128, b128, d128) : gf.divide.w128(&gf, c128, b128, d128);
            if (!GF_W128_EQUAL(a128, d128)) {
              printf("0x%llx 0x%llx * 0x%llx 0x%llx =? 0x%llx 0x%llx (check-a: 0x%llx 0x%llx)\n", a128[0], a128[1], b128[0], b128[1], c128[0], c128[1], d128[0], d128[1]);
              problem("Single multiplication/division failed");
            }
          }
        }
      }
      if (gf.inverse.w128 == NULL) {
        printf("No inverse defined for this method.\n");
      } else {
        if (verbose) { printf("Testing Inversions.\n"); fflush(stdout); }
        for (i = 0; i < 1000; i++) {
          do { MOA_Random_128(a128); } while (GF_W128_IS_ZERO(a128));
          gf.inverse.w128(&gf, a128, b128);
          gf.multiply.w128(&gf, a128, b128, c128);
          if (!(c128[0] == 0 && c128[1] == 1)) problem("Inversion failed.\n"); 
        }
      }

    } else {
      problem("Value of w not implemented yet");
    }
  }

  if (region) {
    
    if (w == 4) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r8b = (uint8_t *) malloc(REGION_SIZE);
        r8c = (uint8_t *) malloc(REGION_SIZE);
        r8d = (uint8_t *) malloc(REGION_SIZE);
        fill_random_region(r8b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 16; a++) {
            fill_random_region(r8c, REGION_SIZE);
            memcpy(r8d, r8c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8c+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = (r8b[i] >> 4);
              c = (r8c[i] >> 4);
              d = (r8d[i] >> 4);
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d * %d = %d, but should equal %d\n", a, b, c, gf.multiply.w32(&gf, a, b) );
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d %d %d %d\n", a, b, c, d);
                printf("       %d %d %d %d\n", a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
              b = (r8b[i] & 0xf);
              c = (r8c[i] & 0xf);
              d = (r8d[i] & 0xf);
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d * %d = %d, but should equal %d\n", a, b, c, gf.multiply.w32(&gf, a, b) );
                printf("i=%d.  0x%x 0x%x 0x%x 0x%x\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       (%d * %d ^ %d) should equal %d - equals %d\n", 
                    a, b, d, (gf.multiply.w32(&gf, a, b) ^ d), c);
                printf("       %d %d %d %d\n", a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 16; a++) {
            fill_random_region(r8b, REGION_SIZE);
            memcpy(r8d, r8b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8b+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = (r8b[i] >> 4);
              d = (r8d[i] >> 4);
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("i=%d.  %d %d %d\n", i, a, b, d);
                printf("i=%d.  %d %d %d\n", i, a, r8b[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
              b = (r8b[i] & 0xf);
              d = (r8d[i] & 0xf);
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("%d %d %d\n", a, b, d);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r8b);
        free(r8c);
        free(r8d);
      }
    } else if (w == 8) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r8b = (uint8_t *) malloc(REGION_SIZE);
        r8c = (uint8_t *) malloc(REGION_SIZE);
        r8d = (uint8_t *) malloc(REGION_SIZE);
        fill_random_region(r8b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 256; a++) {
            fill_random_region(r8c, REGION_SIZE);
            memcpy(r8d, r8c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE;
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8c+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w8(i, (uint8_t*)r8b, REGION_SIZE / 2);
                c = get_alt_map_2w8(i, (uint8_t*)r8c, REGION_SIZE / 2);
                d = get_alt_map_2w8(i, (uint8_t*)r8d, REGION_SIZE / 2);
              } else {
                b = r8b[i];
                c = r8c[i];
                d = r8d[i];
              }
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  %d %d %d %d\n", i, a, b, c, d);
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                printf("%llx.  Sindex: %d\n", r8b+i, sindex);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  %d %d %d %d\n", i, a, b, c, d);
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;  
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 256; a++) {
            fill_random_region(r8b, REGION_SIZE);
            memcpy(r8d, r8b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8b+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = r8b[i];
              d = r8d[i];
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("i=%d.  %d %d %d\n", i, a, b, d);
                printf("i=%d.  %d %d %d\n", i, a, r8b[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r8b);
        free(r8c);
        free(r8d);
      }
    } else if (w == 16) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r16b = (uint16_t *) malloc(REGION_SIZE);
        r16c = (uint16_t *) malloc(REGION_SIZE);
        r16d = (uint16_t *) malloc(REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            fill_random_region(r16b, REGION_SIZE);
            a = MOA_Random_W(w, 0);
            fill_random_region(r16c, REGION_SIZE);
            memcpy(r16d, r16c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE / sizeof(uint16_t);
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint16_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint16_t);
            gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16c+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r16c+sindex), (void *) (r16d+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r16c+sindex), (void *) (r16d+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r16d+sindex), (void *) (r16b+sindex), ai, size, xor);
            }
  
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w16(i, (uint8_t*)r16b, size / 2);
                c = get_alt_map_2w16(i, (uint8_t*)r16c, size / 2);
                d = get_alt_map_2w16(i, (uint8_t*)r16d, size / 2);
              } else {
                b = r16b[i];
                c = r16c[i];
                d = r16d[i];
              }
              if (!xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We have %d * %d = %d, and %d * %d = %d.\n", a, b, c, c, ai, d);
                printf("%d is the inverse of %d\n", ai, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && b != 0) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("   b should equal 0, but it doesn't.  Probe into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_W(w, 0);
            fill_random_region(r16b, REGION_SIZE);
            memcpy(r16d, r16b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint16_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint16_t);
            gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r16d+sindex), (void *) (r16b+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), ai, size, 0);
            }
  
            for (i = sindex; i < eindex; i++) {
              b = r16b[i];
              c = r16c[i];
              d = r16d[i];
              if (!xor && (d != b)) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=b; b = ba; b = b(a^-1).\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                printf("b = %d.  d = %d.  a = %d\n", b, d, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=b; b = b + ba; b += d; b = b(a^-1);\n");
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r16b);
        free(r16c);
        free(r16d);
      }
    } else if (w == 32) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r32b = (uint32_t *) malloc(REGION_SIZE);
        r32c = (uint32_t *) malloc(REGION_SIZE);
        r32d = (uint32_t *) malloc(REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_32();
            fill_random_region(r32b, REGION_SIZE);
            fill_random_region(r32c, REGION_SIZE);
            memcpy(r32d, r32c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE / sizeof(uint32_t);
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint32_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint32_t);
            gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32c+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r32c+sindex), (void *) (r32d+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r32c+sindex), (void *) (r32d+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r32d+sindex), (void *) (r32b+sindex), ai, size, xor);
            }
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w32(i, (uint8_t*)r32b, size / 2);
                c = get_alt_map_2w32(i, (uint8_t*)r32c, size / 2);
                d = get_alt_map_2w32(i, (uint8_t*)r32d, size / 2);
                i++;
              } else {
                b = r32b[i];
                c = r32c[i];
                d = r32d[i];
              }
              if (!xor && d != b) {
                printf("i=%d.  Addresses: b: 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We have %d * %d = %d, and %d * %d = %d.\n", a, b, c, c, ai, d);
                printf("%d is the inverse of %d\n", ai, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && b != 0) {
                printf("i=%d.  Addresses: b: 0x%lx   c: 0x%lx   d: 0x%lx\n", i, 
                      (unsigned long) (r32b+i), (unsigned long) (r32c+i), (unsigned long) (r32d+i));
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("   b should equal 0, but it doesn't.  Probe into it.\n");
                printf("a: %8x  b: %8x  c: %8x,  d: %8x\n", a, b, c, d);
                problem("Failed buffer-constant, xor=1");
              }

            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_32();
            fill_random_region(r32b, REGION_SIZE);
            memcpy(r32d, r32b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint32_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint32_t);
            gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r32d+sindex), (void *) (r32b+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), ai, size, 0);
            }

            for (i = sindex; i < eindex; i++) {
              b = r32b[i];
              c = r32c[i];
              d = r32d[i];
              if (!xor && (d != b)) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We did d=b; b = ba; b = b(a^-1).\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                printf("b = %d.  d = %d.  a = %d\n", b, d, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We did d=b; b = b + ba; b += d; b = b(a^-1);\n");
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r32b);
        free(r32c);
        free(r32d);
      }
    } else if (w == 64) {
      if (gf.multiply_region.w64 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r64b = (uint64_t *) malloc(REGION_SIZE);
        r64c = (uint64_t *) malloc(REGION_SIZE);
        r64d = (uint64_t *) malloc(REGION_SIZE);
        fill_random_region(r64b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a64 = MOA_Random_64();
            fill_random_region(r64c, REGION_SIZE);
            memcpy(r64d, r64c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint64_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t);
            gf.multiply_region.w64(&gf, (void *) (r64b+sindex), (void *) (r64c+sindex), a64, size, xor);
            for (i = sindex; i < eindex; i++) {
              b64 = r64b[i];
              c64 = r64c[i];
              d64 = r64d[i];
              if (!xor && gf.multiply.w64(&gf, a64, b64) != c64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx should be 0x%llx\n", i, a64, b64, c64, 
                            gf.multiply.w64(&gf, a64, b64));
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64c[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w64(&gf, a64, b64) ^ d64) != c64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx 0x%llx\n", i, a64, b64, c64, d64);
                printf("i=%d.  0x%llx 0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64c[i], r64d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a64 = MOA_Random_64();
            fill_random_region(r64b, REGION_SIZE);
            memcpy(r64d, r64b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint64_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t);
            gf.multiply_region.w64(&gf, (void *) (r64b+sindex), (void *) (r64b+sindex), a64, size, xor);
            for (i = sindex; i < eindex; i++) {
              b64 = r64b[i];
              d64 = r64d[i];
              if (!xor && gf.multiply.w64(&gf, a64, d64) != b64) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w64(&gf, a64, d64) ^ d64) != b64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, b64, d64);
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r64b);
        free(r64c);
        free(r64d);
      }
    } else if (w == 128) {
      if (gf.multiply_region.w128 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r128b = (uint64_t *) malloc(REGION_SIZE);
        r128c = (uint64_t *) malloc(REGION_SIZE);
        r128d = (uint64_t *) malloc(REGION_SIZE);
        fill_random_region(r128b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            MOA_Random_128(a128);
            fill_random_region(r128c, REGION_SIZE);
            memcpy(r128d, r128c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/(2*sizeof(uint64_t))-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t)*2;
            gf.multiply_region.w128(&gf, (void *) (r128b+sindex*2), (void *) (r128c+sindex*2), a128, size, xor);
            for (i = sindex; i < eindex; i++) {
              b128[0] = r128b[2*i];
              b128[1] = r128b[2*i+1];
              c128[0] = r128c[2*i];
              c128[1] = r128c[2*i+1];
              d128[0] = r128d[2*i];
              d128[1] = r128d[2*i+1];
              gf.multiply.w128(&gf, a128, b128, e128);
              if (xor) {
                e128[0] ^= d128[0];
                e128[1] ^= d128[1];
              }
              if (!xor && !GF_W128_EQUAL(c128, e128)) {
                printf("i=%d.  0x%llx%llx 0x%llx%llx 0x%llx%llx should be 0x%llx%llx\n",
                        i, a128[0], a128[1], b128[0], b128[1], c128[0], c128[1], e128[0], e128[1]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && !GF_W128_EQUAL(e128, c128)) {
                printf("i=%d.  0x%llx%llx 0x%llx%llx 0x%llx%llx 0x%llx%llx\n", i,
                        a128[0], a128[1], b128[0], b128[1], c128[0], c128[1], d128[0], d128[1]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            MOA_Random_128(a128);
            fill_random_region(r128b, REGION_SIZE);
            memcpy(r128d, r128b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            sindex = 0;
            eindex = REGION_SIZE/(2*sizeof(uint64_t))-MOA_Random_W(3, 1);
            eindex = REGION_SIZE/(2*sizeof(uint64_t));
            size = (eindex-sindex)*sizeof(uint64_t)*2;
            gf.multiply_region.w128(&gf, (void *) (r128b+sindex), (void *) (r128b+sindex), a128, size, xor);
            for (i = sindex; i < eindex; i++) {
              b128[0] = r128b[2*i];
              b128[1] = r128b[2*i + 1];
              d128[0] = r128d[2*i];
              d128[1] = r128d[2*i + 1];
              gf.multiply.w128(&gf, a128, d128, e128);
              if (xor) {
                e128[0] ^= d128[0];
                e128[1] ^= d128[1];
              }
              if (!xor && !GF_W128_EQUAL(b128, e128)) problem("Failed buffer-constant, xor=0");
              if (xor && !GF_W128_EQUAL(b128, e128)) {
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r128b);
        free(r128c);
        free(r128d);
      }
    }
  }
  exit(0);
}

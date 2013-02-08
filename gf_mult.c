/*
 * gf_mult.c
 *
 * Multiplies two numbers in gf_2^w
 */

#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "gf_complete.h"
#include "gf_method.h"

void usage(char *s)
{
  fprintf(stderr, "usage: gf_mult a b w [method] - does multiplication of a and b in GF(2^w)\n");
  fprintf(stderr, "       If w has an h on the end, treat a, b and the product as hexadecimal (no 0x)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       legal w are: 1-32, 64 and 128\n");
  fprintf(stderr, "           128 is hex only (i.e. '128' will be an error - do '128h')\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       For method specification, type gf_methods\n");

  if (s != NULL) fprintf(stderr, "%s", s);
  exit(1);
}

int read_128(char *s, uint64_t *v)
{
  int l, t;
  char save;

  l = strlen(s);
  if (l > 32) return 0;

  if (l > 16) {
    if (sscanf(s + (l-16), "%llx", (long long unsigned int *) &(v[1])) == 0) return 0;
    save = s[l-16];
    s[l-16] = '\0';
    t = sscanf(s, "%llx", (long long unsigned int *) &(v[0]));
    s[l-16] = save;
    return t;
  } else {
    v[0] = 0;
    return sscanf(s, "%llx", (long long unsigned int *)&(v[1]));
  }
  return 1;
}

void print_128(uint64_t *v) 
{
  if (v[0] > 0) {
    printf("%llx", (long long unsigned int) v[0]);
    printf("%016llx", (long long unsigned int) v[1]);
  } else {
    printf("%llx", (long long unsigned int) v[1]);
  }
  printf("\n");
}


int main(int argc, char **argv)
{
  int hex, al, bl, w;
  uint32_t a, b, c, top;
  uint64_t a64, b64, c64;
  uint64_t a128[2], b128[2], c128[2];
  char *format;
  gf_t gf;

  if (argc < 4) usage(NULL);
  if (sscanf(argv[3], "%d", &w) == 0) usage("Bad w\n");

  if (w <= 0 || (w > 32 && w != 64 && w != 128)) usage("Bad w");

  hex = (strchr(argv[3], 'h') != NULL);
  if (create_gf_from_argv(&gf, w, argc, argv, 4) == 0) usage("\nBad Method\n");

  if (!hex && w == 128) usage(NULL);
 
  if (w <= 32) {
    format = (hex) ? "%x" : "%u";
    if (sscanf(argv[1], format, &a) == 0) usage("Bad a\n");
    if (sscanf(argv[2], format, &b) == 0) usage("Bad b\n");

    if (w < 32) {
      top = (w == 31) ? 0x80000000 : (1 << w);
      if (w != 32 && a >= top) usage("a is too large\n");
      if (w != 32 && b >= top) usage("b is too large\n");
    }
  
    c = gf.multiply.w32(&gf, a, b);
    printf(format, c);
    printf("\n");

  } else if (w == 64) {
    format = (hex) ? "%llx" : "%llu";
    if (sscanf(argv[1], format, &a64) == 0) usage("Bad a\n");
    if (sscanf(argv[2], format, &b64) == 0) usage("Bad b\n");
    c64 = gf.multiply.w64(&gf, a64, b64);

    printf(format, c64);
    printf("\n");

  } else if (w == 128) {

    if (read_128(argv[1], a128) == 0) usage("Bad a\n");
    if (read_128(argv[2], b128) == 0) usage("Bad b\n");
    gf.multiply.w128(&gf, a128, b128, c128);

    print_128(c128);
  }
  exit(0);
}

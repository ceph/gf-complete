/*
 * Multiplies four and five in GF(2^4).
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "gf_complete.h"

main()
{
  gf_t gf;
  void *scratch;
  int size;

  size = gf_scratch_size(16, GF_MULT_SPLIT_TABLE,
                             GF_REGION_SSE | GF_REGION_ALTMAP,
                             GF_DIVIDE_DEFAULT,
                             16, 4);
  if (size == -1) exit(1); /* It failed. That shouldn't happen*/
  scratch = (void *) malloc(size);
  if (scratch == NULL) { perror("malloc"); exit(1); }
  if (!gf_init_hard(&gf, 16, GF_MULT_SPLIT_TABLE,
                             GF_REGION_SSE | GF_REGION_ALTMAP,
                             GF_DIVIDE_DEFAULT,
                             0, 16, 4, NULL, scratch)) exit(1);
  printf("Yo\n");
}

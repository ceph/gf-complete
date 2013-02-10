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

  gf_init_easy(&gf, 4);
  printf("%d\n", gf.multiply.w32(&gf, 5, 4));
  exit(0);
}

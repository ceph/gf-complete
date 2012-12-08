/*
 * gf_w64.c
 *
 * Routines for 64-bit Galois fields
 */

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

#define GF_FIELD_WIDTH (64)

static
inline
gf_val_64_t gf_w64_inverse_from_divide (gf_t *gf, gf_val_64_t a)
{
  return gf->divide.w64(gf, 1, a);
}

static
inline
gf_val_64_t gf_w64_divide_from_inverse (gf_t *gf, gf_val_64_t a, gf_val_64_t b)
{
  b = gf->inverse.w64(gf, b);
  return gf->multiply.w64(gf, a, b);
}

static
void
gf_w64_multiply_region_from_single(gf_t *gf, void *src, void *dest, gf_val_64_t val, int bytes, int
xor)
{
  int i;
  gf_val_64_t *s64;
  gf_val_64_t *d64;

  s64 = (gf_val_64_t *) src;
  d64 = (gf_val_64_t *) dest;

  if (xor) {
    for (i = 0; i < bytes/sizeof(gf_val_64_t); i++) {
      d64[i] ^= gf->multiply.w64(gf, val, s64[i]);
    }
  } else {
    for (i = 0; i < bytes/sizeof(gf_val_64_t); i++) {
      d64[i] = gf->multiply.w64(gf, val, s64[i]);
    }
  }
}

static
inline
gf_val_64_t gf_w64_euclid (gf_t *gf, gf_val_64_t b)
{
  gf_val_64_t e_i, e_im1, e_ip1;
  gf_val_64_t d_i, d_im1, d_ip1;
  gf_val_64_t y_i, y_im1, y_ip1;
  gf_val_64_t c_i;
  gf_val_64_t one = 1;

  if (b == 0) return -1;
  e_im1 = ((gf_internal_t *) (gf->scratch))->prim_poly;
  e_i = b;
  d_im1 = 64;
  for (d_i = d_im1-1; ((one << d_i) & e_i) == 0; d_i--) ;
  y_i = 1;
  y_im1 = 0;

  while (e_i != 1) {

    e_ip1 = e_im1;
    d_ip1 = d_im1;
    c_i = 0;

    while (d_ip1 >= d_i) {
      c_i ^= (one << (d_ip1 - d_i));
      e_ip1 ^= (e_i << (d_ip1 - d_i));
      d_ip1--;
      while ((e_ip1 & (one << d_ip1)) == 0) d_ip1--;
    }

    y_ip1 = y_im1 ^ gf->multiply.w64(gf, c_i, y_i);
    y_im1 = y_i;
    y_i = y_ip1;

    e_im1 = e_i;
    d_im1 = d_i;
    e_i = e_ip1;
    d_i = d_ip1;
  }

  return y_i;
}

/* JSP: GF_MULT_SHIFT: The world's dumbest multiplication algorithm.  I only
   include it for completeness.  It does have the feature that it requires no
   extra memory.  
*/

static
inline
gf_val_64_t
gf_w64_shift_multiply (gf_t *gf, gf_val_64_t a64, gf_val_64_t b64)
{
  uint64_t pl, pr, ppl, ppr, i, pp, a, bl, br, one, lbit;
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  ppr = h->prim_poly;
  ppl = 1;
  
  a = a64;
  bl = 0;
  br = b64;
  one = 1;
  lbit = (one << 63);

  pl = 0;
  pr = 0;

  for (i = 0; i < GF_FIELD_WIDTH; i++) {
    if (a & (one << i)) {
      pl ^= bl;
      pr ^= br;
    }
    /* printf("P: %016llx %016llx     ", pl, pr); printf("B: %016llx %016llx\n", bl, br); */
    bl <<= 1;
    if (br & lbit) bl ^= 1;
    br <<= 1;
  }

  one = lbit;
  ppl = ((h->prim_poly >> 1) | lbit);
  ppr = lbit;
  while (one != 0) {
    if (pl & one) {
      pl ^= ppl;
      pr ^= ppr;
    }
    one >>= 1;
    ppr >>= 1;
    if (ppl & 1) ppr ^= lbit;
    ppl >>= 1;
  }
  return pr;
}

static 
int gf_w64_shift_init(gf_t *gf)
{
  gf->multiply.w64 = gf_w64_shift_multiply;
  gf->inverse.w64 = gf_w64_euclid;
  gf->multiply_region.w64 = gf_w64_multiply_region_from_single;
  return 1;
}

int gf_w64_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2)
{
  if (divide_type == GF_DIVIDE_MATRIX) return -1;
  switch(mult_type)
  {
    case GF_MULT_DEFAULT:
    case GF_MULT_SHIFT:
      if (arg1 != 0 || arg2 != 0 || region_type != 0) return -1;
      return sizeof(gf_internal_t);
      break;
    default:
      return -1;
   }
}

int gf_w64_init(gf_t *gf)
{
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  if (h->prim_poly == 0) h->prim_poly = 0x1b; /* Omitting the leftmost 1 as in w=32 */

  gf->multiply.w64 = NULL;
  gf->divide.w64 = NULL;
  gf->inverse.w64 = NULL;
  gf->multiply_region.w64 = NULL;

  switch(h->mult_type) {
    case GF_MULT_DEFAULT: 
    case GF_MULT_SHIFT:     if (gf_w64_shift_init(gf) == 0) return 0; break;
    default: return 0;
  }
  if (h->divide_type == GF_DIVIDE_EUCLID) {
    gf->divide.w64 = gf_w64_divide_from_inverse;
    gf->inverse.w64 = gf_w64_euclid;
  } 

/* else if (h->divide_type == GF_DIVIDE_MATRIX) {
    gf->divide.w64 = gf_w64_divide_from_inverse;
    gf->inverse.w64 = gf_w64_matrix;
  } */

  if (gf->inverse.w64 != NULL && gf->divide.w64 == NULL) {
    gf->divide.w64 = gf_w64_divide_from_inverse;
  }
  if (gf->inverse.w64 == NULL && gf->divide.w64 != NULL) {
    gf->inverse.w64 = gf_w64_inverse_from_divide;
  }
  return 1;
}

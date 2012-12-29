/*
 * gf_w64.c
 *
 * Routines for 64-bit Galois fields
 */

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

#define GF_FIELD_WIDTH (64)
#define GF_FIRST_BIT (1 << 63)

#define GF_BASE_FIELD_WIDTH (32)
#define GF_BASE_FIELD_SIZE       (1 << GF_BASE_FIELD_WIDTH)
#define GF_BASE_FIELD_GROUP_SIZE  GF_BASE_FIELD_SIZE-1

// 10000587 is a valid s for 2^16^2
#define GF_S_GF_16_2_2 (1000587)

// 1000012 is a valid s for 2^32
#define GF_S_GF_32_2 (1000012)

typedef struct w64_composite_int_s {
  uint64_t s; // 's' will be different depending on the base field
} w64_composite_int_t;

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

static
gf_val_64_t gf_w64_extract_word(gf_t *gf, void *start, int bytes, int index)
{
  uint64_t *r64, rv;

  r64 = (uint64_t *) start;
  rv = r64[index];
  return rv;
}

static
gf_val_64_t gf_w64_composite_extract_word(gf_t *gf, void *start, int bytes, int index)
{
  int sub_size;
  gf_internal_t *h;
  uint8_t *r8, *top;
  uint64_t a, b, *r64;
  gf_region_data rd;

  h = (gf_internal_t *) gf->scratch;
  gf_set_region_data(&rd, gf, start, start, bytes, 0, 0, 32);
  r64 = (uint64_t *) start;
  if (r64 + index < (uint64_t *) rd.d_start) return r64[index];
  if (r64 + index >= (uint64_t *) rd.d_top) return r64[index];
  index -= (((uint64_t *) rd.d_start) - r64);
  r8 = (uint8_t *) rd.d_start;
  top = (uint8_t *) rd.d_top;
  sub_size = (top-r8)/2;

  a = h->base_gf->extract_word.w32(h->base_gf, r8, sub_size, index);
  b = h->base_gf->extract_word.w32(h->base_gf, r8+sub_size, sub_size, index);
  return (a | ((uint64_t)b << 32));
}

static
gf_val_64_t
gf_w64_composite_multiply(gf_t *gf, gf_val_64_t a, gf_val_64_t b)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint32_t b0 = b & 0x00000000ffffffff;
  uint32_t b1 = (b & 0xffffffff00000000) >> 32;
  uint32_t a0 = a & 0x00000000ffffffff;
  uint32_t a1 = (a & 0xffffffff00000000) >> 32;
  uint32_t a1b1;
  w64_composite_int_t *comp_int = (w64_composite_int_t*)h->private;

  a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

  return ((uint64_t)(base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) | 
         ((uint64_t)(base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, comp_int->s)) << 32));
}

/*
 * Composite field division trick (explained in 2007 tech report)
 *
 * Compute a / b = a*b^-1, where p(x) = x^2 + sx + 1
 *
 * let c = b^-1
 *
 * c*b = (s*b1c1+b1c0+b0c1)x+(b1c1+b0c0)
 *
 * want (s*b1c1+b1c0+b0c1) = 0 and (b1c1+b0c0) = 1
 *
 * let d = b1c1 and d+1 = b0c0
 *
 * solve s*b1c1+b1c0+b0c1 = 0
 *
 * solution: d = (b1b0^-1)(b1b0^-1+b0b1^-1+s)^-1
 *
 * c0 = (d+1)b0^-1
 * c1 = d*b1^-1
 *
 * a / b = a * c
 */
static
gf_val_64_t
gf_w64_composite_inverse(gf_t *gf, gf_val_64_t a)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint32_t a0 = a & 0x00000000ffffffff;
  uint32_t a1 = (a & 0xffffffff00000000) >> 32;
  uint32_t c0, c1, d, tmp;
  uint64_t c;
  uint32_t a0inv, a1inv;
  w64_composite_int_t *comp_int = (w64_composite_int_t*)h->private;

  if (a0 == 0) {
    a1inv = base_gf->inverse.w32(base_gf, a1);
    c0 = base_gf->multiply.w32(base_gf, a1inv, comp_int->s);
    c1 = a1inv;
  } else if (a1 == 0) {
    c0 = base_gf->inverse.w32(base_gf, a0);
    c1 = 0;
  } else {
    a1inv = base_gf->inverse.w32(base_gf, a1);
    a0inv = base_gf->inverse.w32(base_gf, a0);

    d = base_gf->multiply.w32(base_gf, a1, a0inv);

    tmp = (base_gf->multiply.w32(base_gf, a1, a0inv) ^ base_gf->multiply.w32(base_gf, a0, a1inv) ^ comp_int->s);
    tmp = base_gf->inverse.w32(base_gf, tmp);

    d = base_gf->multiply.w32(base_gf, d, tmp);

    c0 = base_gf->multiply.w32(base_gf, (d^1), a0inv);
    c1 = base_gf->multiply.w32(base_gf, d, a1inv);
  }

  c = c0 | ((uint64_t)c1 << 32);

  return c;
}

static
gf_val_64_t
gf_w64_composite_divide(gf_t *gf, gf_val_64_t a, gf_val_64_t b)
{
  gf_val_64_t binv;

  binv = gf_w64_composite_inverse(gf, b);

  return gf_w64_composite_multiply(gf, a, binv);
}

static
void
gf_w64_composite_multiply_region(gf_t *gf, void *src, void *dest, gf_val_64_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  int i=0;
  uint32_t b0 = val & 0x00000000ffffffff;
  uint32_t b1 = (val & 0xffffffff00000000) >> 32;
  uint64_t *s64, *d64;
  uint64_t *top;
  uint64_t a0, a1, a1b1;
  int num_syms = bytes / 8;
  int sym_divisible = bytes % 4;
  gf_region_data rd;
  w64_composite_int_t *comp_int = (w64_composite_int_t*)h->private;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 8);

  s64 = rd.s_start;
  d64 = rd.d_start;
  top = rd.d_top;
  
  if (xor) {
    while (d64 < top) {
      a0 = *s64 & 0x00000000ffffffff;
      a1 = (*s64 & 0xffffffff00000000) >> 32;
      a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

      *d64 ^= ((uint64_t)(base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) |
                ((uint64_t)(base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, comp_int->s)) << 32));
      s64++;
      d64++;
    }
  } else {
    while (d64 < top) {
      a0 = *s64 & 0x00000000ffffffff;
      a1 = (*s64 & 0xffffffff00000000) >> 32;
      a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

      *d64 = ((base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) |
                ((uint64_t)(base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, comp_int->s)) << 32));
      s64++;
      d64++;
    }
  }
}

static
void
gf_w64_composite_multiply_region_alt(gf_t *gf, void *src, void *dest, gf_val_64_t val, int bytes, int xor)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_val_32_t val0 = val & 0x00000000ffffffff;
  gf_val_32_t val1 = (val & 0xffffffff00000000) >> 32;
  uint8_t *slow, *shigh;
  uint8_t *dlow, *dhigh, *top;
  int sub_reg_size;
  gf_region_data rd;
  w64_composite_int_t *comp_int = (w64_composite_int_t*)h->private;

  if (!xor) {
    memset(dest, 0, bytes);
  }
  
  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 32);
  gf_do_initial_region_alignment(&rd);

  slow = (uint8_t *) rd.s_start;
  dlow = (uint8_t *) rd.d_start;
  top = (uint8_t*) rd.d_top;
  sub_reg_size = (top - dlow)/2;
  shigh = slow + sub_reg_size;
  dhigh = dlow + sub_reg_size;

  base_gf->multiply_region.w32(base_gf, slow, dlow, val0, sub_reg_size, xor);
  base_gf->multiply_region.w32(base_gf, shigh, dlow, val1, sub_reg_size, 1);
  base_gf->multiply_region.w32(base_gf, slow, dhigh, val1, sub_reg_size, xor);
  base_gf->multiply_region.w32(base_gf, shigh, dhigh, val0, sub_reg_size, 1);
  base_gf->multiply_region.w32(base_gf, shigh, dhigh, base_gf->multiply.w32(base_gf, comp_int->s, val1), sub_reg_size, 1);

  gf_do_final_region_alignment(&rd);
}



static
int gf_w64_composite_init(gf_t *gf)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;

  if (h->region_type & GF_REGION_ALTMAP) {
    gf->multiply_region.w64 = gf_w64_composite_multiply_region_alt;
  } else {
    gf->multiply_region.w64 = gf_w64_composite_multiply_region;
  }

  if (h->base_gf != NULL) {
    gf_internal_t *base_h = (gf_internal_t *) h->base_gf->scratch;
    w64_composite_int_t *comp_int = (w64_composite_int_t*)h->private;

    if (base_h->mult_type == GF_MULT_COMPOSITE) {
      comp_int->s = GF_S_GF_16_2_2; 
    } else {
      comp_int->s = GF_S_GF_32_2; 
    }
  } 

  gf->multiply.w64 = gf_w64_composite_multiply;
  gf->divide.w64 = gf_w64_composite_divide;
  gf->inverse.w64 = gf_w64_composite_inverse;

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
    case GF_MULT_COMPOSITE:
      if (region_type & ~(GF_REGION_ALTMAP | GF_REGION_STDMAP)) return -1;
      if (arg1 == 2 && arg2 == 0 || arg1 == 2 && arg2 == 1) {
        return sizeof(gf_internal_t) + sizeof(w64_composite_int_t) + 4;
      } else {
        return -1;
      }
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
    case GF_MULT_COMPOSITE: if (gf_w64_composite_init(gf) == 0) return 0; break;
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

  if (h->region_type & GF_REGION_ALTMAP) {
    if (h->mult_type == GF_MULT_COMPOSITE) {
      gf->extract_word.w64 = gf_w64_composite_extract_word;
    }
  } else {
    gf->extract_word.w64 = gf_w64_extract_word;
  }


  return 1;
}

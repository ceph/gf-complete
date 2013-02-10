/*
 * gf_w16.c
 *
 * Routines for 16-bit Galois fields
 */

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

#define GF_FIELD_WIDTH (16)
#define GF_FIELD_SIZE (1 << GF_FIELD_WIDTH)
#define GF_MULT_GROUP_SIZE GF_FIELD_SIZE-1

#define GF_BASE_FIELD_WIDTH (8)
#define GF_BASE_FIELD_SIZE       (1 << GF_BASE_FIELD_WIDTH)
#define GF_S_GF_8_2 (63)

struct gf_logtable_data {
    int              log_tbl[GF_FIELD_SIZE];
    uint16_t      antilog_tbl[GF_FIELD_SIZE * 2];
    uint16_t      inv_tbl[GF_FIELD_SIZE];
};

struct gf_zero_logtable_data {
    int              log_tbl[GF_FIELD_SIZE];
    uint16_t      _antilog_tbl[GF_FIELD_SIZE * 4];
    uint16_t      *antilog_tbl;
    uint16_t      inv_tbl[GF_FIELD_SIZE];
};

struct gf_lazytable_data {
    int           log_tbl[GF_FIELD_SIZE];
    uint16_t      antilog_tbl[GF_FIELD_SIZE * 2];
    uint16_t      inv_tbl[GF_FIELD_SIZE];
    uint16_t      lazytable[GF_FIELD_SIZE];
};

struct gf_w8_logtable_data {
    uint8_t      log_tbl[GF_BASE_FIELD_SIZE];
    uint8_t      antilog_tbl[GF_BASE_FIELD_SIZE * 2];
    uint8_t      *antilog_tbl_div;
};

struct gf_w8_single_table_data {
    uint8_t      mult[GF_BASE_FIELD_SIZE][GF_BASE_FIELD_SIZE];
};

struct gf_w16_bytwo_data {
    uint64_t prim_poly;
    uint64_t mask1;
    uint64_t mask2;
};

struct gf_w16_group_4_4_data {
    uint16_t reduce[16];
    uint16_t shift[16];
};

#define AB2(ip, am1 ,am2, b, t1, t2) {\
  t1 = (b << 1) & am1;\
  t2 = b & am2; \
  t2 = ((t2 << 1) - (t2 >> (GF_FIELD_WIDTH-1))); \
  b = (t1 ^ (t2 & ip));}

#define SSE_AB2(pp, m1 ,m2, va, t1, t2) {\
          t1 = _mm_and_si128(_mm_slli_epi64(va, 1), m1); \
          t2 = _mm_and_si128(va, m2); \
          t2 = _mm_sub_epi64 (_mm_slli_epi64(t2, 1), _mm_srli_epi64(t2, (GF_FIELD_WIDTH-1))); \
          va = _mm_xor_si128(t1, _mm_and_si128(t2, pp)); }

#define MM_PRINT(s, r) { uint8_t blah[16], ii; printf("%-12s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 2) printf("  %02x %02x", blah[15-ii], blah[14-ii]); printf("\n"); }

static
inline
gf_val_32_t gf_w16_inverse_from_divide (gf_t *gf, gf_val_32_t a)
{
  return gf->divide.w32(gf, 1, a);
}

static
inline
gf_val_32_t gf_w16_divide_from_inverse (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  b = gf->inverse.w32(gf, b);
  return gf->multiply.w32(gf, a, b);
}

static
void
gf_w16_multiply_region_from_single(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  gf_region_data rd;
  uint16_t *s16;
  uint16_t *d16;
  
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);
  gf_do_initial_region_alignment(&rd);

  s16 = (uint16_t *) rd.s_start;
  d16 = (uint16_t *) rd.d_start;

  if (xor) {
    while (d16 < ((uint16_t *) rd.d_top)) {
      *d16 ^= gf->multiply.w32(gf, val, *s16);
      d16++;
      s16++;
    } 
  } else {
    while (d16 < ((uint16_t *) rd.d_top)) {
      *d16 = gf->multiply.w32(gf, val, *s16);
      d16++;
      s16++;
    } 
  }
  gf_do_final_region_alignment(&rd);
}

static
inline
gf_val_32_t gf_w16_euclid (gf_t *gf, gf_val_32_t b)
{
  gf_val_32_t e_i, e_im1, e_ip1;
  gf_val_32_t d_i, d_im1, d_ip1;
  gf_val_32_t y_i, y_im1, y_ip1;
  gf_val_32_t c_i;

  if (b == 0) return -1;
  e_im1 = ((gf_internal_t *) (gf->scratch))->prim_poly;
  e_i = b;
  d_im1 = 16;
  for (d_i = d_im1; ((1 << d_i) & e_i) == 0; d_i--) ;
  y_i = 1;
  y_im1 = 0;

  while (e_i != 1) {

    e_ip1 = e_im1;
    d_ip1 = d_im1;
    c_i = 0;

    while (d_ip1 >= d_i) {
      c_i ^= (1 << (d_ip1 - d_i));
      e_ip1 ^= (e_i << (d_ip1 - d_i));
      while ((e_ip1 & (1 << d_ip1)) == 0) d_ip1--;
    }

    y_ip1 = y_im1 ^ gf->multiply.w32(gf, c_i, y_i);
    y_im1 = y_i;
    y_i = y_ip1;

    e_im1 = e_i;
    d_im1 = d_i;
    e_i = e_ip1;
    d_i = d_ip1;
  }

  return y_i;
}

static
gf_val_32_t gf_w16_extract_word(gf_t *gf, void *start, int bytes, int index)
{
  uint16_t *r16, rv;

  r16 = (uint16_t *) start;
  rv = r16[index];
  return rv;
}

static
gf_val_32_t gf_w16_composite_extract_word(gf_t *gf, void *start, int bytes, int index)
{
  int sub_size;
  gf_internal_t *h;
  uint8_t *r8, *top;
  uint16_t a, b, *r16;
  gf_region_data rd;

  h = (gf_internal_t *) gf->scratch;
  gf_set_region_data(&rd, gf, start, start, bytes, 0, 0, 32);
  r16 = (uint16_t *) start;
  if (r16 + index < (uint16_t *) rd.d_start) return r16[index];
  if (r16 + index >= (uint16_t *) rd.d_top) return r16[index];
  index -= (((uint16_t *) rd.d_start) - r16);
  r8 = (uint8_t *) rd.d_start;
  top = (uint8_t *) rd.d_top;
  sub_size = (top-r8)/2;

  a = h->base_gf->extract_word.w32(h->base_gf, r8, sub_size, index);
  b = h->base_gf->extract_word.w32(h->base_gf, r8+sub_size, sub_size, index);
  return (a | (b << 8));
}

static
gf_val_32_t gf_w16_split_extract_word(gf_t *gf, void *start, int bytes, int index)
{
  uint16_t *r16, rv;
  uint8_t *r8;
  gf_region_data rd;

  gf_set_region_data(&rd, gf, start, start, bytes, 0, 0, 32);
  r16 = (uint16_t *) start;
  if (r16 + index < (uint16_t *) rd.d_start) return r16[index];
  if (r16 + index >= (uint16_t *) rd.d_top) return r16[index];
  index -= (((uint16_t *) rd.d_start) - r16);
  r8 = (uint8_t *) rd.d_start;
  r8 += ((index & 0xfffffff0)*2);
  r8 += (index & 0xf);
  rv = (*r8 << 8);
  r8 += 16;
  rv |= *r8;
  return rv;
}

static
inline
gf_val_32_t gf_w16_matrix (gf_t *gf, gf_val_32_t b)
{
  return gf_bitmatrix_inverse(b, 16, ((gf_internal_t *) (gf->scratch))->prim_poly);
}

/* JSP: GF_MULT_SHIFT: The world's dumbest multiplication algorithm.  I only
   include it for completeness.  It does have the feature that it requires no
   extra memory.  
*/

static
inline
gf_val_32_t
gf_w16_shift_multiply (gf_t *gf, gf_val_32_t a16, gf_val_32_t b16)
{
  gf_val_32_t product, i, pp, a, b;
  gf_internal_t *h;
  
  a = a16;
  b = b16;
  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  product = 0;

  for (i = 0; i < GF_FIELD_WIDTH; i++) { 
    if (a & (1 << i)) product ^= (b << i);
  }
  for (i = (GF_FIELD_WIDTH*2-1); i >= GF_FIELD_WIDTH; i--) {
    if (product & (1 << i)) product ^= (pp << (i-GF_FIELD_WIDTH)); 
  }
  return product;
}

static 
int gf_w16_shift_init(gf_t *gf)
{
  gf->multiply.w32 = gf_w16_shift_multiply;
  gf->inverse.w32 = gf_w16_euclid;
  gf->multiply_region.w32 = gf_w16_multiply_region_from_single;
  return 1;
}

/* KMG: GF_MULT_LOGTABLE: */

static
void
gf_w16_log_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint16_t *s16, *d16;
  int lv;
  struct gf_logtable_data *ltd;
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);
  gf_do_initial_region_alignment(&rd);

  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  s16 = (uint16_t *) rd.s_start;
  d16 = (uint16_t *) rd.d_start;

  lv = ltd->log_tbl[val];

  if (xor) {
    while (d16 < (uint16_t *) rd.d_top) {
      *d16 ^= (*s16 == 0 ? 0 : ltd->antilog_tbl[lv + ltd->log_tbl[*s16]]);
      d16++;
      s16++;
    }
  } else {
    while (d16 < (uint16_t *) rd.d_top) {
      *d16 = (*s16 == 0 ? 0 : ltd->antilog_tbl[lv + ltd->log_tbl[*s16]]);
      d16++;
      s16++;
    }
  }
  gf_do_final_region_alignment(&rd);
}

static
inline
gf_val_32_t
gf_w16_log_multiply(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  struct gf_logtable_data *ltd;

  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return (a == 0 || b == 0) ? 0 : ltd->antilog_tbl[ltd->log_tbl[a] + ltd->log_tbl[b]];
}

static
inline
gf_val_32_t
gf_w16_log_divide(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  int log_sum = 0;
  struct gf_logtable_data *ltd;

  if (a == 0 || b == 0) return 0;
  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;

  log_sum = ltd->log_tbl[a] - ltd->log_tbl[b] + (GF_MULT_GROUP_SIZE);
  return (ltd->antilog_tbl[log_sum]);
}

static
gf_val_32_t
gf_w16_log_inverse(gf_t *gf, gf_val_32_t a)
{
  struct gf_logtable_data *ltd;

  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return (ltd->inv_tbl[a]);
}

static
int gf_w16_log_init(gf_t *gf)
{
  gf_internal_t *h;
  struct gf_logtable_data *ltd;
  int i, b;

  h = (gf_internal_t *) gf->scratch;
  ltd = h->private;

  ltd->log_tbl[0] = 0;

  b = 1;
  for (i = 0; i < GF_MULT_GROUP_SIZE; i++) {
      ltd->log_tbl[b] = i;
      ltd->antilog_tbl[i] = b;
      ltd->antilog_tbl[i+GF_MULT_GROUP_SIZE] = b;
      b <<= 1;
      if (b & GF_FIELD_SIZE) {
          b = b ^ h->prim_poly;
      }
  }
  ltd->inv_tbl[0] = 0;  /* Not really, but we need to fill it with something  */
  ltd->inv_tbl[1] = 1;
  for (i = 2; i < GF_FIELD_SIZE; i++) {
    ltd->inv_tbl[i] = ltd->antilog_tbl[GF_MULT_GROUP_SIZE-ltd->log_tbl[i]];
  }

  gf->inverse.w32 = gf_w16_log_inverse;
  gf->divide.w32 = gf_w16_log_divide;
  gf->multiply.w32 = gf_w16_log_multiply;
  gf->multiply_region.w32 = gf_w16_log_multiply_region;

  return 1;
}

/* JSP: GF_MULT_SPLIT_TABLE: Using 8 multiplication tables to leverage SSE instructions.
*/

static
void
gf_w16_split_4_16_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint64_t i, j, a, c, prod;
  uint16_t *s16, *d16, *top;
  gf_internal_t *h;
  uint16_t table[4][16];
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);
  gf_do_initial_region_alignment(&rd);

  h = (gf_internal_t *) gf->scratch;

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      table[i][j] = gf_w16_log_multiply(gf, c, val);
    }
  }

  s16 = (uint16_t *) rd.s_start;
  d16 = (uint16_t *) rd.d_start;
  top = (uint16_t *) rd.d_top;

  while (d16 < top) {
    a = *s16;
    prod = (xor) ? *d16 : 0;
    for (i = 0; i < 4; i++) {
      prod ^= table[i][a&0xf];
      a >>= 4;
    }
    *d16 = prod;
    s16++;
    d16++;
  }
}

static
void
gf_w16_split_8_16_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint64_t j, a, c, prod, *s64, *d64, *top64;
  gf_internal_t *h;
  uint64_t htable[256], ltable[256];
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 8);
  gf_do_initial_region_alignment(&rd);
  
  h = (gf_internal_t *) gf->scratch;

  for (j = 0; j < 256; j++) {
    ltable[j] = gf_w16_log_multiply(gf, j, val);
    htable[j] = gf_w16_log_multiply(gf, (j<<8), val);
  }

  s64 = (uint64_t *) rd.s_start;
  d64 = (uint64_t *) rd.d_start;
  top64 = (uint64_t *) rd.d_top;
  
/* Does Unrolling Matter?  -- Doesn't seem to.
  while (d64 != top64) {
    a = *s64;

    prod = htable[a >> 56];
    a <<= 8;
    prod ^= ltable[a >> 56];
    a <<= 8;
    prod <<= 16;

    prod ^= htable[a >> 56];
    a <<= 8;
    prod ^= ltable[a >> 56];
    a <<= 8;
    prod <<= 16;

    prod ^= htable[a >> 56];
    a <<= 8;
    prod ^= ltable[a >> 56];
    a <<= 8;
    prod <<= 16;

    prod ^= htable[a >> 56];
    a <<= 8;
    prod ^= ltable[a >> 56];
    prod ^= ((xor) ? *d64 : 0); 
    *d64 = prod;
    *s64++;
    *d64++;
  }
*/
  
  while (d64 != top64) {
    a = *s64;

    prod = 0;
    for (j = 0; j < 4; j++) {
      prod <<= 16;
      prod ^= htable[a >> 56];
      a <<= 8;
      prod ^= ltable[a >> 56];
      a <<= 8;
    }

    prod ^= ((xor) ? *d64 : 0); 
    *d64 = prod;
    *s64++;
    *d64++;
  }
  gf_do_final_region_alignment(&rd);
}

static void
gf_w16_table_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint64_t j, a, c, pp;
  gf_internal_t *h;
  struct gf_lazytable_data *ltd;
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 8);
  gf_do_initial_region_alignment(&rd);

  h = (gf_internal_t *) gf->scratch;
  ltd = (struct gf_lazytable_data *) h->private;

  ltd->lazytable[0] = 0;

  /*
  a = val;
  c = 1;
  pp = h->prim_poly;

  do {
    ltd->lazytable[c] = a;
    c <<= 1;
    if (c & (1 << GF_FIELD_WIDTH)) c ^= pp;
    a <<= 1;
    if (a & (1 << GF_FIELD_WIDTH)) a ^= pp;
  } while (c != 1);
  */

  a = ltd->log_tbl[val];
  for (c = 1; c < GF_FIELD_SIZE; c++) {
    ltd->lazytable[c] = ltd->antilog_tbl[ltd->log_tbl[c]+a];
  }
   
  gf_two_byte_region_table_multiply(&rd, ltd->lazytable);
  gf_do_final_region_alignment(&rd);
}

static
void
gf_w16_split_4_16_lazy_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  uint64_t i, j, *s64, *d64, *top64;;
  uint64_t a, c, prod;
  uint8_t low[4][16];
  uint8_t high[4][16];
  gf_region_data rd;

  __m128i  mask, ta, tb, ti, tpl, tph, tlow[4], thigh[4], tta, ttb, shuffler, unshuffler, lmask;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 32);
  gf_do_initial_region_alignment(&rd);

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      prod = gf_w16_log_multiply(gf, c, val);
      low[i][j] = (prod & 0xff);
      high[i][j] = (prod >> 8);
    }
  }

  for (i = 0; i < 4; i++) {
    tlow[i] = _mm_loadu_si128((__m128i *)low[i]);
    thigh[i] = _mm_loadu_si128((__m128i *)high[i]);
  }

  s64 = (uint64_t *) rd.s_start;
  d64 = (uint64_t *) rd.d_start;
  top64 = (uint64_t *) rd.d_top;

  mask = _mm_set1_epi8 (0x0f);
  lmask = _mm_set1_epi16 (0xff);

  if (xor) {
    while (d64 != top64) {
      
      ta = _mm_load_si128((__m128i *) s64);
      tb = _mm_load_si128((__m128i *) (s64+2));

      tta = _mm_srli_epi16(ta, 8);
      ttb = _mm_srli_epi16(tb, 8);
      tpl = _mm_and_si128(tb, lmask);
      tph = _mm_and_si128(ta, lmask);

      tb = _mm_packus_epi16(tpl, tph);
      ta = _mm_packus_epi16(ttb, tta);

      ti = _mm_and_si128 (mask, tb);
      tph = _mm_shuffle_epi8 (thigh[0], ti);
      tpl = _mm_shuffle_epi8 (tlow[0], ti);
  
      tb = _mm_srli_epi16(tb, 4);
      ti = _mm_and_si128 (mask, tb);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[1], ti), tph);

      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[2], ti), tph);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[3], ti), tph);

      ta = _mm_unpackhi_epi8(tpl, tph);
      tb = _mm_unpacklo_epi8(tpl, tph);

      tta = _mm_load_si128((__m128i *) d64);
      ta = _mm_xor_si128(ta, tta);
      ttb = _mm_load_si128((__m128i *) (d64+2));
      tb = _mm_xor_si128(tb, ttb); 
      _mm_store_si128 ((__m128i *)d64, ta);
      _mm_store_si128 ((__m128i *)(d64+2), tb);

      d64 += 4;
      s64 += 4;
      
    }
  } else {
    while (d64 != top64) {
      
      ta = _mm_load_si128((__m128i *) s64);
      tb = _mm_load_si128((__m128i *) (s64+2));

      tta = _mm_srli_epi16(ta, 8);
      ttb = _mm_srli_epi16(tb, 8);
      tpl = _mm_and_si128(tb, lmask);
      tph = _mm_and_si128(ta, lmask);

      tb = _mm_packus_epi16(tpl, tph);
      ta = _mm_packus_epi16(ttb, tta);

      ti = _mm_and_si128 (mask, tb);
      tph = _mm_shuffle_epi8 (thigh[0], ti);
      tpl = _mm_shuffle_epi8 (tlow[0], ti);
  
      tb = _mm_srli_epi16(tb, 4);
      ti = _mm_and_si128 (mask, tb);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[1], ti), tph);

      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[2], ti), tph);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[3], ti), tph);

      ta = _mm_unpackhi_epi8(tpl, tph);
      tb = _mm_unpacklo_epi8(tpl, tph);

      _mm_store_si128 ((__m128i *)d64, ta);
      _mm_store_si128 ((__m128i *)(d64+2), tb);

      d64 += 4;
      s64 += 4;
    }
  }

  gf_do_final_region_alignment(&rd);
#endif
}

static
void
gf_w16_split_4_16_lazy_sse_altmap_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  uint64_t i, j, *s64, *d64, *top64;;
  uint64_t c, prod;
  uint8_t low[4][16];
  uint8_t high[4][16];
  gf_region_data rd;
  struct gf_single_table_data *std;
  __m128i  mask, ta, tb, ti, tpl, tph, tlow[4], thigh[4];

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 32);
  gf_do_initial_region_alignment(&rd);

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      prod = gf_w16_log_multiply(gf, c, val);
      low[i][j] = (prod & 0xff);
      high[i][j] = (prod >> 8);
    }
  }

  for (i = 0; i < 4; i++) {
    tlow[i] = _mm_loadu_si128((__m128i *)low[i]);
    thigh[i] = _mm_loadu_si128((__m128i *)high[i]);
  }

  s64 = (uint64_t *) rd.s_start;
  d64 = (uint64_t *) rd.d_start;
  top64 = (uint64_t *) rd.d_top;

  mask = _mm_set1_epi8 (0x0f);

  if (xor) {
    while (d64 != top64) {

      ta = _mm_load_si128((__m128i *) s64);
      tb = _mm_load_si128((__m128i *) (s64+2));

      ti = _mm_and_si128 (mask, tb);
      tph = _mm_shuffle_epi8 (thigh[0], ti);
      tpl = _mm_shuffle_epi8 (tlow[0], ti);
  
      tb = _mm_srli_epi16(tb, 4);
      ti = _mm_and_si128 (mask, tb);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[1], ti), tph);

      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[2], ti), tph);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[3], ti), tph);

      ta = _mm_load_si128((__m128i *) d64);
      tph = _mm_xor_si128(tph, ta);
      _mm_store_si128 ((__m128i *)d64, tph);
      tb = _mm_load_si128((__m128i *) (d64+2));
      tpl = _mm_xor_si128(tpl, tb);
      _mm_store_si128 ((__m128i *)(d64+2), tpl);

      d64 += 4;
      s64 += 4;
    }
  } else {
    while (d64 != top64) {

      ta = _mm_load_si128((__m128i *) s64);
      tb = _mm_load_si128((__m128i *) (s64+2));

      ti = _mm_and_si128 (mask, tb);
      tph = _mm_shuffle_epi8 (thigh[0], ti);
      tpl = _mm_shuffle_epi8 (tlow[0], ti);
  
      tb = _mm_srli_epi16(tb, 4);
      ti = _mm_and_si128 (mask, tb);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[1], ti), tph);

      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[2], ti), tph);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tpl = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tpl);
      tph = _mm_xor_si128(_mm_shuffle_epi8 (thigh[3], ti), tph);

      _mm_store_si128 ((__m128i *)d64, tph);
      _mm_store_si128 ((__m128i *)(d64+2), tpl);

      d64 += 4;
      s64 += 4;
      
    }
  }
  gf_do_final_region_alignment(&rd);

#endif
}

static 
int gf_w16_split_init(gf_t *gf)
{
  gf_internal_t *h;
  gf_w16_log_init(gf);

  h = (gf_internal_t *) gf->scratch;

  if (h->mult_type == GF_MULT_DEFAULT) {
    gf->multiply_region.w32 = gf_w16_split_8_16_lazy_multiply_region;
#ifdef INTEL_SSE4
    gf->multiply_region.w32 = gf_w16_split_4_16_lazy_sse_multiply_region;
#endif
  } else if ((h->arg1 == 8 && h->arg2 == 16) || (h->arg2 == 8 && h->arg1 == 16)) {
    gf->multiply_region.w32 = gf_w16_split_8_16_lazy_multiply_region;
  } else if ((h->arg1 == 4 && h->arg2 == 16) || (h->arg2 == 4 && h->arg1 == 16)) {
    if (h->region_type & GF_REGION_SSE) {
      if (h->region_type & GF_REGION_ALTMAP) {
        gf->multiply_region.w32 = gf_w16_split_4_16_lazy_sse_altmap_multiply_region;
      } else {
        gf->multiply_region.w32 = gf_w16_split_4_16_lazy_sse_multiply_region;
      }
    } else {
      gf->multiply_region.w32 = gf_w16_split_4_16_lazy_multiply_region;
    }
  }
  return 1;
}

static 
int gf_w16_table_init(gf_t *gf)
{
  gf_internal_t *h;
  gf_w16_log_init(gf);

  h = (gf_internal_t *) gf->scratch;
  gf->multiply_region.w32 = NULL;
  gf->multiply_region.w32 = gf_w16_table_lazy_multiply_region; 
  return 1;
}

static
void
gf_w16_log_zero_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint16_t lv;
  int i;
  uint16_t *s16, *d16, *top16;
  struct gf_zero_logtable_data *ltd;
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);
  gf_do_initial_region_alignment(&rd);

  ltd = (struct gf_zero_logtable_data*) ((gf_internal_t *) gf->scratch)->private;
  s16 = (uint16_t *) rd.s_start;
  d16 = (uint16_t *) rd.d_start;
  top16 = (uint16_t *) rd.d_top;
  bytes = top16 - d16;

  lv = ltd->log_tbl[val];

  if (xor) {
    for (i = 0; i < bytes; i++) {
      d16[i] ^= (ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]]);
    }
  } else {
    for (i = 0; i < bytes; i++) {
      d16[i] = (ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]]);
    }
  }

  /* This isn't necessary. */
  gf_do_final_region_alignment(&rd);
}

/* Here -- double-check Kevin */
static
inline
gf_val_32_t
gf_w16_log_zero_multiply (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  struct gf_zero_logtable_data *ltd;

  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return ltd->antilog_tbl[ltd->log_tbl[a] + ltd->log_tbl[b]];
}

static
inline
gf_val_32_t
gf_w16_log_zero_divide (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  int log_sum = 0;
  struct gf_zero_logtable_data *ltd;

  if (a == 0 || b == 0) return 0;
  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;

  log_sum = ltd->log_tbl[a] - ltd->log_tbl[b] + (GF_MULT_GROUP_SIZE);
  return (ltd->antilog_tbl[log_sum]);
}

static
gf_val_32_t
gf_w16_log_zero_inverse (gf_t *gf, gf_val_32_t a)
{
  struct gf_zero_logtable_data *ltd;

  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return (ltd->inv_tbl[a]);
}

static
inline
gf_val_32_t
gf_w16_bytwo_p_multiply (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  uint32_t prod, pp, pmask, amask;
  gf_internal_t *h;
  
  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  
  prod = 0;
  pmask = 0x8000;
  amask = 0x8000;

  while (amask != 0) {
    if (prod & pmask) {
      prod = ((prod << 1) ^ pp);
    } else {
      prod <<= 1;
    }
    if (a & amask) prod ^= b;
    amask >>= 1;
  }
  return prod;
}

static
inline
gf_val_32_t
gf_w16_bytwo_b_multiply (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  uint32_t prod, pp, bmask;
  gf_internal_t *h;
  
  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  prod = 0;
  bmask = 0x8000;

  while (1) {
    if (a & 1) prod ^= b;
    a >>= 1;
    if (a == 0) return prod;
    if (b & bmask) {
      b = ((b << 1) ^ pp);
    } else {
      b <<= 1;
    }
  }
}

static
void 
gf_w16_bytwo_p_nosse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint64_t *s64, *d64, t1, t2, ta, prod, amask;
  gf_region_data rd;
  struct gf_w16_bytwo_data *btd;
    
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  btd = (struct gf_w16_bytwo_data *) ((gf_internal_t *) (gf->scratch))->private;

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 8);
  gf_do_initial_region_alignment(&rd);

  s64 = (uint64_t *) rd.s_start;
  d64 = (uint64_t *) rd.d_start;

  if (xor) {
    while (s64 < (uint64_t *) rd.s_top) {
      prod = 0;
      amask = 0x8000;
      ta = *s64;
      while (amask != 0) {
        AB2(btd->prim_poly, btd->mask1, btd->mask2, prod, t1, t2);
        if (val & amask) prod ^= ta;
        amask >>= 1;
      }
      *d64 ^= prod;
      d64++;
      s64++;
    }
  } else { 
    while (s64 < (uint64_t *) rd.s_top) {
      prod = 0;
      amask = 0x8000;
      ta = *s64;
      while (amask != 0) {
        AB2(btd->prim_poly, btd->mask1, btd->mask2, prod, t1, t2);
        if (val & amask) prod ^= ta;
        amask >>= 1;
      }
      *d64 = prod;
      d64++;
      s64++;
    }
  }
  gf_do_final_region_alignment(&rd);
}

#define BYTWO_P_ONESTEP {\
      SSE_AB2(pp, m1 ,m2, prod, t1, t2); \
      t1 = _mm_and_si128(v, one); \
      t1 = _mm_sub_epi16(t1, one); \
      t1 = _mm_and_si128(t1, ta); \
      prod = _mm_xor_si128(prod, t1); \
      v = _mm_srli_epi64(v, 1); }

static
void 
gf_w16_bytwo_p_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  int i;
  uint8_t *s8, *d8;
  uint32_t vrev;
  uint64_t amask;
  __m128i pp, m1, m2, ta, prod, t1, t2, tp, one, v;
  struct gf_w16_bytwo_data *btd;
  gf_region_data rd;
    
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  btd = (struct gf_w16_bytwo_data *) ((gf_internal_t *) (gf->scratch))->private;

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 16);
  gf_do_initial_region_alignment(&rd);

  vrev = 0;
  for (i = 0; i < 16; i++) {
    vrev <<= 1;
    if (!(val & (1 << i))) vrev |= 1;
  }

  s8 = (uint8_t *) rd.s_start;
  d8 = (uint8_t *) rd.d_start;

  pp = _mm_set1_epi16(btd->prim_poly&0xffff);
  m1 = _mm_set1_epi16((btd->mask1)&0xffff);
  m2 = _mm_set1_epi16((btd->mask2)&0xffff);
  one = _mm_set1_epi16(1);

  while (d8 < (uint8_t *) rd.d_top) {
    prod = _mm_setzero_si128();
    v = _mm_set1_epi16(vrev);
    ta = _mm_load_si128((__m128i *) s8);
    tp = (!xor) ? _mm_setzero_si128() : _mm_load_si128((__m128i *) d8);
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    BYTWO_P_ONESTEP;
    _mm_store_si128((__m128i *) d8, _mm_xor_si128(prod, tp));
    d8 += 16;
    s8 += 16;
  }
  gf_do_final_region_alignment(&rd);
#endif
}

static
void
gf_w16_bytwo_b_sse_region_2_noxor(gf_region_data *rd, struct gf_w16_bytwo_data *btd)
{
#ifdef   INTEL_SSE4
  int i;
  uint8_t *d8, *s8, tb;
  __m128i pp, m1, m2, t1, t2, va, vb;

  s8 = (uint8_t *) rd->s_start;
  d8 = (uint8_t *) rd->d_start;

  pp = _mm_set1_epi16(btd->prim_poly&0xffff);
  m1 = _mm_set1_epi16((btd->mask1)&0xffff);
  m2 = _mm_set1_epi16((btd->mask2)&0xffff);

  while (d8 < (uint8_t *) rd->d_top) {
    va = _mm_load_si128 ((__m128i *)(s8));
    SSE_AB2(pp, m1, m2, va, t1, t2);
    _mm_store_si128((__m128i *)d8, va);
    d8 += 16;
    s8 += 16;
  }
#endif
}

static
void
gf_w16_bytwo_b_sse_region_2_xor(gf_region_data *rd, struct gf_w16_bytwo_data *btd)
{
#ifdef   INTEL_SSE4
  int i;
  uint8_t *d8, *s8, tb;
  __m128i pp, m1, m2, t1, t2, va, vb;

  s8 = (uint8_t *) rd->s_start;
  d8 = (uint8_t *) rd->d_start;

  pp = _mm_set1_epi16(btd->prim_poly&0xffff);
  m1 = _mm_set1_epi16((btd->mask1)&0xffff);
  m2 = _mm_set1_epi16((btd->mask2)&0xffff);

  while (d8 < (uint8_t *) rd->d_top) {
    va = _mm_load_si128 ((__m128i *)(s8));
    SSE_AB2(pp, m1, m2, va, t1, t2);
    vb = _mm_load_si128 ((__m128i *)(d8));
    vb = _mm_xor_si128(vb, va);
    _mm_store_si128((__m128i *)d8, vb);
    d8 += 16;
    s8 += 16;
  }
#endif
}


static
void 
gf_w16_bytwo_b_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  int itb;
  uint8_t *d8, *s8;
  __m128i pp, m1, m2, t1, t2, va, vb;
  struct gf_w16_bytwo_data *btd;
  gf_region_data rd;
    
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 16);
  gf_do_initial_region_alignment(&rd);

  btd = (struct gf_w16_bytwo_data *) ((gf_internal_t *) (gf->scratch))->private;

  if (val == 2) {
    if (xor) {
      gf_w16_bytwo_b_sse_region_2_xor(&rd, btd);
    } else {
      gf_w16_bytwo_b_sse_region_2_noxor(&rd, btd);
    }
    gf_do_final_region_alignment(&rd);
    return;
  }

  s8 = (uint8_t *) rd.s_start;
  d8 = (uint8_t *) rd.d_start;

  pp = _mm_set1_epi16(btd->prim_poly&0xffff);
  m1 = _mm_set1_epi16((btd->mask1)&0xffff);
  m2 = _mm_set1_epi16((btd->mask2)&0xffff);

  while (d8 < (uint8_t *) rd.d_top) {
    va = _mm_load_si128 ((__m128i *)(s8));
    vb = (!xor) ? _mm_setzero_si128() : _mm_load_si128 ((__m128i *)(d8));
    itb = val;
    while (1) {
      if (itb & 1) vb = _mm_xor_si128(vb, va);
      itb >>= 1;
      if (itb == 0) break;
      SSE_AB2(pp, m1, m2, va, t1, t2);
    }
    _mm_store_si128((__m128i *)d8, vb);
    d8 += 16;
    s8 += 16;
  }

  gf_do_final_region_alignment(&rd);
#endif
}

static
void 
gf_w16_bytwo_b_nosse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  int i;
  uint64_t *s64, *d64, t1, t2, ta, tb, prod;
  struct gf_w16_bytwo_data *btd;
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 16);
  gf_do_initial_region_alignment(&rd);

  btd = (struct gf_w16_bytwo_data *) ((gf_internal_t *) (gf->scratch))->private;
  s64 = (uint64_t *) rd.s_start;
  d64 = (uint64_t *) rd.d_start;

  switch (val) {
  case 2:
    if (xor) {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 ^= ta;
        d64++;
        s64++;
      }
    } else {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 = ta;
        d64++;
        s64++;
      }
    }
    break; 
  case 3:
    if (xor) {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        prod = ta;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 ^= (ta ^ prod);
        d64++;
        s64++;
      }
    } else {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        prod = ta;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 = (ta ^ prod);
        d64++;
        s64++;
      }
    }
    break; 
  case 4:
    if (xor) {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 ^= ta;
        d64++;
        s64++;
      }
    } else {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 = ta;
        d64++;
        s64++;
      }
    }
    break; 
  case 5:
    if (xor) {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        prod = ta;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 ^= (ta ^ prod);
        d64++;
        s64++;
      }
    } else {
      while (d64 < (uint64_t *) rd.d_top) {
        ta = *s64;
        prod = ta;
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        *d64 = ta ^ prod;
        d64++;
        s64++;
      }
    }
  default:
    if (xor) {
      while (d64 < (uint64_t *) rd.d_top) {
        prod = *d64 ;
        ta = *s64;
        tb = val;
        while (1) {
          if (tb & 1) prod ^= ta;
          tb >>= 1;
          if (tb == 0) break;
          AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        }
        *d64 = prod;
        d64++;
        s64++;
      }
    } else {
      while (d64 < (uint64_t *) rd.d_top) {
        prod = 0 ;
        ta = *s64;
        tb = val;
        while (1) {
          if (tb & 1) prod ^= ta;
          tb >>= 1;
          if (tb == 0) break;
          AB2(btd->prim_poly, btd->mask1, btd->mask2, ta, t1, t2);
        }
        *d64 = prod;
        d64++;
        s64++;
      }
    }
    break;
  }
  gf_do_final_region_alignment(&rd);
}

static
int gf_w16_bytwo_init(gf_t *gf)
{
  gf_internal_t *h;
  uint64_t ip, m1, m2;
  struct gf_w16_bytwo_data *btd;

  h = (gf_internal_t *) gf->scratch;
  btd = (struct gf_w16_bytwo_data *) (h->private);
  ip = h->prim_poly & 0xffff;
  m1 = 0xfffe;
  m2 = 0x8000;
  btd->prim_poly = 0;
  btd->mask1 = 0;
  btd->mask2 = 0;

  while (ip != 0) {
    btd->prim_poly |= ip;
    btd->mask1 |= m1;
    btd->mask2 |= m2;
    ip <<= GF_FIELD_WIDTH;
    m1 <<= GF_FIELD_WIDTH;
    m2 <<= GF_FIELD_WIDTH;
  }

  if (h->mult_type == GF_MULT_BYTWO_p) {
    gf->multiply.w32 = gf_w16_bytwo_p_multiply;
    if (h->region_type == GF_REGION_SSE) {
      gf->multiply_region.w32 = gf_w16_bytwo_p_sse_multiply_region;
    } else {
      gf->multiply_region.w32 = gf_w16_bytwo_p_nosse_multiply_region;
    }
  } else {
    gf->multiply.w32 = gf_w16_bytwo_b_multiply;
    if (h->region_type == GF_REGION_SSE) {
      gf->multiply_region.w32 = gf_w16_bytwo_b_sse_multiply_region;
    } else {
      gf->multiply_region.w32 = gf_w16_bytwo_b_nosse_multiply_region;
    }
  }
  gf->inverse.w32 = gf_w16_euclid;
  return 1;
}

static
int gf_w16_log_zero_init(gf_t *gf)
{
  gf_internal_t *h;
  struct gf_zero_logtable_data *ltd;
  int i, b;

  h = (gf_internal_t *) gf->scratch;
  ltd = h->private;

  ltd->log_tbl[0] = (-GF_MULT_GROUP_SIZE) + 1;

  bzero(&(ltd->_antilog_tbl[0]), sizeof(ltd->_antilog_tbl));

  ltd->antilog_tbl = &(ltd->_antilog_tbl[GF_FIELD_SIZE * 2]);

  b = 1;
  for (i = 0; i < GF_MULT_GROUP_SIZE; i++) {
      ltd->log_tbl[b] = (uint16_t)i;
      ltd->antilog_tbl[i] = (uint16_t)b;
      ltd->antilog_tbl[i+GF_MULT_GROUP_SIZE] = (uint16_t)b;
      b <<= 1;
      if (b & GF_FIELD_SIZE) {
          b = b ^ h->prim_poly;
      }
  }
  ltd->inv_tbl[0] = 0;  /* Not really, but we need to fill it with something  */
  ltd->inv_tbl[1] = 1;
  for (i = 2; i < GF_FIELD_SIZE; i++) {
    ltd->inv_tbl[i] = ltd->antilog_tbl[GF_MULT_GROUP_SIZE-ltd->log_tbl[i]];
  }

  gf->inverse.w32 = gf_w16_log_zero_inverse;
  gf->divide.w32 = gf_w16_log_zero_divide;
  gf->multiply.w32 = gf_w16_log_zero_multiply;
  gf->multiply_region.w32 = gf_w16_log_zero_multiply_region;
  return 1;
}

static
gf_val_32_t
gf_w16_composite_multiply_recursive(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t b0 = b & 0x00ff;
  uint8_t b1 = (b & 0xff00) >> 8;
  uint8_t a0 = a & 0x00ff;
  uint8_t a1 = (a & 0xff00) >> 8;
  uint8_t a1b1;
  uint16_t rv;

  a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

  rv = ((base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) | ((base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, GF_S_GF_8_2)) << 8));
  return rv;
}

static
gf_val_32_t
gf_w16_composite_multiply_table(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  struct gf_w8_single_table_data * std;

  uint8_t b0 = b & 0x00ff;
  uint8_t b1 = (b & 0xff00) >> 8;
  uint8_t a0 = a & 0x00ff;
  uint8_t a1 = (a & 0xff00) >> 8;
  uint8_t a1b1;
  uint16_t rv;

  std = (struct gf_w8_single_table_data *) h->private;

  a1b1 = std->mult[a1][b1];

  rv = ((std->mult[a0][b0] ^ a1b1) | 
        ((std->mult[a1][b0] ^ std->mult[a0][b1] ^ std->mult[a1b1][GF_S_GF_8_2]) << 8));
  return rv;
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
gf_val_32_t
gf_w16_composite_inverse(gf_t *gf, gf_val_32_t a)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t a0 = a & 0x00ff;
  uint8_t a1 = (a & 0xff00) >> 8;
  uint8_t c0, c1, d, tmp;
  uint16_t c;
  uint8_t a0inv, a1inv;

  if (a0 == 0) {
    a1inv = base_gf->inverse.w32(base_gf, a1);
    c0 = base_gf->multiply.w32(base_gf, a1inv, GF_S_GF_8_2);
    c1 = a1inv;
  } else if (a1 == 0) {
    c0 = base_gf->inverse.w32(base_gf, a0);
    c1 = 0;
  } else {
    a1inv = base_gf->inverse.w32(base_gf, a1);
    a0inv = base_gf->inverse.w32(base_gf, a0);

    d = base_gf->multiply.w32(base_gf, a1, a0inv);

    tmp = (base_gf->multiply.w32(base_gf, a1, a0inv) ^ base_gf->multiply.w32(base_gf, a0, a1inv) ^ GF_S_GF_8_2);
    tmp = base_gf->inverse.w32(base_gf, tmp);

    d = base_gf->multiply.w32(base_gf, d, tmp);

    c0 = base_gf->multiply.w32(base_gf, (d^1), a0inv);
    c1 = base_gf->multiply.w32(base_gf, d, a1inv);
  }

  c = c0 | (c1 << 8);

  return c;
}

static
gf_val_32_t
gf_w16_composite_divide(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  uint16_t binv;

  binv = gf->inverse.w32(gf, b);
  return gf->multiply.w32(gf, a, binv);
}

static
void
gf_w16_composite_multiply_region_inline(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  struct gf_w8_single_table_data * std;
  uint8_t b0 = val & 0x00ff;
  uint8_t b1 = (val & 0xff00) >> 8;
  uint16_t *s16, *d16, *top;
  uint8_t a0, a1, a1b1;
  struct gf_logtable_data *ltd;
  gf_region_data rd;
  
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }

  std = (struct gf_w8_single_table_data *) h->private;
  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);

  s16 = rd.s_start;
  d16 = rd.d_start;
  top = rd.d_top;

  if (xor) {
    while (d16 < top) {
      a0 = (*s16) & 0x00ff;
      a1 = ((*s16) & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      *d16 ^= ((std->mult[a0][b0] ^ a1b1) | ((std->mult[a1][b0] ^ std->mult[a0][b1] ^ std->mult[a1b1][GF_S_GF_8_2]) << 8));
      s16++;
      d16++;
    }
  } else {
    while (d16 < top) {
      a0 = (*s16) & 0x00ff;
      a1 = ((*s16) & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      *d16 = ((std->mult[a0][b0] ^ a1b1) | ((std->mult[a1][b0] ^ std->mult[a0][b1] ^ std->mult[a1b1][GF_S_GF_8_2]) << 8));
      s16++;
      d16++;
    }
  }
}

static
void
gf_w16_composite_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t b0 = val & 0x00ff;
  uint8_t b1 = (val & 0xff00) >> 8;
  uint16_t *s16, *d16, *top;
  uint8_t a0, a1, a1b1;
  gf_region_data rd;
  struct gf_logtable_data *ltd;
  
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);

  s16 = rd.s_start;
  d16 = rd.d_start;
  top = rd.d_top;

  if (xor) {
    while (d16 < top) {
      a0 = (*s16) & 0x00ff;
      a1 = ((*s16) & 0xff00) >> 8;
      a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

      (*d16) ^= ((base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) |
                ((base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, GF_S_GF_8_2)) << 8));
      s16++;
      d16++;
    }
  } else {
    while (d16 < top) {
      a0 = (*s16) & 0x00ff;
      a1 = ((*s16) & 0xff00) >> 8;
      a1b1 = base_gf->multiply.w32(base_gf, a1, b1);

      (*d16) = ((base_gf->multiply.w32(base_gf, a0, b0) ^ a1b1) |
               ((base_gf->multiply.w32(base_gf, a1, b0) ^ base_gf->multiply.w32(base_gf, a0, b1) ^ base_gf->multiply.w32(base_gf, a1b1, GF_S_GF_8_2)) << 8));
      s16++;
      d16++;
    }
  }
}

static
void
gf_w16_composite_multiply_region_alt(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t val0 = val & 0x00ff;
  uint8_t val1 = (val & 0xff00) >> 8;
  gf_region_data rd;
  int sub_reg_size;
  uint8_t *slow, *shigh;
  uint8_t *dlow, *dhigh, *top;;

  /* JSP: I want the two pointers aligned wrt each other on 16 byte 
     boundaries.  So I'm going to make sure that the area on 
     which the two operate is a multiple of 32. Of course, that 
     junks up the mapping, but so be it -- that's why we have extract_word.... */

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 32);
  gf_do_initial_region_alignment(&rd);

  slow = (uint8_t *) rd.s_start;
  dlow = (uint8_t *) rd.d_start;
  top = (uint8_t *)  rd.d_top;
  sub_reg_size = (top - dlow)/2;
  shigh = slow + sub_reg_size;
  dhigh = dlow + sub_reg_size;

  base_gf->multiply_region.w32(base_gf, slow, dlow, val0, sub_reg_size, xor);
  base_gf->multiply_region.w32(base_gf, shigh, dlow, val1, sub_reg_size, 1);
  base_gf->multiply_region.w32(base_gf, slow, dhigh, val1, sub_reg_size, xor);
  base_gf->multiply_region.w32(base_gf, shigh, dhigh, val0, sub_reg_size, 1);
  base_gf->multiply_region.w32(base_gf, shigh, dhigh, base_gf->multiply.w32(base_gf, GF_S_GF_8_2, val1), sub_reg_size, 1);

  gf_do_final_region_alignment(&rd);
}

static
int gf_w16_composite_init(gf_t *gf)
{
  struct gf_w8_single_table_data * std;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_internal_t *base_h = (gf_internal_t *) base_gf->scratch;
  uint16_t a, b;

  if (h->region_type & GF_REGION_ALTMAP) {
    gf->multiply_region.w32 = gf_w16_composite_multiply_region_alt;
  } else if (h->arg2 == 0 && base_h->mult_type == GF_MULT_TABLE && 
                             base_h->region_type == GF_REGION_DEFAULT) {
    gf->multiply_region.w32 = gf_w16_composite_multiply_region_inline;
  } else {
    gf->multiply_region.w32 = gf_w16_composite_multiply_region;
  }
  
  if (h->arg2 == 0) {
    std = (struct gf_w8_single_table_data *) h->private;
    for (a = 0; a < 256; a++) {
      for (b = 0; b < 256; b++) {
        std->mult[a][b] = base_gf->multiply.w32(base_gf, a, b);
      }
    }
    gf->multiply.w32 = gf_w16_composite_multiply_table;
  } else {
    gf->multiply.w32 = gf_w16_composite_multiply_recursive;
  }

  gf->divide.w32 = gf_w16_composite_divide;
  gf->inverse.w32 = gf_w16_composite_inverse;

  return 1;
}

static
void
gf_w16_group_4_set_shift_tables(uint16_t *shift, uint16_t val, gf_internal_t *h)
{
  int i, j;

  shift[0] = 0;
  for (i = 0; i < 16; i += 2) {
    j = (shift[i>>1] << 1);
    if (j & (1 << 16)) j ^= h->prim_poly;
    shift[i] = j;
    shift[i^1] = j^val;
  }
}

static
inline
gf_val_32_t
gf_w16_group_4_4_multiply(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  int i;
  uint16_t p, l, ind, r, a16;

  struct gf_w16_group_4_4_data *d44;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;

  d44 = (struct gf_w16_group_4_4_data *) h->private;
  gf_w16_group_4_set_shift_tables(d44->shift, b, h);

  a16 = a;
  ind = a16 >> 12;
  a16 <<= 4;
  p = d44->shift[ind];
  r = p & 0xfff;
  l = p >> 12;
  ind = a16 >> 12;
  a16 <<= 4;
  p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
  r = p & 0xfff;
  l = p >> 12;
  ind = a16 >> 12;
  a16 <<= 4;
  p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
  r = p & 0xfff;
  l = p >> 12;
  ind = a16 >> 12;
  p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
  return p;
}

static
void gf_w16_group_4_4_region_multiply(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  int i;
  uint16_t p, l, ind, r, a16, p16;
  struct gf_w16_group_4_4_data *d44;
  gf_region_data rd;
  uint16_t *s16, *d16, *top;
  
  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  d44 = (struct gf_w16_group_4_4_data *) h->private;
  gf_w16_group_4_set_shift_tables(d44->shift, val, h);

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 2);
  gf_do_initial_region_alignment(&rd);

  s16 = (uint16_t *) rd.s_start;
  d16 = (uint16_t *) rd.d_start;
  top = (uint16_t *) rd.d_top;

  while (d16 < top) {
    p = 0;
    a16 = *s16;
    p16 = (xor) ? *d16 : 0;
    ind = a16 >> 12;
    a16 <<= 4;
    p = d44->shift[ind];
    r = p & 0xfff;
    l = p >> 12;
    ind = a16 >> 12;
    a16 <<= 4;
    p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
    r = p & 0xfff;
    l = p >> 12;
    ind = a16 >> 12;
    a16 <<= 4;
    p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
    r = p & 0xfff;
    l = p >> 12;
    ind = a16 >> 12;
    p = (d44->shift[ind] ^ d44->reduce[l] ^ (r << 4));
    p ^= p16;
    *d16 = p;
    d16++;
    s16++;
  }
  gf_do_final_region_alignment(&rd);
}

static
int gf_w16_group_init(gf_t *gf)
{
  int i, j, p;
  struct gf_w16_group_4_4_data *d44;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;

  d44 = (struct gf_w16_group_4_4_data *) h->private;
  d44->reduce[0] = 0;
  for (i = 0; i < 16; i++) {
    p = 0;
    for (j = 0; j < 4; j++) {
      if (i & (1 << j)) p ^= (h->prim_poly << j);
    }
    d44->reduce[p>>16] = (p&0xffff);
  }

  gf->multiply.w32 = gf_w16_group_4_4_multiply;
  gf->divide.w32 = NULL;
  gf->inverse.w32 = NULL;
  gf->multiply_region.w32 = gf_w16_group_4_4_region_multiply;

  return 1;
}

int gf_w16_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2)
{
  int ss;
  int sa;

  ss = (GF_REGION_SSE | GF_REGION_NOSSE);
  sa = (GF_REGION_STDMAP | GF_REGION_ALTMAP);

  switch(mult_type)
  {
    case GF_MULT_TABLE:
      region_type |= GF_REGION_LAZY;
      if (arg1 != 0 || arg2 != 0 || region_type != GF_REGION_LAZY) return -1;
      return sizeof(gf_internal_t) + sizeof(struct gf_lazytable_data) + 64;
      break;
    case GF_MULT_BYTWO_p:
    case GF_MULT_BYTWO_b:
      if (arg1 != 0 || arg2 != 0 || (region_type | ss) != ss ||
         (region_type & ss) == ss) return -1;
      return sizeof(gf_internal_t) + sizeof(struct gf_w16_bytwo_data);
      break;
    case GF_MULT_DEFAULT:
    case GF_MULT_LOG_TABLE:
      if (arg2 != 0) return -1;
      if (region_type != GF_REGION_DEFAULT) return -1;
      if (arg1 == 1) {
        return sizeof(gf_internal_t) + sizeof(struct gf_zero_logtable_data) + 64;
      } else if (arg1 == 0) {
        return sizeof(gf_internal_t) + sizeof(struct gf_logtable_data) + 64;
      } else {
        return -1;
      }
      break;
    case GF_MULT_SPLIT_TABLE: 
        if ((arg1 == 8 && arg2 == 16) || (arg2 == 8 && arg1 == 16)) {
          region_type |= GF_REGION_LAZY;
          if (region_type != GF_REGION_LAZY) return -1;
          return sizeof(gf_internal_t) + sizeof(struct gf_logtable_data) + 64;
        } else if ((arg1 == 4 && arg2 == 16) || (arg2 == 4 && arg1 == 16)) {
          region_type &= (~GF_REGION_LAZY);    /* Ignore GF_REGION_LAZY */
          if ((region_type & ss) == ss) return -1;
          if ((region_type & sa) == sa) return -1;
          if ((region_type & ss) == 0) region_type |= GF_REGION_SSE;
          if (region_type & GF_REGION_NOSSE) {
            if (region_type != GF_REGION_NOSSE) return -1;
            return sizeof(gf_internal_t) + sizeof(struct gf_logtable_data) + 64;
          } else {
            if ((region_type | ss | sa) != (ss|sa)) return -1;
            return sizeof(gf_internal_t) + sizeof(struct gf_logtable_data) + 64;
          }
        }
        return -1;
        break;
    case GF_MULT_GROUP:     
      if (arg1 == 4 && arg2 == 4) {
            return sizeof(gf_internal_t) + sizeof(struct gf_w16_group_4_4_data) + 64;
      }
      return -1;
    case GF_MULT_SHIFT:
      if (arg1 != 0 || arg2 != 0 || region_type != 0) return -1;
      return sizeof(gf_internal_t);
      break;
    case GF_MULT_COMPOSITE:
      if (region_type & ~(GF_REGION_ALTMAP | GF_REGION_STDMAP)) return -1;
      if (arg1 == 2 && arg2 == 0) {
        return sizeof(gf_internal_t) + sizeof(struct gf_w8_single_table_data) + 64;
      } else if (arg1 == 2 && arg2 == 1) {
        return sizeof(gf_internal_t) + 64;
      } else {
        return -1;
      }

    default:
      return -1;
   }
}

int gf_w16_init(gf_t *gf)
{
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  if (h->prim_poly == 0) h->prim_poly = 0x1100b;

  gf->multiply.w32 = NULL;
  gf->divide.w32 = NULL;
  gf->inverse.w32 = NULL;
  gf->multiply_region.w32 = NULL;

  switch(h->mult_type) {
    case GF_MULT_LOG_TABLE:        
      if (h->arg1 == 1) {
        if (gf_w16_log_zero_init(gf) == 0) return 0;
      } else {
        if (gf_w16_log_init(gf) == 0) return 0;
      }
      break;
    case GF_MULT_DEFAULT: 
    case GF_MULT_SPLIT_TABLE: if (gf_w16_split_init(gf) == 0) return 0; break;
    case GF_MULT_TABLE:       if (gf_w16_table_init(gf) == 0) return 0; break;
    case GF_MULT_SHIFT:     if (gf_w16_shift_init(gf) == 0) return 0; break;
    case GF_MULT_COMPOSITE: if (gf_w16_composite_init(gf) == 0) return 0; break;
    case GF_MULT_BYTWO_p: 
    case GF_MULT_BYTWO_b:   if (gf_w16_bytwo_init(gf) == 0) return 0; break;
    case GF_MULT_GROUP:     if (gf_w16_group_init(gf) == 0) return 0; break;
    default: return 0;
  }
  if (h->divide_type == GF_DIVIDE_EUCLID) {
    gf->divide.w32 = gf_w16_divide_from_inverse;
    gf->inverse.w32 = gf_w16_euclid;
  } else if (h->divide_type == GF_DIVIDE_MATRIX) {
    gf->divide.w32 = gf_w16_divide_from_inverse;
    gf->inverse.w32 = gf_w16_matrix;
  }

  if (gf->inverse.w32== NULL && gf->divide.w32 == NULL) gf->inverse.w32 = gf_w16_euclid;

  if (gf->inverse.w32 != NULL && gf->divide.w32 == NULL) {
    gf->divide.w32 = gf_w16_divide_from_inverse;
  }
  if (gf->inverse.w32 == NULL && gf->divide.w32 != NULL) {
    gf->inverse.w32 = gf_w16_inverse_from_divide;
  }
  if (h->region_type & GF_REGION_ALTMAP) {
    if (h->mult_type == GF_MULT_COMPOSITE) {
      gf->extract_word.w32 = gf_w16_composite_extract_word;
    } else {
      gf->extract_word.w32 = gf_w16_split_extract_word;
    }
  } else {
    gf->extract_word.w32 = gf_w16_extract_word;
  }
  return 1;
}

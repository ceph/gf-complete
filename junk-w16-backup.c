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
    gf_val_16_t      antilog_tbl[GF_FIELD_SIZE * 2];
    gf_val_16_t      inv_tbl[GF_FIELD_SIZE];
};

struct gf_zero_logtable_data {
    int              log_tbl[GF_FIELD_SIZE];
    gf_val_16_t      _antilog_tbl[GF_FIELD_SIZE * 4];
    gf_val_16_t      *antilog_tbl;
    gf_val_16_t      inv_tbl[GF_FIELD_SIZE];
};

struct gf_lazytable_data {
    int              log_tbl[GF_FIELD_SIZE];
    gf_val_16_t      antilog_tbl[GF_FIELD_SIZE * 2];
    gf_val_16_t      inv_tbl[GF_FIELD_SIZE];
    gf_val_16_t      lazytable[GF_FIELD_SIZE];
};

struct gf_w8_logtable_data {
    gf_val_8_t      log_tbl[GF_BASE_FIELD_SIZE];
    gf_val_8_t      antilog_tbl[GF_BASE_FIELD_SIZE * 2];
    gf_val_8_t      *antilog_tbl_div;
};

struct gf_w8_single_table_data {
    gf_val_8_t      mult[GF_BASE_FIELD_SIZE][GF_BASE_FIELD_SIZE];
    gf_val_8_t      div[GF_BASE_FIELD_SIZE][GF_BASE_FIELD_SIZE];
};

struct gf_w8_double_table_data {
    gf_val_8_t      div[GF_BASE_FIELD_SIZE][GF_BASE_FIELD_SIZE];
    gf_val_8_t      mult[GF_BASE_FIELD_SIZE][GF_BASE_FIELD_SIZE*GF_BASE_FIELD_SIZE];
};


#define MM_PRINT(s, r) { uint8_t blah[16], ii; printf("%-12s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 2) printf("  %02x %02x", blah[15-ii], blah[14-ii]); printf("\n"); }

static
inline
gf_val_16_t gf_w16_inverse_from_divide (gf_t *gf, gf_val_16_t a)
{
  return gf->divide.w16(gf, 1, a);
}

static
inline
gf_val_16_t gf_w16_divide_from_inverse (gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  b = gf->inverse.w16(gf, b);
  return gf->multiply.w16(gf, a, b);
}

static
void
gf_w16_multiply_region_from_single(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  int i;
  gf_val_16_t *s16;
  gf_val_16_t *d16;
  
  s16 = (gf_val_16_t *) src;
  d16 = (gf_val_16_t *) dest;

  if (xor) {
    for (i = 0; i < bytes/2; i++) {
      d16[i] ^= gf->multiply.w16(gf, val, s16[i]);
    } 
  } else {
    for (i = 0; i < bytes/2; i++) {
      d16[i] = gf->multiply.w16(gf, val, s16[i]);
    } 
  }
}

static
inline
gf_val_16_t gf_w16_euclid (gf_t *gf, gf_val_16_t b)
{
  gf_val_32_t e_i, e_im1, e_ip1;
  gf_val_32_t d_i, d_im1, d_ip1;
  gf_val_16_t y_i, y_im1, y_ip1;
  gf_val_16_t c_i;

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

    y_ip1 = y_im1 ^ gf->multiply.w16(gf, c_i, y_i);
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
inline
gf_val_16_t gf_w16_matrix (gf_t *gf, gf_val_16_t b)
{
  return gf_bitmatrix_inverse(b, 16, ((gf_internal_t *) (gf->scratch))->prim_poly);
}

/* JSP: GF_MULT_SHIFT: The world's dumbest multiplication algorithm.  I only
   include it for completeness.  It does have the feature that it requires no
   extra memory.  
*/

static
inline
gf_val_16_t
gf_w16_shift_multiply (gf_t *gf, gf_val_16_t a16, gf_val_16_t b16)
{
  uint32_t product, i, pp, a, b;
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
  gf->multiply.w16 = gf_w16_shift_multiply;
  gf->inverse.w16 = gf_w16_euclid;
  gf->multiply_region.w16 = gf_w16_multiply_region_from_single;
  return 1;
}

/* KMG: GF_MULT_LOGTABLE: */

static
void
gf_w16_log_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  int i;
  uint16_t lv, b, c;
  uint16_t *s16, *d16;
  int num_syms = bytes >> 1;
  int sym_divisible = bytes % 2;

  struct gf_logtable_data *ltd;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w16_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w16_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;

  lv = ltd->log_tbl[val];

  if (xor) {
    for (i = 0; i < num_syms; i++) {
      d16[i] ^= (s16[i] == 0 ? 0 : ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]]);
    }
  } else {
    for (i = 0; i < num_syms; i++) {
      d16[i] = (s16[i] == 0 ? 0 : ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]]);
    }
  }
}

static
inline
gf_val_16_t
gf_w16_log_multiply(gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  struct gf_logtable_data *ltd;

  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return (a == 0 || b == 0) ? 0 : ltd->antilog_tbl[ltd->log_tbl[a] + ltd->log_tbl[b]];
}

static
inline
gf_val_16_t
gf_w16_log_divide(gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  int log_sum = 0;
  struct gf_logtable_data *ltd;

  if (a == 0 || b == 0) return 0;
  ltd = (struct gf_logtable_data *) ((gf_internal_t *) gf->scratch)->private;

  log_sum = ltd->log_tbl[a] - ltd->log_tbl[b] + (GF_MULT_GROUP_SIZE);
  return (ltd->antilog_tbl[log_sum]);
}

static
gf_val_16_t
gf_w16_log_inverse(gf_t *gf, gf_val_16_t a)
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
      ltd->log_tbl[b] = (gf_val_16_t)i;
      ltd->antilog_tbl[i] = (gf_val_16_t)b;
      ltd->antilog_tbl[i+GF_MULT_GROUP_SIZE] = (gf_val_16_t)b;
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

  gf->inverse.w16 = gf_w16_log_inverse;
  gf->divide.w16 = gf_w16_log_divide;
  gf->multiply.w16 = gf_w16_log_multiply;
  gf->multiply_region.w16 = gf_w16_log_multiply_region;

  return 1;
}

/* JSP: GF_MULT_SPLIT_TABLE: Using 8 multiplication tables to leverage SSE instructions.
*/

static
void
gf_w16_split_4_16_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  uint64_t i, j, a, c, prod;
  uint16_t *s16, *d16, *top;
  gf_internal_t *h;
  uint16_t table[4][16];
  
  h = (gf_internal_t *) gf->scratch;

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      table[i][j] = gf_w16_log_multiply(gf, c, val);
    }
  }

  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

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
gf_w16_split_8_16_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  uint64_t j, a, c, prod, *s64, *d64, *top64;
  uint16_t *s16, *d16, *top;
  gf_internal_t *h;
  uint64_t htable[256], ltable[256];
  unsigned long uls, uld;
  
  h = (gf_internal_t *) gf->scratch;

  uls = ((unsigned long) src) & 0xf;
  uld = ((unsigned long) dest) & 0xf;
  if (uls != uld || uls % 2 != 0 || bytes % 2 != 0) gf_alignment_error("gf_w16_split_8_16_lazy_multiply_region", 2);

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  for (j = 0; j < 256; j++) {
    ltable[j] = gf_w16_log_multiply(gf, j, val);
    htable[j] = gf_w16_log_multiply(gf, (j<<8), val);
  }

  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

  if (uls != 0) {
    while (uls != 16 && d16 < top) {
      a = *s16;
      prod = (xor) ? *d16 : 0;
      prod ^= ltable[a&0xff];
      a >>= 8;
      prod ^= htable[a];
      *d16 = prod;
      s16++;
      d16++;
      uls += 2;
    }
    if (d16 == top) return;
  }

  uls = ((unsigned long) top) & 0xf;
  uld = ((unsigned long) top) ^ uls;
  top64 = (uint64_t *) uld;
  s64 = (uint64_t *) s16;
  d64 = (uint64_t *) d16;
  
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

  
  if (uls != 0) {
    d16 = (uint16_t *) d64;
    s16 = (uint16_t *) s64;
    while (d16 < top) {
      a = *s16;
      prod = (xor) ? *d16 : 0;
      prod ^= ltable[a&0xff];
      a >>= 8;
      prod ^= htable[a];
      *d16 = prod;
      s16++;
      d16++;
    }
  }
  return;
}

static
void
gf_w16_table_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  uint64_t j, a, c, prod, *s64, *d64, *top64, pp;
  uint16_t *s16, *d16, *top;
  gf_internal_t *h;
  struct gf_lazytable_data *ltd;
  unsigned long uls, uld;
  
  h = (gf_internal_t *) gf->scratch;

  uls = ((unsigned long) src) & 0xf;
  uld = ((unsigned long) dest) & 0xf;
  if (uls != uld || uls % 2 != 0 || bytes % 2 != 0) gf_alignment_error("gf_w16_table_lazy_multiply_region", 2);

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  ltd = (struct gf_lazytable_data *) h->private;

  ltd->lazytable[0] = 0;
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
   
  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

  if (uls != 0) {
    while (uls != 16 && d16 < top) {
      prod = (xor) ? *d16 : 0;
      prod ^= ltd->lazytable[*s16];
      *d16 = prod;
      s16++;
      d16++;
      uls += 2;
    }
    if (d16 == top) return;
  }

  uls = ((unsigned long) top) & 0xf;
  uld = ((unsigned long) top) ^ uls;
  top64 = (uint64_t *) uld;
  s64 = (uint64_t *) s16;
  d64 = (uint64_t *) d16;
  
  /* Unrolling doesn't seem to matter 
  while (d64 != top64) {
    a = *s64;

    prod = ltd->lazytable[a >> 48];
    a <<= 16;
    prod <<= 16;

    prod ^= ltd->lazytable[a >> 48];
    a <<= 16;
    prod <<= 16;

    prod ^= ltd->lazytable[a >> 48];
    a <<= 16;
    prod <<= 16;

    prod ^= ltd->lazytable[a >> 48];

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
      prod ^= ltd->lazytable[a >> 48];
      a <<= 16;
    }
    prod ^= ((xor) ? *d64 : 0); 
    *d64 = prod;
    *s64++;
    *d64++;
  } 
  
  if (uls != 0) {
    d16 = (uint16_t *) d64;
    s16 = (uint16_t *) s64;
    while (d16 < top) {
      prod = (xor) ? *d16 : 0;
      prod ^= ltd->lazytable[*s16];
      *d16 = prod;
      s16++;
      d16++;
    }
  }
  return;
}

static
void
gf_w16_split_4_16_lazy_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  uint64_t i, j, *s64, *d64, *top64;;
  uint64_t a, c, prod;
  uint16_t *s16, *d16, *top;
  uint8_t low[4][16];
  uint8_t high[4][16];
  unsigned long uls, uld;

  __m128i  mask, ta, tb, ti, tpl, tph, tlow[4], thigh[4], shuffler, unshuffler, tta, ttb;

  struct gf_single_table_data *std;

  uls = ((unsigned long) src) & 0xf;
  uld = ((unsigned long) dest) & 0xf;
  if (uls != uld || uls % 2 != 0 || bytes % 2 != 0) gf_alignment_error("gf_w16_split_4_16_lazy_sse_altmap_multiply_region", 2);

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      prod = gf_w16_log_multiply(gf, c, val);
      low[i][j] = (prod & 0xff);
      high[i][j] = (prod >> 8);
    }
  }

  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

  if (uls != 0) {
    while (uls != 16 && d16 < top) {
      a = *s16;
      prod = (xor) ? *d16 : 0;
      for (i = 0; i < 4; i++) {
        c = a & 0xf;
        prod ^= low[i][c];
        prod ^= (high[i][c] << 8);
        a >>= 4;
      }
      *d16 = prod;
      s16++;
      d16++;
      uls += 2;
    }
    if (d16 == top) return;
  }

  for (i = 0; i < 4; i++) {
    tlow[i] = _mm_loadu_si128((__m128i *)low[i]);
    thigh[i] = _mm_loadu_si128((__m128i *)high[i]);
  }

  uls = ((unsigned long) top);
  uld = ((unsigned long) d16);
  bytes = (uls - uld);
  if ((bytes & 0x1f) != 0) bytes -= (bytes & 0x1f);

  top64 = (uint64_t *) (uld + bytes);
  s64 = (uint64_t *) s16;
  d64 = (uint64_t *) d16;
  mask = _mm_set1_epi8 (0x0f);
  shuffler = _mm_set_epi8(0xf, 0xd, 0xb, 0x9, 7, 5, 3, 1, 0xe, 0xc, 0xa, 8, 6, 4, 2, 0);
  unshuffler = _mm_set_epi8(0xf, 7, 0xe, 6, 0xd, 5, 0xc, 4, 0xb, 3, 0xa, 2, 9, 1, 8, 0);

  if (xor) {
    while (d64 != top64) {
      
      ta = _mm_load_si128((__m128i *) s64);
      MM_PRINT("Ta", ta);
      tb = _mm_load_si128((__m128i *) (s64+2));
      MM_PRINT("Tb", tb);
      tta = _mm_shuffle_epi8(ta, shuffler);
      ttb = _mm_shuffle_epi8(tb, shuffler);
      ta = _mm_unpackhi_epi64(ttb, tta);
      MM_PRINT("New ta", ta);
      tb = _mm_unpacklo_epi64(ttb, tta);
      MM_PRINT("New tb", tb);
      exit(0);
      

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

      tta = _mm_unpackhi_epi64(tpl, tph);
      ttb = _mm_unpacklo_epi64(tpl, tph);
      ta = _mm_shuffle_epi8(tta, unshuffler);
      tb = _mm_shuffle_epi8(ttb, unshuffler);
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
      tta = _mm_shuffle_epi8(ta, shuffler);
      ttb = _mm_shuffle_epi8(tb, shuffler);
      ta = _mm_unpackhi_epi64(ttb, tta);
      tb = _mm_unpacklo_epi64(ttb, tta);

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

      tta = _mm_unpackhi_epi64(tpl, tph);
      ttb = _mm_unpacklo_epi64(tpl, tph);
      ta = _mm_shuffle_epi8(tta, unshuffler);
      tb = _mm_shuffle_epi8(ttb, unshuffler);
      _mm_store_si128 ((__m128i *)d64, ta);
      _mm_store_si128 ((__m128i *)(d64+2), tb);

      d64 += 4;
      s64 += 4;
    }
  }

  d16 = (uint16_t *) d64;
  s16 = (uint16_t *) s64;

  while (d16 != top) {
    a = *s16;
    prod = (xor) ? *d16 : 0;
    for (i = 0; i < 4; i++) {
      c = a & 0xf;
      prod ^= low[i][c];
      prod ^= (high[i][c] << 8);
      a >>= 4;
    }
    *d16 = prod;
    s16++;
    d16++;
  }
#endif
}

/*
static
void
gf_w16_split_4_16_lazy_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  uint64_t i, j, *s64, *d64, *top64;;
  uint64_t a, c, prod;
  uint16_t *s16, *d16, *top;
  uint8_t low[4][16];
  uint8_t high[4][16];
  unsigned long uls, uld;

  __m128i  mask, ta, ti, tp, tlow[4], thigh[4];

  struct gf_single_table_data *std;

  uls = ((unsigned long) src) & 0xf;
  uld = ((unsigned long) dest) & 0xf;
  if (uls != uld || uls % 2 != 0 || bytes % 2 != 0) gf_alignment_error("gf_w16_split_4_16_lazy_sse_multiply_region", 2);

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      prod = gf_w16_log_multiply(gf, c, val);
      low[i][j] = (prod & 0xff);
      high[i][j] = (prod >> 8);
    }
  }

  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

  if (uls != 0) {
    while (uls != 16 && d16 < top) {
      a = *s16;
      prod = (xor) ? *d16 : 0;
      for (i = 0; i < 4; i++) {
        c = a & 0xf;
        prod ^= low[i][c];
        prod ^= (high[i][c] << 8);
        a >>= 4;
      }
      *d16 = prod;
      s16++;
      d16++;
      uls += 2;
    }
    if (d16 == top) return;
  }

  for (i = 0; i < 4; i++) {
    tlow[i] = _mm_loadu_si128((__m128i *)low[i]);
    thigh[i] = _mm_loadu_si128((__m128i *)high[i]);
  }

  uls = ((unsigned long) top) & 0xf;
  uld = ((unsigned long) top) ^ uls;
  top64 = (uint64_t *) uld;
  s64 = (uint64_t *) s16;
  d64 = (uint64_t *) d16;
  mask = _mm_set1_epi16 (0x0f);

  if (xor) {
    while (d64 != top64) {
      ta = _mm_load_si128((__m128i *) s64);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_shuffle_epi8 (tlow[0], ti);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[0], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[1], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[2], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ti = _mm_srli_epi16(ta, 4);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[3], ti);
      tp = _mm_xor_si128 (tp, ti);
      ti = _mm_load_si128((__m128i *)d64);
      tp = _mm_xor_si128 (tp, ti);
      _mm_store_si128 ((__m128i *)d64, tp);
      s64 += 2;
      d64 += 2;
    }
  } else {
    while (d64 != top64) {
      ta = _mm_load_si128((__m128i *) s64);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_shuffle_epi8 (tlow[0], ti);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[0], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[1], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[1], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ta = _mm_srli_epi16(ta, 4);
      ti = _mm_and_si128 (mask, ta);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[2], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[2], ti);
      tp = _mm_xor_si128 (tp, ti);
  
      ti = _mm_srli_epi16(ta, 4);
      tp = _mm_xor_si128(_mm_shuffle_epi8 (tlow[3], ti), tp);
      ti = _mm_slli_epi16 (ti, 8);
      ti = _mm_shuffle_epi8 (thigh[3], ti);
      tp = _mm_xor_si128 (tp, ti);
      _mm_store_si128 ((__m128i *)d64, tp);
      s64 += 2;
      d64 += 2;
    }
  }

  d16 = (uint16_t *) d64;
  s16 = (uint16_t *) s64;

  while (d16 != top) {
    a = *s16;
    prod = (xor) ? *d16 : 0;
    for (i = 0; i < 4; i++) {
      c = a & 0xf;
      prod ^= low[i][c];
      prod ^= (high[i][c] << 8);
      a >>= 4;
    }
    *d16 = prod;
    s16++;
    d16++;
  }
#endif
}
*/


static
void
gf_w16_split_4_16_lazy_sse_altmap_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
#ifdef   INTEL_SSE4
  uint64_t i, j, *s64, *d64, *top64;;
  uint64_t a, c, prod;
  uint16_t *s16, *d16, *top;
  uint8_t low[4][16];
  uint8_t high[4][16];
  unsigned long uls, uld;

  __m128i  mask, ta, tb, ti, tpl, tph, tlow[4], thigh[4];

  struct gf_single_table_data *std;

  uls = ((unsigned long) src) & 0xf;
  uld = ((unsigned long) dest) & 0xf;
  if (uls != uld || uls % 2 != 0 || bytes % 2 != 0) gf_alignment_error("gf_w16_split_4_16_lazy_sse_altmap_multiply_region", 2);

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  for (j = 0; j < 16; j++) {
    for (i = 0; i < 4; i++) {
      c = (j << (i*4));
      prod = gf_w16_log_multiply(gf, c, val);
      low[i][j] = (prod & 0xff);
      high[i][j] = (prod >> 8);
    }
  }

  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;
  top = (uint16_t *) (dest+bytes);

  if (uls != 0) {
    while (uls != 16 && d16 < top) {
      a = *s16;
      prod = (xor) ? *d16 : 0;
      for (i = 0; i < 4; i++) {
        c = a & 0xf;
        prod ^= low[i][c];
        prod ^= (high[i][c] << 8);
        a >>= 4;
      }
      *d16 = prod;
      s16++;
      d16++;
      uls += 2;
    }
    if (d16 == top) return;
  }

  for (i = 0; i < 4; i++) {
    tlow[i] = _mm_loadu_si128((__m128i *)low[i]);
    thigh[i] = _mm_loadu_si128((__m128i *)high[i]);
  }

  uls = ((unsigned long) top);
  uld = ((unsigned long) d16);
  bytes = (uls - uld);
  if ((bytes & 0x1f) != 0) bytes -= (bytes & 0x1f);

  top64 = (uint64_t *) (uld + bytes);
  s64 = (uint64_t *) s16;
  d64 = (uint64_t *) d16;
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

  d16 = (uint16_t *) d64;
  s16 = (uint16_t *) s64;

  while (d16 != top) {
    a = *s16;
    prod = (xor) ? *d16 : 0;
    for (i = 0; i < 4; i++) {
      c = a & 0xf;
      prod ^= low[i][c];
      prod ^= (high[i][c] << 8);
      a >>= 4;
    }
    *d16 = prod;
    s16++;
    d16++;
  }
#endif
}

static 
int gf_w16_split_init(gf_t *gf)
{
  gf_internal_t *h;
  gf_w16_log_init(gf);

  h = (gf_internal_t *) gf->scratch;
  if (h->arg1 == 8 || h->arg2 == 8) {
    gf->multiply_region.w16 = gf_w16_split_8_16_lazy_multiply_region;
  } else if (h->arg1 == 4 || h->arg2 == 4) {
    if (h->region_type & GF_REGION_SSE) {
      if (h->region_type & GF_REGION_ALTMAP) {
        gf->multiply_region.w16 = gf_w16_split_4_16_lazy_sse_altmap_multiply_region;
      } else {
        gf->multiply_region.w16 = gf_w16_split_4_16_lazy_sse_multiply_region;
      }
    } else {
      gf->multiply_region.w16 = gf_w16_split_4_16_lazy_multiply_region;
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
  gf->multiply_region.w16 = NULL;
  gf->multiply_region.w16 = gf_w16_table_lazy_multiply_region; 
  return 1;
}

static
void
gf_w16_log_zero_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  int i;
  uint16_t lv, b, c;
  uint16_t *s16, *d16;
  int num_syms = bytes >> 1;
  int sym_divisible = bytes % 2;

  struct gf_zero_logtable_data *ltd;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w16_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w16_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  ltd = (struct gf_zero_logtable_data*) ((gf_internal_t *) gf->scratch)->private;
  s16 = (uint16_t *) src;
  d16 = (uint16_t *) dest;

  lv = ltd->log_tbl[val];

  if (xor) {
    for (i = 0; i < num_syms; i++) {
      d16[i] ^= ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]];
    }
  } else {
    for (i = 0; i < num_syms; i++) {
      d16[i] = ltd->antilog_tbl[lv + ltd->log_tbl[s16[i]]];
    }
  }
}

static
inline
gf_val_16_t
gf_w16_log_zero_multiply (gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  struct gf_zero_logtable_data *ltd;

  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return ltd->antilog_tbl[ltd->log_tbl[a] + ltd->log_tbl[b]];
}

static
inline
gf_val_16_t
gf_w16_log_zero_divide (gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  int log_sum = 0;
  struct gf_zero_logtable_data *ltd;

  if (a == 0 || b == 0) return 0;
  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;

  log_sum = ltd->log_tbl[a] - ltd->log_tbl[b] + (GF_MULT_GROUP_SIZE);
  return (ltd->antilog_tbl[log_sum]);
}

static
gf_val_16_t
gf_w16_log_zero_inverse (gf_t *gf, gf_val_16_t a)
{
  struct gf_zero_logtable_data *ltd;

  ltd = (struct gf_zero_logtable_data *) ((gf_internal_t *) gf->scratch)->private;
  return (ltd->inv_tbl[a]);
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
      ltd->log_tbl[b] = (gf_val_16_t)i;
      ltd->antilog_tbl[i] = (gf_val_16_t)b;
      ltd->antilog_tbl[i+GF_MULT_GROUP_SIZE] = (gf_val_16_t)b;
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

  gf->inverse.w16 = gf_w16_log_zero_inverse;
  gf->divide.w16 = gf_w16_log_zero_divide;
  gf->multiply.w16 = gf_w16_log_zero_multiply;
  gf->multiply_region.w16 = gf_w16_log_zero_multiply_region;
  return 1;
}

static
gf_val_16_t
gf_w16_composite_multiply(gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t b0 = b & 0x00ff;
  uint8_t b1 = (b & 0xff00) >> 8;
  uint8_t a0 = a & 0x00ff;
  uint8_t a1 = (a & 0xff00) >> 8;
  uint8_t a1b1;

  a1b1 = base_gf->multiply.w8(base_gf, a1, b1);

  return ((base_gf->multiply.w8(base_gf, a0, b0) ^ a1b1) | ((base_gf->multiply.w8(base_gf, a1, b0) ^ base_gf->multiply.w8(base_gf, a0, b1) ^ base_gf->multiply.w8(base_gf, a1b1, GF_S_GF_8_2)) << 8));
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
gf_val_16_t
gf_w16_composite_inverse(gf_t *gf, gf_val_16_t a)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint8_t a0 = a & 0x00ff;
  uint8_t a1 = (a & 0xff00) >> 8;
  uint8_t c0, c1, d, tmp;
  uint16_t c;
  uint8_t a0inv, a1inv;

  if (a0 == 0) {
    a1inv = base_gf->inverse.w8(base_gf, a1);
    c0 = base_gf->multiply.w8(base_gf, a1inv, GF_S_GF_8_2);
    c1 = a1inv;
  } else if (a1 == 0) {
    c0 = base_gf->inverse.w8(base_gf, a0);
    c1 = 0;
  } else {
    a1inv = base_gf->inverse.w8(base_gf, a1);
    a0inv = base_gf->inverse.w8(base_gf, a0);

    d = base_gf->multiply.w8(base_gf, a1, a0inv);

    tmp = (base_gf->multiply.w8(base_gf, a1, a0inv) ^ base_gf->multiply.w8(base_gf, a0, a1inv) ^ GF_S_GF_8_2);
    tmp = base_gf->inverse.w8(base_gf, tmp);

    d = base_gf->multiply.w8(base_gf, d, tmp);

    c0 = base_gf->multiply.w8(base_gf, (d^1), a0inv);
    c1 = base_gf->multiply.w8(base_gf, d, a1inv);
  }

  c = c0 | (c1 << 8);

  return c;
}

static
gf_val_16_t
gf_w16_composite_divide(gf_t *gf, gf_val_16_t a, gf_val_16_t b)
{
  gf_val_16_t binv;

  binv = gf_w16_composite_inverse(gf, b);

  return gf_w16_composite_multiply(gf, a, binv);
}

static
void
gf_w16_composite_multiply_region_table(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  int i=0;
  struct gf_w8_single_table_data * std;
  uint8_t b0 = val & 0x00ff;
  uint8_t b1 = (val & 0xff00) >> 8;
  uint16_t *s16 = (uint16_t *) src;
  uint16_t *d16 = (uint16_t *) dest;
  uint8_t a0, a1, a1b1;
  int num_syms = bytes >> 1;
  int sym_divisible = bytes % 2;

  struct gf_logtable_data *ltd;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w16_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w16_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  std = (struct gf_w8_single_table_data *) h->private;

  if (xor) {
    for (i = 0;i < num_syms; i++) {
      a0 = s16[i] & 0x00ff;
      a1 = (s16[i] & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      d16[i] ^= ((std->mult[a0][b0] ^ a1b1) | ((std->mult[a1][b0] ^ std->mult[a0][b1] ^ std->mult[a1b1][GF_S_GF_8_2]) << 8));

    }
  } else {
    for (i = 0;i < num_syms; i++) {
      a0 = s16[i] & 0x00ff;
      a1 = (s16[i] & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      d16[i] = ((std->mult[a0][b0] ^ a1b1) | ((std->mult[a1][b0] ^ std->mult[a0][b1] ^ std->mult[a1b1][GF_S_GF_8_2]) << 8));
    }
  }
}

static
void
gf_w16_composite_multiply_region(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  int i=0;
  struct gf_w8_single_table_data * std;
  uint8_t b0 = val & 0x00ff;
  uint8_t b1 = (val & 0xff00) >> 8;
  uint16_t *s16 = (uint16_t *) src;
  uint16_t *d16 = (uint16_t *) dest;
  uint8_t a0, a1, a1b1;
  int num_syms = bytes >> 1;
  int sym_divisible = bytes % 2;

  struct gf_logtable_data *ltd;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w16_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w16_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  std = (struct gf_w8_single_table_data *) h->private;

  if (xor) {
    for (i = 0;i < num_syms; i++) {
      a0 = s16[i] & 0x00ff;
      a1 = (s16[i] & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      d16[i] ^= ((base_gf->multiply.w8(base_gf, a0, b0) ^ a1b1) |
                ((base_gf->multiply.w8(base_gf, a1, b0) ^ base_gf->multiply.w8(base_gf, a0, b1) ^ base_gf->multiply.w8(base_gf, a1b1, GF_S_GF_8_2)) << 8));

    }
  } else {
    for (i = 0;i < num_syms; i++) {
      a0 = s16[i] & 0x00ff;
      a1 = (s16[i] & 0xff00) >> 8;
      a1b1 = std->mult[a1][b1];

      d16[i] = ((base_gf->multiply.w8(base_gf, a0, b0) ^ a1b1) |
               ((base_gf->multiply.w8(base_gf, a1, b0) ^ base_gf->multiply.w8(base_gf, a0, b1) ^ base_gf->multiply.w8(base_gf, a1b1, GF_S_GF_8_2)) << 8));
    }
  }
}

static
void
gf_w16_composite_multiply_region_alt(gf_t *gf, void *src, void *dest, gf_val_16_t val, int bytes, int xor)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_val_8_t val0 = val & 0x00ff;
  gf_val_8_t val1 = (val & 0xff00) >> 8;
  int sub_reg_size = bytes / 2;

  if (!xor) {
    memset(dest, 0, bytes);
  }

  if (bytes % 2 != 0) gf_alignment_error("gf_w8_composite_multiply_region_alt", 1);

  base_gf->multiply_region.w8(base_gf, src, dest, val0, sub_reg_size, xor);
  base_gf->multiply_region.w8(base_gf, src+sub_reg_size, dest, val1, sub_reg_size, 1);
  base_gf->multiply_region.w8(base_gf, src, dest+sub_reg_size, val1, sub_reg_size, xor);
  base_gf->multiply_region.w8(base_gf, src+sub_reg_size, dest+sub_reg_size, val0, sub_reg_size, 1);
  base_gf->multiply_region.w8(base_gf, src+sub_reg_size, dest+sub_reg_size, base_gf->multiply.w8(base_gf, GF_S_GF_8_2, val1), sub_reg_size, 1);
}

static
int gf_w16_composite_init(gf_t *gf)
{
  struct gf_w8_single_table_data * std;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_val_16_t a, b;

  std = (struct gf_w8_single_table_data *) h->private;

  for (a = 0; a < 256; a++) {
    for (b = 0; b < 256; b++) {
      std->mult[a][b] = base_gf->multiply.w8(base_gf, a, b);
    }
  }

  if (h->region_type & GF_REGION_ALTMAP) {
    gf->multiply_region.w16 = gf_w16_composite_multiply_region_alt;
  } else {
    if (h->region_type & GF_REGION_SINGLE_TABLE) {
      gf->multiply_region.w16 = gf_w16_composite_multiply_region_table;
    } else {
      gf->multiply_region.w16 = gf_w16_composite_multiply_region;
    }
  }

  gf->multiply.w16 = gf_w16_composite_multiply;
  gf->divide.w16 = gf_w16_composite_divide;
  gf->inverse.w16 = gf_w16_composite_inverse;

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
    case GF_MULT_DEFAULT:
    case GF_MULT_SHIFT:
      if (arg1 != 0 || arg2 != 0 || region_type != 0) return -1;
      return sizeof(gf_internal_t);
      break;
    case GF_MULT_COMPOSITE:
      if (region_type & ~(GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP | GF_REGION_STDMAP)) return -1;
      if ((region_type & (GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP)) == (GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP)) return -1;
      if (arg1 == 2 && arg2 == 8) {
        return sizeof(gf_internal_t) + sizeof(struct gf_w8_single_table_data) + 64;
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

  gf->multiply.w16 = NULL;
  gf->divide.w16 = NULL;
  gf->inverse.w16 = NULL;
  gf->multiply_region.w16 = NULL;

  switch(h->mult_type) {
    case GF_MULT_LOG_TABLE:        
      if (h->arg1 == 1) {
        return gf_w16_log_zero_init(gf);
      } else {
        return gf_w16_log_init(gf);
      }
    case GF_MULT_SPLIT_TABLE: return gf_w16_split_init(gf);
    case GF_MULT_TABLE:       return gf_w16_table_init(gf);
    case GF_MULT_DEFAULT: 
    case GF_MULT_SHIFT:     if (gf_w16_shift_init(gf) == 0) return 0; break;
    case GF_MULT_COMPOSITE: if (gf_w16_composite_init(gf) == 0) return 0; break;
    default: return 0;
  }
  if (h->divide_type == GF_DIVIDE_EUCLID) {
    gf->divide.w16 = gf_w16_divide_from_inverse;
    gf->inverse.w16 = gf_w16_euclid;
  } else if (h->divide_type == GF_DIVIDE_MATRIX) {
    gf->divide.w16 = gf_w16_divide_from_inverse;
    gf->inverse.w16 = gf_w16_matrix;
  }

  if (gf->inverse.w16 != NULL && gf->divide.w16 == NULL) {
    gf->divide.w16 = gf_w16_divide_from_inverse;
  }
  if (gf->inverse.w16 == NULL && gf->divide.w16 != NULL) {
    gf->inverse.w16 = gf_w16_inverse_from_divide;
  }
  return 1;
}

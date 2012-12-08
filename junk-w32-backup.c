/*
 * gf_w32.c
 *
 * Routines for 32-bit Galois fields
 */

#define MM_PRINT32(s, r) { uint8_t blah[16], ii; printf("%-12s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 4) printf(" %02x%02x%02x%02x", blah[15-ii], blah[14-ii], blah[13-ii], blah[12-ii]); printf("\n"); }

#define MM_PRINT8(s, r) { uint8_t blah[16], ii; printf("%-12s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 1) printf("%s%02x", (ii%4==0) ? "   " : " ", blah[15-ii]); printf("\n"); }

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

#define GF_FIELD_WIDTH (32)
#define GF_FIRST_BIT (1 << 31)

#define GF_BASE_FIELD_WIDTH (16)
#define GF_BASE_FIELD_SIZE       (1 << GF_BASE_FIELD_WIDTH)
#define GF_BASE_FIELD_GROUP_SIZE  GF_BASE_FIELD_SIZE-1
#define GF_S_GF_16_2 (40188)
#define GF_MULTBY_TWO(p) (((p) & GF_FIRST_BIT) ? (((p) << 1) ^ h->prim_poly) : (p) << 1);


struct gf_w16_logtable_data {
    int              log_tbl[GF_BASE_FIELD_SIZE];
    gf_val_16_t      _antilog_tbl[GF_BASE_FIELD_SIZE * 4];
    gf_val_16_t      *antilog_tbl;
    gf_val_16_t      inv_tbl[GF_BASE_FIELD_SIZE];
};

struct gf_split_2_32_lazy_data {
    gf_val_32_t      last_value;
    gf_val_32_t      tables[16][4];
};

struct gf_split_8_8_data {
    gf_val_32_t      tables[7][256][256];
};

struct gf_split_4_32_lazy_data {
    gf_val_32_t      last_value;
    gf_val_32_t      tables[8][16];
};

static
inline
gf_val_32_t gf_w32_inverse_from_divide (gf_t *gf, gf_val_32_t a)
{
  return gf->divide.w32(gf, 1, a);
}

static
inline
gf_val_32_t gf_w32_divide_from_inverse (gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  b = gf->inverse.w32(gf, b);
  return gf->multiply.w32(gf, a, b);
}

static
void
gf_w32_multiply_region_from_single(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int 
xor)
{
  int i;
  gf_val_32_t *s32;
  gf_val_32_t *d32;
   
  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest; 
 
  if (xor) {
    for (i = 0; i < bytes/sizeof(gf_val_32_t); i++) {
      d32[i] ^= gf->multiply.w32(gf, val, s32[i]);
    } 
  } else {
    for (i = 0; i < bytes/sizeof(gf_val_32_t); i++) {
      d32[i] = gf->multiply.w32(gf, val, s32[i]);
    } 
  }
}

static
inline
gf_val_32_t gf_w32_euclid (gf_t *gf, gf_val_32_t b)
{
  gf_val_32_t e_i, e_im1, e_ip1;
  gf_val_32_t d_i, d_im1, d_ip1;
  gf_val_32_t y_i, y_im1, y_ip1;
  gf_val_32_t c_i;

  if (b == 0) return -1;
  e_im1 = ((gf_internal_t *) (gf->scratch))->prim_poly;
  e_i = b;
  d_im1 = 32;
  for (d_i = d_im1-1; ((1 << d_i) & e_i) == 0; d_i--) ;
  y_i = 1;
  y_im1 = 0;

  while (e_i != 1) {

    e_ip1 = e_im1;
    d_ip1 = d_im1;
    c_i = 0;

    while (d_ip1 >= d_i) {
      c_i ^= (1 << (d_ip1 - d_i));
      e_ip1 ^= (e_i << (d_ip1 - d_i));
      d_ip1--;
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
inline
gf_val_32_t gf_w32_matrix (gf_t *gf, gf_val_32_t b)
{
  return gf_bitmatrix_inverse(b, 32, ((gf_internal_t *) (gf->scratch))->prim_poly);
}

/* JSP: GF_MULT_SHIFT: The world's dumbest multiplication algorithm.  I only
   include it for completeness.  It does have the feature that it requires no
   extra memory.  
*/

static
inline
gf_val_32_t
gf_w32_shift_multiply (gf_t *gf, gf_val_32_t a32, gf_val_32_t b32)
{
  uint64_t product, i, pp, a, b, one;
  gf_internal_t *h;
  
  a = a32;
  b = b32;
  h = (gf_internal_t *) gf->scratch;
  one = 1;
  pp = h->prim_poly | (one << 32);

  product = 0;

  for (i = 0; i < GF_FIELD_WIDTH; i++) { 
    if (a & (one << i)) product ^= (b << i);
  }
  for (i = (GF_FIELD_WIDTH*2-1); i >= GF_FIELD_WIDTH; i--) {
    if (product & (one << i)) product ^= (pp << (i-GF_FIELD_WIDTH)); 
  }
  return product;
}

static 
int gf_w32_shift_init(gf_t *gf)
{
  gf->multiply.w32 = gf_w32_shift_multiply;
  gf->inverse.w32 = gf_w32_euclid;
  gf->multiply_region.w32 = gf_w32_multiply_region_from_single;
  return 1;
}

static
inline
gf_val_32_t
gf_w32_split_8_8_multiply (gf_t *gf, gf_val_32_t a32, gf_val_32_t b32)
{
  uint32_t product, i, j, mask, tb;
  gf_internal_t *h;
  struct gf_split_8_8_data *d8;
  
  h = (gf_internal_t *) gf->scratch;
  d8 = (struct gf_split_8_8_data *) h->private;
  product = 0;
  mask = 0xff;

  for (i = 0; i < 4; i++) {
    tb = b32;
    for (j = 0; j < 4; j++) {
      product ^= d8->tables[i+j][a32&mask][tb&mask];
      tb >>= 8;
    }
    a32 >>= 8;
  }
  return product;
}

static
inline
void
gf_w32_split_8_8_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  uint32_t product, mask, tb, tv, tp;
  gf_internal_t *h;
  struct gf_split_8_8_data *d8;
  uint32_t *p00, *p01, *p02, *p03;
  uint32_t *p10, *p11, *p12, *p13;
  uint32_t *p20, *p21, *p22, *p23;
  uint32_t *p30, *p31, *p32, *p33;
  uint32_t *s32, *d32, *top;
  unsigned long uls, uld;
  
  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0x7) != (uld & 0x7))) gf_alignment_error("gf_w32_split_8_8_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_8_8_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  tv = val;
  h = (gf_internal_t *) gf->scratch;
  d8 = (struct gf_split_8_8_data *) h->private;
  mask = 0xff;

  p00 = &(d8->tables[0][val&mask][0]);
  p01 = &(d8->tables[1][val&mask][0]);
  p02 = &(d8->tables[2][val&mask][0]);
  p03 = &(d8->tables[3][val&mask][0]);
  val >>= 8;
  p10 = &(d8->tables[1][val&mask][0]);
  p11 = &(d8->tables[2][val&mask][0]);
  p12 = &(d8->tables[3][val&mask][0]);
  p13 = &(d8->tables[4][val&mask][0]);
  val >>= 8;
  p20 = &(d8->tables[2][val&mask][0]);
  p21 = &(d8->tables[3][val&mask][0]);
  p22 = &(d8->tables[4][val&mask][0]);
  p23 = &(d8->tables[5][val&mask][0]);
  val >>= 8;
  p30 = &(d8->tables[3][val&mask][0]);
  p31 = &(d8->tables[4][val&mask][0]);
  p32 = &(d8->tables[5][val&mask][0]);
  p33 = &(d8->tables[6][val&mask][0]);

  s32 = (uint32_t *) src;
  d32 = (uint32_t *) dest;
  top = (d32 + (bytes/4));

  while (d32 < top) {
    tb = *s32;
    tp = *d32;
    product = (xor) ? (*d32) : 0;
    product ^= p00[tb&mask];
    product ^= p10[tb&mask];
    product ^= p20[tb&mask];
    product ^= p30[tb&mask];

    tb >>= 8;
    product ^= p01[tb&mask];
    product ^= p11[tb&mask];
    product ^= p21[tb&mask];
    product ^= p31[tb&mask];

    tb >>= 8;
    product ^= p02[tb&mask];
    product ^= p12[tb&mask];
    product ^= p22[tb&mask];
    product ^= p32[tb&mask];

    tb >>= 8;
    product ^= p03[tb&mask];
    product ^= p13[tb&mask];
    product ^= p23[tb&mask];
    product ^= p33[tb&mask];
    *d32 = product;
    s32++;
    d32++;
  }
}

static
void
gf_w32_split_2_32_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h;
  struct gf_split_2_32_lazy_data *ld;
  int i;
  gf_val_32_t pp, v, v2, s, *s32, *d32, *top;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0x7) != (uld & 0x7))) gf_alignment_error("gf_w32_split_2_32_lazy_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_2_32_lazy_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  ld = (struct gf_split_2_32_lazy_data *) h->private;
  
  if (ld->last_value != val) {
    v = val;
    for (i = 0; i < 16; i++) {
      v2 = (v << 1);
      if (v & GF_FIRST_BIT) v2 ^= pp;
      ld->tables[i][0] = 0;
      ld->tables[i][1] = v;
      ld->tables[i][2] = v2;
      ld->tables[i][3] = (v2 ^ v);
      v = (v2 << 1);
      if (v2 & GF_FIRST_BIT) v ^= pp;
    }
  }
  ld->last_value = val;

  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest;
  top = d32 + (bytes/4);

  while (d32 != top) {
    v = (xor) ? *d32 : 0;
    s = *s32;
    i = 0;
    while (s != 0) {
      v ^= ld->tables[i][s&3];
      s >>= 2;
      i++;
    }
    *d32 = v;
    d32++;
    s32++;
  }
}

static
void
gf_w32_split_2_32_lazy_sse_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef INTEL_SSE4
  unsigned long uls, uld;
  gf_internal_t *h;
  int i, m, j, tindex;
  gf_val_32_t pp, v, v2, s, *s32, *d32, *top;
  __m128i vi, si, pi, shuffler, tables[16], adder, xi, mask1, mask2;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0xf) != (uld & 0xf))) gf_alignment_error("gf_w32_split_2_32_lazy_sse_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_2_32_lazy_sse_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;
  
  uls &= 0xf;

  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest;
  top = d32 + (bytes/4);
  
  if (uls != 0) {
    while (uls != 16) {
      if (xor) {
        *d32 ^= gf->multiply.w32(gf, *s32, val);
      } else {
        *d32 = gf->multiply.w32(gf, *s32, val);
      }
      *s32++;
      *d32++;
      if (d32 == top) return;
      uls += 4;
    }
  }

  uld = (unsigned long) top;
  top = (gf_val_32_t *) (uld - (uld & 0xf));
  uld &= 0xf;
  
  v = val;
  for (i = 0; i < 16; i++) {
    v2 = (v << 1);
    if (v & GF_FIRST_BIT) v2 ^= pp;
    tables[i] = _mm_set_epi32(v2 ^ v, v2, v, 0);
    v = (v2 << 1);
    if (v2 & GF_FIRST_BIT) v ^= pp;
  }

  shuffler = _mm_set_epi8(0xc, 0xc, 0xc, 0xc, 8, 8, 8, 8, 4, 4, 4, 4, 0, 0, 0, 0);
  adder = _mm_set_epi8(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
  mask1 = _mm_set1_epi8(0x3);
  mask2 = _mm_set1_epi8(0xc);

  while (d32 != top) {
    pi = (xor) ? _mm_load_si128 ((__m128i *) d32) : _mm_setzero_si128();
    vi = _mm_load_si128((__m128i *) s32);
 
    tindex = 0;
    for (i = 0; i < 4; i++) {
      si = _mm_shuffle_epi8(vi, shuffler);

      xi = _mm_and_si128(si, mask1);
      xi = _mm_slli_epi16(xi, 2);
      xi = _mm_xor_si128(xi, adder);
      pi = _mm_xor_si128(pi, _mm_shuffle_epi8(tables[tindex], xi));
      tindex++;

      xi = _mm_and_si128(si, mask2);
      xi = _mm_xor_si128(xi, adder);
      pi = _mm_xor_si128(pi, _mm_shuffle_epi8(tables[tindex], xi));
      si = _mm_srli_epi16(si, 2);
      tindex++;

      xi = _mm_and_si128(si, mask2);
      xi = _mm_xor_si128(xi, adder);
      pi = _mm_xor_si128(pi, _mm_shuffle_epi8(tables[tindex], xi));
      si = _mm_srli_epi16(si, 2);
      tindex++;

      xi = _mm_and_si128(si, mask2);
      xi = _mm_xor_si128(xi, adder);
      pi = _mm_xor_si128(pi, _mm_shuffle_epi8(tables[tindex], xi));
      si = _mm_srli_epi16(si, 2);
      tindex++;
      
      vi = _mm_srli_epi32(vi, 8);
    }
    _mm_store_si128((__m128i *) d32, pi);
    d32 += 4;
    s32 += 4;
  }

  while (uld > 0) {
    if (xor) {
      *d32 ^= gf->multiply.w32(gf, *s32, val);
    } else {
      *d32 = gf->multiply.w32(gf, *s32, val);
    }
    *s32++;
    *d32++;
    uld -= 4;
  }


#endif
}

static
void
gf_w32_split_4_32_lazy_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h;
  struct gf_split_4_32_lazy_data *ld;
  int i, j, k;
  gf_val_32_t pp, v, s, *s32, *d32, *top;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0x7) != (uld & 0x7))) gf_alignment_error("gf_w32_split_4_32_lazy_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_4_32_lazy_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  ld = (struct gf_split_4_32_lazy_data *) h->private;
  
  if (ld->last_value != val) {
    v = val;
    for (i = 0; i < 8; i++) {
      ld->tables[i][0] = 0;
      for (j = 1; j < 16; j <<= 1) {
        for (k = 0; k < j; k++) {
          ld->tables[i][k^j] = (v ^ ld->tables[i][k]);
        }
        v = (v & GF_FIRST_BIT) ? ((v << 1) ^ pp) : (v << 1);
      }
    }
  }
  ld->last_value = val;

  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest;
  top = d32 + (bytes/4);

  while (d32 != top) {
    v = (xor) ? *d32 : 0;
    s = *s32;
    i = 0;
    while (s != 0) {
      v ^= ld->tables[i][s&0xf];
      s >>= 4;
      i++;
    }
    *d32 = v;
    d32++;
    s32++;
  }
}

static
void
gf_w32_split_4_32_lazy_sse_altmap_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef INTEL_SSE4
  unsigned long uls, uld;
  gf_internal_t *h;
  int i, m, j, k, tindex;
  gf_val_32_t pp, v, s, *s32, *d32, *top, *realtop;
  __m128i vi, si, tables[8][4], p0, p1, p2, p3, mask1, v0, v1, v2, v3;
  __m128i tv1, tv2, tv3, tv0;
  struct gf_split_4_32_lazy_data *ld;
  uint8_t btable[16];

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0xf) != (uld & 0xf))) gf_alignment_error("gf_w32_split_4_32_lazy_sse_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_4_32_lazy_sse_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;
  
  uls &= 0xf;

  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest;
  top = d32 + (bytes/4);
  
  if (uls != 0) {
    while (uls != 16) {
      if (xor) {
        *d32 ^= gf->multiply.w32(gf, *s32, val);
      } else {
        *d32 = gf->multiply.w32(gf, *s32, val);
      }
      *s32++;
      *d32++;
      if (d32 == top) return;
      uls += 4;
    }
  }

  uld = (unsigned long) top;
  realtop = top;
  
  /* You need the size of this region to be a multiple of 64 bytes */
  bytes = (top - d32);
  bytes -= (bytes & 0xf);
  top = (d32 + bytes);

  ld = (struct gf_split_4_32_lazy_data *) h->private;
 
  v = val;
  for (i = 0; i < 8; i++) {
    ld->tables[i][0] = 0;
    for (j = 1; j < 16; j <<= 1) {
      for (k = 0; k < j; k++) {
        ld->tables[i][k^j] = (v ^ ld->tables[i][k]);
      }
      v = (v & GF_FIRST_BIT) ? ((v << 1) ^ pp) : (v << 1);
    }
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 16; k++) {
        btable[k] = (uint8_t) ld->tables[i][k];
        ld->tables[i][k] >>= 8;
      }
      tables[i][j] = _mm_loadu_si128((__m128i *) btable);
    }
  }

  mask1 = _mm_set1_epi8(0xf);

  if (xor) {
    while (d32 != top) {
      p0 = _mm_load_si128 ((__m128i *) d32);
      p1 = _mm_load_si128 ((__m128i *) (d32+4));
      p2 = _mm_load_si128 ((__m128i *) (d32+8));
      p3 = _mm_load_si128 ((__m128i *) (d32+12));
  
      v0 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v1 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v2 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v3 = _mm_load_si128((__m128i *) s32); s32 += 4;
  
      si = _mm_and_si128(v0, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[0][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[0][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[0][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[0][3], si));
      
      v0 = _mm_srli_epi32(v0, 4);
      si = _mm_and_si128(v0, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[1][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[1][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[1][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[1][3], si));
  
      si = _mm_and_si128(v1, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[2][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[2][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[2][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[2][3], si));
      
      v1 = _mm_srli_epi32(v1, 4);
      si = _mm_and_si128(v1, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[3][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[3][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[3][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[3][3], si));
  
      si = _mm_and_si128(v2, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[4][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[4][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[4][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[4][3], si));
      
      v2 = _mm_srli_epi32(v2, 4);
      si = _mm_and_si128(v2, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[5][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[5][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[5][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[5][3], si));
  
      si = _mm_and_si128(v3, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[6][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[6][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[6][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[6][3], si));
      
      v3 = _mm_srli_epi32(v3, 4);
      si = _mm_and_si128(v3, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[7][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[7][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[7][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[7][3], si));
  
      _mm_store_si128((__m128i *) d32, p0);
      _mm_store_si128((__m128i *) (d32+4), p1);
      _mm_store_si128((__m128i *) (d32+8), p2);
      _mm_store_si128((__m128i *) (d32+12), p3);
      d32 += 16;
    } 
  } else {
    while (d32 != top) {
  
      v0 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v1 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v2 = _mm_load_si128((__m128i *) s32); s32 += 4;
      v3 = _mm_load_si128((__m128i *) s32); s32 += 4;

     
  
      si = _mm_and_si128(v0, mask1);
      p0 = _mm_shuffle_epi8(tables[0][0], si);
      p1 = _mm_shuffle_epi8(tables[0][1], si);
      p2 = _mm_shuffle_epi8(tables[0][2], si);
      p3 = _mm_shuffle_epi8(tables[0][3], si);
      
      v0 = _mm_srli_epi32(v0, 4);
      si = _mm_and_si128(v0, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[1][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[1][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[1][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[1][3], si));
  
      si = _mm_and_si128(v1, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[2][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[2][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[2][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[2][3], si));
      
      v1 = _mm_srli_epi32(v1, 4);
      si = _mm_and_si128(v1, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[3][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[3][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[3][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[3][3], si));
  
      si = _mm_and_si128(v2, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[4][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[4][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[4][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[4][3], si));
      
      v2 = _mm_srli_epi32(v2, 4);
      si = _mm_and_si128(v2, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[5][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[5][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[5][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[5][3], si));
  
      si = _mm_and_si128(v3, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[6][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[6][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[6][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[6][3], si));
      
      v3 = _mm_srli_epi32(v3, 4);
      si = _mm_and_si128(v3, mask1);
      p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[7][0], si));
      p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[7][1], si));
      p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[7][2], si));
      p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[7][3], si));
  
      _mm_store_si128((__m128i *) d32, p0);
      _mm_store_si128((__m128i *) (d32+4), p1);
      _mm_store_si128((__m128i *) (d32+8), p2);
      _mm_store_si128((__m128i *) (d32+12), p3);
      d32 += 16;
    } 
  }

  while (d32 < realtop) {
    if (xor) {
      *d32 ^= gf->multiply.w32(gf, *s32, val);
    } else {
      *d32 = gf->multiply.w32(gf, *s32, val);
    }
    *s32++;
    *d32++;
  }


#endif
}

/*
static
void
gf_w32_split_4_32_lazy_sse_altmap_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
#ifdef INTEL_SSE4
  unsigned long uls, uld;
  gf_internal_t *h;
  int i, m, j, k, tindex;
  gf_val_32_t pp, v, s, *s32, *d32, *top, *realtop;
  __m128i vi, si, tables[8][4], p0, p1, p2, p3, mask1;
  struct gf_split_4_32_lazy_data *ld;
  uint8_t btable[16];

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if (uls %4 != 0 || ((uls & 0xf) != (uld & 0xf))) gf_alignment_error("gf_w32_split_4_32_lazy_sse_multiply_region", 4);
  if (bytes % 4 != 0) {
    gf_alignment_error("gf_w32_split_4_32_lazy_sse_multiply_region: buffer size not divisible by symbol size = 4 bytes", 4);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;
  
  uls &= 0xf;

  s32 = (gf_val_32_t *) src;
  d32 = (gf_val_32_t *) dest;
  top = d32 + (bytes/4);
  
  if (uls != 0) {
    while (uls != 16) {
      if (xor) {
        *d32 ^= gf->multiply.w32(gf, *s32, val);
      } else {
        *d32 = gf->multiply.w32(gf, *s32, val);
      }
      *s32++;
      *d32++;
      if (d32 == top) return;
      uls += 4;
    }
  }

  uld = (unsigned long) top;
  realtop = top;
  
  bytes = (top - d32);
  bytes -= (bytes & 0xf);
  top = (d32 + bytes);

  ld = (struct gf_split_4_32_lazy_data *) h->private;
 
  v = val;
  for (i = 0; i < 8; i++) {
    ld->tables[i][0] = 0;
    for (j = 1; j < 16; j <<= 1) {
      for (k = 0; k < j; k++) {
        ld->tables[i][k^j] = (v ^ ld->tables[i][k]);
      }
      v = (v & GF_FIRST_BIT) ? ((v << 1) ^ pp) : (v << 1);
    }
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 16; k++) {
        btable[k] = (uint8_t) ld->tables[i][k];
        ld->tables[i][k] >>= 8;
      }
      tables[i][j] = _mm_loadu_si128((__m128i *) btable);
    }
  }

  mask1 = _mm_set1_epi8(0xf);

  if (xor) {
    while (d32 != top) {
        p0 = _mm_load_si128 ((__m128i *) d32);
        p1 = _mm_load_si128 ((__m128i *) (d32+4));
        p2 = _mm_load_si128 ((__m128i *) (d32+8));
        p3 = _mm_load_si128 ((__m128i *) (d32+12));
  
      for (i = 0; i < 8; i++) {
        vi = _mm_load_si128((__m128i *) s32);
  
        si = _mm_and_si128(vi, mask1);
        p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[i][0], si));
        p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[i][1], si));
        p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[i][2], si));
        p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[i][3], si));
  
        i++;
        vi = _mm_srli_epi32(vi, 4);
        si = _mm_and_si128(vi, mask1);
        p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[i][0], si));
        p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[i][1], si));
        p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[i][2], si));
        p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[i][3], si));
        s32 += 4;
      }
      _mm_store_si128((__m128i *) d32, p0);
      _mm_store_si128((__m128i *) (d32+4), p1);
      _mm_store_si128((__m128i *) (d32+8), p2);
      _mm_store_si128((__m128i *) (d32+12), p3);
      d32 += 16;
    } 
  } else {
    while (d32 != top) {
      for (i = 0; i < 8; i++) {
        vi = _mm_load_si128((__m128i *) s32);
  
        si = _mm_and_si128(vi, mask1);
        p0 = _mm_shuffle_epi8(tables[i][0], si);
        p1 = _mm_shuffle_epi8(tables[i][1], si);
        p2 = _mm_shuffle_epi8(tables[i][2], si);
        p3 = _mm_shuffle_epi8(tables[i][3], si);
  
        i++;
        vi = _mm_srli_epi32(vi, 4);
        si = _mm_and_si128(vi, mask1);
        p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(tables[i][0], si));
        p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(tables[i][1], si));
        p2 = _mm_xor_si128(p2, _mm_shuffle_epi8(tables[i][2], si));
        p3 = _mm_xor_si128(p3, _mm_shuffle_epi8(tables[i][3], si));
        s32 += 4;
      }
      _mm_store_si128((__m128i *) d32, p0);
      _mm_store_si128((__m128i *) (d32+4), p1);
      _mm_store_si128((__m128i *) (d32+8), p2);
      _mm_store_si128((__m128i *) (d32+12), p3);
      d32 += 16;
    } 
  }

  while (d32 < realtop) {
    if (xor) {
      *d32 ^= gf->multiply.w32(gf, *s32, val);
    } else {
      *d32 = gf->multiply.w32(gf, *s32, val);
    }
    *s32++;
    *d32++;
  }


#endif
}
*/

static 
int gf_w32_split_init(gf_t *gf)
{
  gf_internal_t *h;
  struct gf_split_2_32_lazy_data *ld2;
  struct gf_split_4_32_lazy_data *ld4;
  struct gf_split_8_8_data *d8;
  uint32_t p, basep;
  int i, j, exp;

  h = (gf_internal_t *) gf->scratch;

  /* Defaults */
  gf->multiply_region.w32 = gf_w32_multiply_region_from_single;
  gf->multiply.w32 = gf_w32_shift_multiply;
  gf->inverse.w32 = gf_w32_euclid;

  if (h->arg1 == 8 && h->arg2 == 8) {
    gf->multiply.w32 = gf_w32_split_8_8_multiply;
    gf->multiply_region.w32 = gf_w32_split_8_8_multiply_region;
    d8 = (struct gf_split_8_8_data *) h->private;
    basep = 1;
    for (exp = 0; exp < 7; exp++) {
      for (j = 0; j < 256; j++) d8->tables[exp][0][j] = 0;
      for (i = 0; i < 256; i++) d8->tables[exp][i][0] = 0;
      d8->tables[exp][1][1] = basep;
      for (i = 2; i < 256; i++) {
        if (i&1) {
          p = d8->tables[exp][i^1][1];
          d8->tables[exp][i][1] = p ^ basep;
        } else {
          p = d8->tables[exp][i>>1][1];
          d8->tables[exp][i][1] = GF_MULTBY_TWO(p);
        }
      }
      for (i = 1; i < 256; i++) {
        p = d8->tables[exp][i][1];
        for (j = 1; j < 256; j++) {
          if (j&1) {
            d8->tables[exp][i][j] = d8->tables[exp][i][j^1] ^ p;
          } else {
            d8->tables[exp][i][j] = GF_MULTBY_TWO(d8->tables[exp][i][j>>1]);
          }
        }
      }
      for (i = 0; i < 8; i++) basep = GF_MULTBY_TWO(basep);
    }
  }
  if ((h->arg1 == 2 && h->arg2 == 32) || (h->arg1 == 32 && h->arg2 == 2)) {
    ld2 = (struct gf_split_2_32_lazy_data *) h->private;
    ld2->last_value = 0;
    if (h->region_type & GF_REGION_SSE) {
      gf->multiply_region.w32 = gf_w32_split_2_32_lazy_sse_multiply_region;
    } else {
      gf->multiply_region.w32 = gf_w32_split_2_32_lazy_multiply_region;
    }
  } 
  if ((h->arg1 == 4 && h->arg2 == 32) || (h->arg1 == 32 && h->arg2 == 4)) {
    ld4 = (struct gf_split_4_32_lazy_data *) h->private;
    ld4->last_value = 0;
    if (h->region_type & GF_REGION_SSE) {
      if (h->region_type & GF_REGION_ALTMAP) {
        gf->multiply_region.w32 = gf_w32_split_4_32_lazy_sse_altmap_multiply_region;
      }
    } else {
      gf->multiply_region.w32 = gf_w32_split_4_32_lazy_multiply_region;
    }
  } 
  return 1;
}

static
gf_val_32_t
gf_w32_composite_multiply(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint16_t b0 = b & 0x0000ffff;
  uint16_t b1 = (b & 0xffff0000) >> 16;
  uint16_t a0 = a & 0x0000ffff;
  uint16_t a1 = (a & 0xffff0000) >> 16;
  uint16_t a1b1;

  a1b1 = base_gf->multiply.w16(base_gf, a1, b1);

  return ((base_gf->multiply.w16(base_gf, a0, b0) ^ a1b1) | ((base_gf->multiply.w16(base_gf, a1, b0) ^ base_gf->multiply.w16(base_gf, a0, b1) ^ base_gf->multiply.w16(base_gf, a1b1, GF_S_GF_16_2)) << 16));
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
gf_w32_composite_inverse(gf_t *gf, gf_val_32_t a)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  uint16_t a0 = a & 0x0000ffff;
  uint16_t a1 = (a & 0xffff0000) >> 16;
  uint16_t c0, c1, d, tmp;
  uint32_t c;
  uint16_t a0inv, a1inv;

  if (a0 == 0) {
    a1inv = base_gf->inverse.w16(base_gf, a1);
    c0 = base_gf->multiply.w16(base_gf, a1inv, GF_S_GF_16_2);
    c1 = a1inv;
  } else if (a1 == 0) {
    c0 = base_gf->inverse.w16(base_gf, a0);
    c1 = 0;
  } else {
    a1inv = base_gf->inverse.w16(base_gf, a1);
    a0inv = base_gf->inverse.w16(base_gf, a0);

    d = base_gf->multiply.w16(base_gf, a1, a0inv);

    tmp = (base_gf->multiply.w16(base_gf, a1, a0inv) ^ base_gf->multiply.w16(base_gf, a0, a1inv) ^ GF_S_GF_16_2);
    tmp = base_gf->inverse.w16(base_gf, tmp);

    d = base_gf->multiply.w16(base_gf, d, tmp);

    c0 = base_gf->multiply.w16(base_gf, (d^1), a0inv);
    c1 = base_gf->multiply.w16(base_gf, d, a1inv);
  }

  c = c0 | (c1 << 16);

  return c;
}

static
gf_val_32_t
gf_w32_composite_divide(gf_t *gf, gf_val_32_t a, gf_val_32_t b)
{
  gf_val_32_t binv;

  binv = gf_w32_composite_inverse(gf, b);

  return gf_w32_composite_multiply(gf, a, binv);
}

static
void
gf_w32_composite_multiply_region_table(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  int i=0;
  struct gf_w16_logtable_data * ltd;
  uint16_t b0 = val & 0x0000ffff;
  uint16_t b1 = (val & 0xffff0000) >> 16;
  uint32_t *s32 = (uint32_t *) src;
  uint32_t *d32 = (uint32_t *) dest;
  uint16_t a0, a1, a1b1;
  int num_syms = bytes >> 2;
  int sym_divisible = bytes % 4;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w32_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w32_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  ltd = (struct gf_w16_logtable_data *) h->private;

  if (xor) {
    for (i = 0;i < num_syms; i++) {
      a0 = s32[i] & 0x0000ffff;
      a1 = (s32[i] & 0xffff0000) >> 16;
      a1b1 = ltd->antilog_tbl[ltd->log_tbl[a1] + ltd->log_tbl[b1]];

      d32[i] ^= ((ltd->antilog_tbl[ltd->log_tbl[a0] + ltd->log_tbl[b0]] ^ a1b1) | 
                 ((ltd->antilog_tbl[ltd->log_tbl[a1] + ltd->log_tbl[b0]] ^ ltd->antilog_tbl[ltd->log_tbl[a0] + ltd->log_tbl[b1]] ^ 
                   ltd->antilog_tbl[ltd->log_tbl[a1b1] + ltd->log_tbl[GF_S_GF_16_2]]) << 16));

    }
  } else {
    for (i = 0;i < num_syms; i++) {
      a0 = s32[i] & 0x0000ffff;
      a1 = (s32[i] & 0xffff0000) >> 16;
      a1b1 = ltd->antilog_tbl[ltd->log_tbl[a1] + ltd->log_tbl[b1]];

      d32[i] = ((ltd->antilog_tbl[ltd->log_tbl[a0] + ltd->log_tbl[b0]] ^ a1b1) | 
                 ((ltd->antilog_tbl[ltd->log_tbl[a1] + ltd->log_tbl[b0]] ^ ltd->antilog_tbl[ltd->log_tbl[a0] + ltd->log_tbl[b1]] ^ 
                   ltd->antilog_tbl[ltd->log_tbl[a1b1] + ltd->log_tbl[GF_S_GF_16_2]]) << 16));
    }
  }
}

static
void
gf_w32_composite_multiply_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  unsigned long uls, uld;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  int i=0;
  struct gf_w16_logtable_data * ltd;
  uint16_t b0 = val & 0x0000ffff;
  uint16_t b1 = (val & 0xffff0000) >> 16;
  uint32_t *s32 = (uint32_t *) src;
  uint32_t *d32 = (uint32_t *) dest;
  uint16_t a0, a1, a1b1;
  int num_syms = bytes >> 2;
  int sym_divisible = bytes % 4;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;
  if ((uls & 0x7) != (uld & 0x7)) gf_alignment_error("gf_w32_buf_const_log", 2);
  if (sym_divisible) {
    gf_alignment_error("gf_w32_buf_const_log: buffer size not divisible by symbol size = 2 bytes", 2);
  }

  if (val == 0) {
    if (xor) return;
    bzero(dest, bytes);
    return;
  }

  ltd = (struct gf_w16_logtable_data *) h->private;

  if (xor) {
    for (i = 0;i < num_syms; i++) {
      a0 = s32[i] & 0x0000ffff;
      a1 = (s32[i] & 0xffff0000) >> 16;
      a1b1 = base_gf->multiply.w16(base_gf, a1, b1);

      d32[i] ^= ((base_gf->multiply.w16(base_gf, a0, b0) ^ a1b1) |
                ((base_gf->multiply.w16(base_gf, a1, b0) ^ base_gf->multiply.w16(base_gf, a0, b1) ^ base_gf->multiply.w16(base_gf, a1b1, GF_S_GF_16_2)) << 16)); 

    }
  } else {
    for (i = 0;i < num_syms; i++) {
      a0 = s32[i] & 0x0000ffff;
      a1 = (s32[i] & 0xffff0000) >> 16;
      a1b1 = base_gf->multiply.w16(base_gf, a1, b1);

      d32[i] = ((base_gf->multiply.w16(base_gf, a0, b0) ^ a1b1) |
                ((base_gf->multiply.w16(base_gf, a1, b0) ^ base_gf->multiply.w16(base_gf, a0, b1) ^ base_gf->multiply.w16(base_gf, a1b1, GF_S_GF_16_2)) << 16)); 
    }
  }
}



static
void
gf_w32_composite_multiply_region_alt(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_val_16_t val0 = val & 0x0000ffff;
  gf_val_16_t val1 = (val & 0xffff0000) >> 16;
  int sub_reg_size = bytes / 2;

  if (bytes % 2 != 0) gf_alignment_error("gf_w32_composite_multiply_region_alt", 1);
  if (sub_reg_size % 2 != 0) gf_alignment_error("gf_w32_composite_multiply_region_alt", 1);
  
  if (!xor) {
    memset(dest, 0, bytes);
  }

  base_gf->multiply_region.w16(base_gf, src, dest, val0, sub_reg_size, xor);
  base_gf->multiply_region.w16(base_gf, src+sub_reg_size, dest, val1, sub_reg_size, 1);
  base_gf->multiply_region.w16(base_gf, src, dest+sub_reg_size, val1, sub_reg_size, xor);
  base_gf->multiply_region.w16(base_gf, src+sub_reg_size, dest+sub_reg_size, val0, sub_reg_size, 1);
  base_gf->multiply_region.w16(base_gf, src+sub_reg_size, dest+sub_reg_size, base_gf->multiply.w16(base_gf, GF_S_GF_16_2, val1), sub_reg_size, 1);
}

static
int gf_w32_composite_init(gf_t *gf)
{
  struct gf_w16_logtable_data *ltd;
  gf_internal_t *h = (gf_internal_t *) gf->scratch;
  gf_t *base_gf = h->base_gf;
  gf_val_32_t a, b;
  uint64_t prim_poly = ((gf_internal_t *) base_gf->scratch)->prim_poly;
  int i;

  ltd = (struct gf_w16_logtable_data *) h->private;

  ltd->log_tbl[0] = 0;

  bzero(&(ltd->_antilog_tbl[0]), sizeof(ltd->_antilog_tbl));

  ltd->antilog_tbl = &(ltd->_antilog_tbl[GF_BASE_FIELD_SIZE * 2]);

  b = 1;
  for (i = 0; i < GF_BASE_FIELD_GROUP_SIZE; i++) {
      ltd->log_tbl[b] = (gf_val_16_t)i;
      ltd->antilog_tbl[i] = (gf_val_16_t)b;
      ltd->antilog_tbl[i+GF_BASE_FIELD_GROUP_SIZE] = (gf_val_16_t)b;
      b <<= 1;
      if (b & GF_BASE_FIELD_SIZE) {
          b = b ^ prim_poly;
      }
  }
  ltd->inv_tbl[0] = 0;  /* Not really, but we need to fill it with something  */
  ltd->inv_tbl[1] = 1;
  for (i = 2; i < GF_BASE_FIELD_SIZE; i++) {
    ltd->inv_tbl[i] = ltd->antilog_tbl[GF_BASE_FIELD_GROUP_SIZE-ltd->log_tbl[i]];
  }

  if (h->region_type & GF_REGION_ALTMAP) {
    gf->multiply_region.w32 = gf_w32_composite_multiply_region_alt;
  } else {
    if (h->region_type & GF_REGION_SINGLE_TABLE) {
      gf->multiply_region.w32 = gf_w32_composite_multiply_region_table;
    } else {
      gf->multiply_region.w32 = gf_w32_composite_multiply_region;
    }
  }

  gf->multiply.w32 = gf_w32_composite_multiply;
  gf->divide.w32 = gf_w32_composite_divide;
  gf->inverse.w32 = gf_w32_composite_inverse;

  return 1;
}

int gf_w32_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2)
{
  int ss;

  ss = (GF_REGION_SSE | GF_REGION_NOSSE);
  switch(mult_type)
  {
    case GF_MULT_SPLIT_TABLE: 
        if (arg1 == 8 && arg2 == 8){
          if (region_type != GF_REGION_DEFAULT) return -1;
          return sizeof(gf_internal_t) + sizeof(struct gf_split_8_8_data) + 64;
        }
        if ((arg1 == 2 && arg2 == 32) || (arg2 == 2 && arg1 == 32)) {
          region_type &= (~GF_REGION_LAZY);
          if ((region_type & ss) == ss) return -1;
          if ((region_type | ss) != ss) return -1;
          return sizeof(gf_internal_t) + sizeof(struct gf_split_2_32_lazy_data) + 64;
        }
        if ((arg1 == 4 && arg2 == 32) || (arg2 == 4 && arg1 == 32)) {
          region_type &= (~GF_REGION_LAZY);
          if (region_type & GF_REGION_ALTMAP) {
            region_type &= (~GF_REGION_ALTMAP);
            if ((region_type & ss) == ss) return -1;
            if ((region_type | ss) != ss) return -1;
            return sizeof(gf_internal_t) + sizeof(struct gf_split_4_32_lazy_data) + 64;
          } else return -1;
        }
        return -1;
    case GF_MULT_DEFAULT:
    case GF_MULT_SHIFT:
      if (arg1 != 0 || arg2 != 0 || region_type != 0) return -1;
      return sizeof(gf_internal_t);
      break;
    case GF_MULT_COMPOSITE:
      if (region_type & ~(GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP | GF_REGION_STDMAP)) return -1;
      if ((region_type & (GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP)) == (GF_REGION_SINGLE_TABLE | GF_REGION_ALTMAP)) return -1;
      if (arg1 == 2 && arg2 == 16 || arg2 == 2 && arg1 == 16) {
        return sizeof(gf_internal_t) + sizeof(struct gf_w16_logtable_data) + 64;
      } else {
        return -1;
      }
    default:
      return -1;
   }
}

int gf_w32_init(gf_t *gf)
{
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  if (h->prim_poly == 0) h->prim_poly = 0x400007;

  gf->multiply.w32 = NULL;
  gf->divide.w32 = NULL;
  gf->inverse.w32 = NULL;
  gf->multiply_region.w32 = NULL;

  switch(h->mult_type) {
    case GF_MULT_DEFAULT: 
    case GF_MULT_SHIFT:     if (gf_w32_shift_init(gf) == 0) return 0; break;
    case GF_MULT_COMPOSITE: if (gf_w32_composite_init(gf) == 0) return 0; break;
    case GF_MULT_SPLIT_TABLE: if (gf_w32_split_init(gf) == 0) return 0; break;
    default: return 0;
  }
  if (h->divide_type == GF_DIVIDE_EUCLID) {
    gf->divide.w32 = gf_w32_divide_from_inverse;
    gf->inverse.w32 = gf_w32_euclid;
  } else if (h->divide_type == GF_DIVIDE_MATRIX) {
    gf->divide.w32 = gf_w32_divide_from_inverse;
    gf->inverse.w32 = gf_w32_matrix;
  }

  if (gf->inverse.w32 != NULL && gf->divide.w32 == NULL) {
    gf->divide.w32 = gf_w32_divide_from_inverse;
  }
  if (gf->inverse.w32 == NULL && gf->divide.w32 != NULL) {
    gf->inverse.w32 = gf_w32_inverse_from_divide;
  }
  return 1;
}

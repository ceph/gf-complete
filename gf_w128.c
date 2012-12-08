/*
 * gf_w128.c
 *
 * Routines for 128-bit Galois fields
 */

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

#define GF_FIELD_WIDTH (128)

#define two_x(a) {\
  a[0] <<= 1; \
  if (a[1] & (uint64_t) 1 << 63) a[0] ^= 1; \
  a[1] <<= 1; }
  
#define a_get_b(a, i, b, j) {\
  a[i] = b[j]; \
  a[i + 1] = b[j + 1];}

#define set_zero(a, i) {\
  a[i] = 0; \
  a[i + 1] = 0;}

typedef struct gf_group_tables_s {
  gf_val_128_t m_table;
  gf_val_128_t r_table;
} gf_group_tables_t;

static
void
gf_w128_multiply_region_from_single(gf_t *gf, void *src, void *dest, gf_val_128_t val, int bytes,
int xor)
{
    int i;
    gf_val_128_t s128;
    gf_val_128_t d128;
    uint64_t c128[2];

    set_zero(c128, 0);

    s128 = (gf_val_128_t) src;
    d128 = (gf_val_128_t) dest;

    if (xor) {
      for (i = 0; i < bytes/sizeof(gf_val_64_t); i += 2) {
        gf->multiply.w128(gf, &s128[i], val, c128);
        d128[i] ^= c128[0];
        d128[i+1] ^= c128[1];
      }
    } else {
      for (i = 0; i < bytes/sizeof(gf_val_64_t); i += 2) {
        gf->multiply.w128(gf, &s128[i], val, &d128[i]);
      }
    }
}

/*
 * Some w128 notes:
 * --Big Endian
 * --return values allocated beforehand
 */
void
gf_w128_shift_multiply(gf_t *gf, gf_val_128_t a128, gf_val_128_t b128, gf_val_128_t c128)
{
  /* ordered highest bit to lowest l[0] l[1] r[0] r[1] */
  uint64_t pl[2], pr[2], ppl[2], ppr[2], i, a[2], bl[2], br[2], one, lbit;
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;

  if (GF_W128_IS_ZERO(a128) || GF_W128_IS_ZERO(b128)) {
    set_zero(c128, 0);
    return;
  }

  a_get_b(a, 0, a128, 0);
  a_get_b(br, 0, b128, 0);
  set_zero(bl, 0);

  one = 1;
  lbit = (one << 63);

  set_zero(pl, 0);
  set_zero(pr, 0);

  for (i = 0; i < GF_FIELD_WIDTH/2; i++) {
    if (a[1] & (one << i)) {
      pl[1] ^= bl[1];
      pr[0] ^= br[0];
      pr[1] ^= br[1];
    }
    bl[1] <<= 1;
    if (br[0] & lbit) bl[1] ^= 1;
    br[0] <<= 1;
    if (br[1] & lbit) br[0] ^= 1;
    br[1] <<= 1;
  }

  for (i = 0; i < GF_FIELD_WIDTH/2; i++) {
    if (a[0] & (one << i)) {
      pl[0] ^= bl[0];
      pl[1] ^= bl[1];
      pr[0] ^= br[0];
    }
    bl[0] <<= 1;
    if (bl[1] & lbit) bl[0] ^= 1;
    bl[1] <<= 1;
    if (br[0] & lbit) bl[1] ^= 1;
    br[0] <<= 1;
  }

  one = lbit;
  ppl[0] = lbit;
  ppl[1] = h->prim_poly >> 1;
  ppr[0] = lbit;
  ppr[1] = 0;
  while (one != 0) {
    if (pl[0] & one) {
      pl[0] ^= ppl[0];
      pl[1] ^= ppl[1];
      pr[0] ^= ppr[0];
      pr[1] ^= ppr[1];
    }
    one >>= 1;
    ppr[1] >>= 1;
    if (ppr[0] & 1) ppr[1] ^= lbit;
    ppr[0] >>= 1;
    if (ppl[1] & 1) ppr[0] ^= lbit;
    ppl[1] >>= 1;
    if (ppl[0] & 1) ppl[1] ^= lbit;
    ppl[0] >>= 1;
  }

  one = lbit;
  while (one != 0) {
    if (pl[1] & one) {
      pl[1] ^= ppl[1];
      pr[0] ^= ppr[0];
      pr[1] ^= ppr[1];
    }
    one >>= 1;
    ppr[1] >>= 1;
    if (ppr[0] & 1) ppr[1] ^= lbit;
    ppr[0] >>= 1;
    if (ppl[1] & 1) ppr[0] ^= lbit;
    ppl[1] >>= 1;
  }

  c128[0] = pr[0];
  c128[1] = pr[1];

  return;
}

static
void gf_w128_group_m_init(gf_t *gf, gf_val_128_t b128)
{
  int i, j;
  int g_m;
  uint64_t prim_poly, lbit;
  gf_internal_t *scratch;
  gf_group_tables_t *gt;
  uint64_t a128[2];
  scratch = (gf_internal_t *) gf->scratch;
  gt = scratch->private;
  g_m = scratch->arg1;
  prim_poly = scratch->prim_poly;

  set_zero(gt->m_table, 0);
  a_get_b(gt->m_table, 2, b128, 0);
  lbit = 1;
  lbit <<= 63;

  for (i = 2; i < (1 << g_m); i <<= 1) {
    a_get_b(a128, 0, gt->m_table, 2 * (i >> 1));
    two_x(a128);
    a_get_b(gt->m_table, 2 * i, a128, 0);
    if (gt->m_table[2 * (i >> 1)] & lbit) gt->m_table[(2 * i) + 1] ^= prim_poly;
    for (j = 0; j < i; j++) {
      gt->m_table[(2 * i) + (2 * j)] = gt->m_table[(2 * i)] ^ gt->m_table[(2 * j)];
      gt->m_table[(2 * i) + (2 * j) + 1] = gt->m_table[(2 * i) + 1] ^ gt->m_table[(2 * j) + 1];
    }
  }
  return;
}

void
gf_w128_group_multiply(GFP gf, gf_val_128_t a128, gf_val_128_t b128, gf_val_128_t c128)
{
  int i;
  /* index_r, index_m, total_m (if g_r > g_m) */
  int i_r, i_m, t_m;
  int mask_m, mask_r;
  int g_m, g_r;
  uint64_t p_i[2], a[2];
  gf_internal_t *scratch;
  gf_group_tables_t *gt;

  scratch = (gf_internal_t *) gf->scratch;
  gt = scratch->private;
  g_m = scratch->arg1;
  g_r = scratch->arg2;

  mask_m = (1 << g_m) - 1;
  mask_r = (1 << g_r) - 1;

  if (b128[0] != gt->m_table[2] || b128[1] != gt->m_table[3]) {
    gf_w128_group_m_init(gf, b128);
  }

  p_i[0] = 0;
  p_i[1] = 0;
  a[0] = a128[0];
  a[1] = a128[1];

  t_m = 0;
  i_r = 0;

  /* Top 64 bits */
  for (i = ((GF_FIELD_WIDTH / 2) / g_m) - 1; i >= 0; i--) {
    i_m = (a[0] >> (i * g_m)) & mask_m;
    i_r ^= (p_i[0] >> (64 - g_m)) & mask_r;
    p_i[0] <<= g_m;
    p_i[0] ^= (p_i[1] >> (64-g_m));
    p_i[1] <<= g_m;
    p_i[0] ^= gt->m_table[2 * i_m];
    p_i[1] ^= gt->m_table[(2 * i_m) + 1];
    t_m += g_m;
    if (t_m == g_r) {
      p_i[1] ^= gt->r_table[i_r];
      t_m = 0;
      i_r = 0;
    } else {
      i_r <<= g_m;
    }
  }

  for (i = ((GF_FIELD_WIDTH / 2) / g_m) - 1; i >= 0; i--) {
    i_m = (a[1] >> (i * g_m)) & mask_m;
    i_r ^= (p_i[0] >> (64 - g_m)) & mask_r;
    p_i[0] <<= g_m;
    p_i[0] ^= (p_i[1] >> (64-g_m));
    p_i[1] <<= g_m;
    p_i[0] ^= gt->m_table[2 * i_m];
    p_i[1] ^= gt->m_table[(2 * i_m) + 1];
    t_m += g_m;
    if (t_m == g_r) {
      p_i[1] ^= gt->r_table[i_r];
      t_m = 0;
      i_r = 0;
    } else {
      i_r <<= g_m;
    }
  }

  c128[0] = p_i[0];
  c128[1] = p_i[1];
}

/* a^-1 -> b */
void
gf_w128_euclid(GFP gf, gf_val_128_t a128, gf_val_128_t b128)
{
  uint64_t e_i[2], e_im1[2], e_ip1[2];
  uint64_t d_i, d_im1, d_ip1;
  uint64_t y_i[2], y_im1[2], y_ip1[2];
  uint64_t c_i[2];
  uint64_t *b;
  uint64_t one = 1;
  uint64_t buf, buf1;

  /* This needs to return some sort of error (in b128?) */
  if (a128[0] == 0 && a128[1] == 0) return;

  e_im1[0] = 0;
  e_im1[1] = ((gf_internal_t *) (gf->scratch))->prim_poly;
  e_i[0] = a128[0];
  e_i[1] = a128[1];
  d_im1 = 128;
  for (d_i = (d_im1-1) % 64; ((one << d_i) & e_i[0]) == 0 && d_i > 0; d_i--) ;
  if (!((one << d_i) & e_i[0])) {
    for (d_i = (d_im1-1) % 64; ((one << d_i) & e_i[1] == 0); d_i--) ;
  } else {
    d_i += 64;
  }
  y_i[0] = 0;
  y_i[1] = 1;
  y_im1[0] = 0;
  y_im1[1] = 0;

  while (!(e_i[0] == 0 && e_i[1] == 1)) {

    e_ip1[0] = e_im1[0];
    e_ip1[1] = e_im1[1];
    d_ip1 = d_im1;
    c_i[0] = 0;
    c_i[1] = 0;

    while (d_ip1 >= d_i) {
      if ((d_ip1 - d_i) >= 64) {
        c_i[0] ^= (one << ((d_ip1 - d_i) - 64));
        e_ip1[0] ^= (e_i[1] << ((d_ip1 - d_i) - 64));
      } else {
        c_i[1] ^= (one << (d_ip1 - d_i));
        e_ip1[0] ^= (e_i[0] << (d_ip1 - d_i));
        if (d_ip1 - d_i > 0) e_ip1[0] ^= (e_i[1] >> (64 - (d_ip1 - d_i)));
        e_ip1[1] ^= (e_i[1] << (d_ip1 - d_i));
      }
        d_ip1--;
      while (d_ip1 >= 64 && (e_ip1[0] & (one << (d_ip1 - 64))) == 0) d_ip1--;
      while (d_ip1 <  64 && (e_ip1[1] & (one << d_ip1)) == 0) d_ip1--;
    }

    gf->multiply.w128(gf, c_i, y_i, y_ip1);
    y_ip1[0] ^= y_im1[0];
    y_ip1[1] ^= y_im1[1];

    y_im1[0] = y_i[0];
    y_im1[1] = y_i[1];

    y_i[0] = y_ip1[0];
    y_i[1] = y_ip1[1];

    e_im1[0] = e_i[0];
    e_im1[1] = e_i[1];
    d_im1 = d_i;
    e_i[0] = e_ip1[0];
    e_i[1] = e_ip1[1];
    d_i = d_ip1;
  }

  b = (uint64_t *) b128;
  b[0] = y_i[0];
  b[1] = y_i[1];

  return;
}

void
gf_w128_divide_from_inverse(GFP gf, gf_val_128_t a128, gf_val_128_t b128, gf_val_128_t c128)
{
  uint64_t d[2];
  gf->inverse.w128(gf, b128, d);
  gf->multiply.w128(gf, a128, d, c128);
  return;
}

void
gf_w128_inverse_from_divide(GFP gf, gf_val_128_t a128, gf_val_128_t b128)
{
  uint64_t one128[2];
  one128[0] = 0;
  one128[1] = 1;
  gf->divide.w128(gf, one128, a128, b128);
  return;
}

static
int gf_w128_shift_init(gf_t *gf)
{
  gf->multiply.w128 = gf_w128_shift_multiply;
  gf->inverse.w128 = gf_w128_euclid;
  gf->multiply_region.w128 = gf_w128_multiply_region_from_single;
  return 1;
}

/*
 * Because the prim poly is only 8 bits and we are limiting g_r to 16, I do not need the high 64
 * bits in all of these numbers.
 */
static
void gf_w128_group_r_init(gf_t *gf)
{
  int i, j;
  int g_r;
  uint64_t pp;
  gf_internal_t *scratch;
  gf_group_tables_t *gt;
  scratch = (gf_internal_t *) gf->scratch;
  gt = scratch->private;
  g_r = scratch->arg2;
  pp = scratch->prim_poly;

  gt->r_table[0] = 0;
  for (i = 1; i < (1 << g_r); i++) {
    gt->r_table[i] = 0;
    for (j = 0; j < g_r; j++) {
      if (i & (1 << j)) {
        gt->r_table[i] ^= (pp << j);
      }
    }
  }
  return;
}

static
int gf_w128_group_init(gf_t *gf)
{
  gf_internal_t *scratch;
  gf_group_tables_t *gt;
  int g_m, g_r, size_r;

  scratch = (gf_internal_t *) gf->scratch;
  gt = scratch->private;
  g_m = scratch->arg1;
  g_r = scratch->arg2;
  size_r = (1 << g_r);

  gt->r_table = scratch->private + (2 * sizeof(uint64_t *));
  gt->m_table = gt->r_table + size_r;
  gt->m_table[2] = 0;
  gt->m_table[3] = 0;

  gf_w128_group_r_init(gf);

  gf->multiply.w128 = gf_w128_group_multiply;
  gf->inverse.w128 = gf_w128_euclid;
  gf->multiply_region.w128 = gf_w128_multiply_region_from_single; /* This needs to change */
  return 1;
}

int gf_w128_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2)
{
  int size_m, size_r;
  int w = 128;
  switch(mult_type)
  {
    case GF_MULT_DEFAULT:
    case GF_MULT_SHIFT:
      if (arg1 != 0 || arg2 != 0 || region_type != 0) return -1;
      return sizeof(gf_internal_t);
      break;
    case GF_MULT_GROUP:

      /* arg1 == mult size, arg2 == reduce size */
      /* Should prevent anything over arg1 > 16 || arg2 > 16 */
      if (region_type != 0) return -1;
      if (arg1 <= 0 || arg2 <= 0 || arg1 > 16 || arg2 > 16) return -1;
      if (GF_FIELD_WIDTH % arg1 != 0 || GF_FIELD_WIDTH % arg2 != 0) return -1;
      /*
       * Currently implementing code where g_m and g_r are the same or where g_r is larger, as
       * these it is more efficient to have g_r as large as possible (but still not > 16)
       */
      if (arg1 > arg2) return -1;

      /* size of each group, 128 bits */
      size_m = (1 << arg1) * 2 * sizeof(uint64_t);
      /* The PP is only 8 bits and we are limiting g_r to 16, so only uint64_t */
      size_r = (1 << arg2) * sizeof(uint64_t);

      /* 
       * two pointers prepend the table data for structure
       * because the tables are of dynamic size
       */
      return sizeof(gf_internal_t) + size_m + size_r + 2 * sizeof(uint64_t *);
    default:
      return -1;
   }
}

int gf_w128_init(gf_t *gf)
{
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  if (h->prim_poly == 0) h->prim_poly = 0x87; /* Omitting the leftmost 1 as in w=32 */

  gf->multiply.w128 = NULL;
  gf->divide.w128 = NULL;
  gf->inverse.w128 = NULL;
  gf->multiply_region.w128 = NULL;

  switch(h->mult_type) {
    case GF_MULT_DEFAULT: 
    case GF_MULT_SHIFT:        if (gf_w128_shift_init(gf) == 0) return 0; break;
    case GF_MULT_GROUP:        if (gf_w128_group_init(gf) == 0) return 0; break;
    default: return 0;
  }
  if (h->divide_type == GF_DIVIDE_EUCLID) {
    gf->divide.w128 = gf_w128_divide_from_inverse;
    gf->inverse.w128 = gf_w128_euclid;
  } /* } else if (h->divide_type == GF_DIVIDE_MATRIX) {
    gf->divide.w128 = gf_w128_divide_from_inverse;
    gf->inverse.w128 = gf_w128_matrix;
  } */

  if (gf->inverse.w128 != NULL && gf->divide.w128 == NULL) {
    gf->divide.w128 = gf_w128_divide_from_inverse;
  }
  if (gf->inverse.w128 == NULL && gf->divide.w128 != NULL) {
    gf->inverse.w128 = gf_w128_inverse_from_divide;
  }
  return 1;
}

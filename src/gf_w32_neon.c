/*
 * GF-Complete: A Comprehensive Open Source Library for Galois Field Arithmetic
 * James S. Plank, Ethan L. Miller, Kevin M. Greenan,
 * Benjamin A. Arnold, John A. Burnum, Adam W. Disney, Allen C. McBride.
 *
 * gf_w32_neon.c
 *
 * Neon routines for 32-bit Galois fields
 *
 * Author: Janne Grunau <j@jannau.net>
 *
 */


#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>
#include "gf_w32.h"

#ifndef ARCH_AARCH64
#define vqtbl1q_u8(tbl, v) vcombine_u8(vtbl2_u8(tbl, vget_low_u8(v)),   \
                                       vtbl2_u8(tbl, vget_high_u8(v)))
#endif

static
void
neon_w32_split_4_32_multiply_region(gf_t *gf, uint32_t *src, uint32_t *dst,
                                    uint32_t *d_end, uint8_t btable[8][4][16],
                                    uint32_t val, int xor, int altmap)
{
  int i, j;
#ifdef ARCH_AARCH64
  uint8x16_t tables[8][4];
#else
  uint8x8x2_t tables[8][4];
#endif
  uint32x4_t v0, v1, v2, v3, s0, s1, s2, s3;
  uint8x16_t p0, p1, p2, p3, si, mask1;
  uint16x8x2_t r0, r1;
  uint8x16x2_t q0, q1;

  for (i = 0; i < 8; i++) {
    for (j = 0; j < 4; j++) {
#ifdef ARCH_AARCH64
      tables[i][j] = vld1q_u8(btable[i][j]);
#else
      tables[i][j].val[0] = vld1_u8(btable[i][j]);
      tables[i][j].val[1] = vld1_u8(btable[i][j] + 8);
#endif
    }
  }

  mask1 = vdupq_n_u8(0xf);

  while (dst < d_end) {

      v0 = vld1q_u32(src); src += 4;
      v1 = vld1q_u32(src); src += 4;
      v2 = vld1q_u32(src); src += 4;
      v3 = vld1q_u32(src); src += 4;

      if (altmap) {
          q0.val[0] = vreinterpretq_u8_u32(v0);
          q0.val[1] = vreinterpretq_u8_u32(v1);
          q1.val[0] = vreinterpretq_u8_u32(v2);
          q1.val[1] = vreinterpretq_u8_u32(v3);
      } else {
          r0 = vtrnq_u16(vreinterpretq_u16_u32(v0), vreinterpretq_u16_u32(v2));
          r1 = vtrnq_u16(vreinterpretq_u16_u32(v1), vreinterpretq_u16_u32(v3));

          q0 = vtrnq_u8(vreinterpretq_u8_u16(r0.val[0]),
                        vreinterpretq_u8_u16(r1.val[0]));
          q1 = vtrnq_u8(vreinterpretq_u8_u16(r0.val[1]),
                        vreinterpretq_u8_u16(r1.val[1]));
      }

      si = vandq_u8(q0.val[0], mask1);
      p0 = vqtbl1q_u8(tables[0][0], si);
      p1 = vqtbl1q_u8(tables[0][1], si);
      p2 = vqtbl1q_u8(tables[0][2], si);
      p3 = vqtbl1q_u8(tables[0][3], si);

      si = vshrq_n_u8(q0.val[0], 4);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[1][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[1][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[1][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[1][3], si));

      si = vandq_u8(q0.val[1], mask1);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[2][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[2][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[2][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[2][3], si));

      si = vshrq_n_u8(q0.val[1], 4);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[3][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[3][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[3][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[3][3], si));

      si = vandq_u8(q1.val[0], mask1);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[4][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[4][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[4][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[4][3], si));

      si = vshrq_n_u8(q1.val[0], 4);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[5][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[5][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[5][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[5][3], si));

      si = vandq_u8(q1.val[1], mask1);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[6][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[6][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[6][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[6][3], si));

      si = vshrq_n_u8(q1.val[1], 4);
      p0 = veorq_u8(p0, vqtbl1q_u8(tables[7][0], si));
      p1 = veorq_u8(p1, vqtbl1q_u8(tables[7][1], si));
      p2 = veorq_u8(p2, vqtbl1q_u8(tables[7][2], si));
      p3 = veorq_u8(p3, vqtbl1q_u8(tables[7][3], si));

      if (altmap) {
          s0 = vreinterpretq_u32_u8(p0);
          s1 = vreinterpretq_u32_u8(p1);
          s2 = vreinterpretq_u32_u8(p2);
          s3 = vreinterpretq_u32_u8(p3);
      } else {
          q0 = vtrnq_u8(p0, p1);
          q1 = vtrnq_u8(p2, p3);

          r0 = vtrnq_u16(vreinterpretq_u16_u8(q0.val[0]),
                         vreinterpretq_u16_u8(q1.val[0]));
          r1 = vtrnq_u16(vreinterpretq_u16_u8(q0.val[1]),
                         vreinterpretq_u16_u8(q1.val[1]));

          s0 = vreinterpretq_u32_u16(r0.val[0]);
          s1 = vreinterpretq_u32_u16(r1.val[0]);
          s2 = vreinterpretq_u32_u16(r0.val[1]);
          s3 = vreinterpretq_u32_u16(r1.val[1]);
      }

      if (xor) {
          v0 = vld1q_u32(dst);
          v1 = vld1q_u32(dst + 4);
          v2 = vld1q_u32(dst + 8);
          v3 = vld1q_u32(dst + 12);
          s0 = veorq_u32(s0, v0);
          s1 = veorq_u32(s1, v1);
          s2 = veorq_u32(s2, v2);
          s3 = veorq_u32(s3, v3);
      }

      vst1q_u32(dst,      s0);
      vst1q_u32(dst + 4,  s1);
      vst1q_u32(dst + 8,  s2);
      vst1q_u32(dst + 12, s3);

      dst += 16;
  }
}

static
inline
void
neon_w32_split_4_32_lazy_multiply_region(gf_t *gf, void *src, void *dest, uint32_t val, int bytes, int xor, int altmap)
{
  gf_internal_t *h;
  int i, j, k;
  uint32_t pp, v, *s32, *d32, *top, tmp_table[16];
  uint8_t btable[8][4][16];
  gf_region_data rd;

  if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
  if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

  h = (gf_internal_t *) gf->scratch;
  pp = h->prim_poly;

  gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 64);
  gf_do_initial_region_alignment(&rd);

  s32 = (uint32_t *) rd.s_start;
  d32 = (uint32_t *) rd.d_start;
  top = (uint32_t *) rd.d_top;

  v = val;
  for (i = 0; i < 8; i++) {
    tmp_table[0] = 0;
    for (j = 1; j < 16; j <<= 1) {
      for (k = 0; k < j; k++) {
        tmp_table[k^j] = (v ^ tmp_table[k]);
      }
      v = (v & GF_FIRST_BIT) ? ((v << 1) ^ pp) : (v << 1);
    }
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 16; k++) {
        btable[i][j][k] = (uint8_t) tmp_table[k];
        tmp_table[k] >>= 8;
      }
    }
  }

  if (xor)
    neon_w32_split_4_32_multiply_region(gf, s32, d32, top, btable, val, 1, altmap);
  else
    neon_w32_split_4_32_multiply_region(gf, s32, d32, top, btable, val, 0, altmap);

  gf_do_final_region_alignment(&rd);
}

static
void
gf_w32_split_4_32_lazy_multiply_region_neon(gf_t *gf, void *src, void *dest,
                                            gf_val_32_t val, int bytes, int xor)
{
  neon_w32_split_4_32_lazy_multiply_region(gf, src, dest, val, bytes, xor, 0);
}

static
void
gf_w32_split_4_32_lazy_altmap_multiply_region_neon(gf_t *gf, void *src,
                                                   void *dest, gf_val_32_t val,
                                                   int bytes, int xor)
{
  neon_w32_split_4_32_lazy_multiply_region(gf, src, dest, val, bytes, xor, 1);
}

void gf_w32_neon_split_init(gf_t *gf)
{
  gf_internal_t *h = (gf_internal_t *) gf->scratch;

  if (h->region_type & GF_REGION_ALTMAP)
      gf->multiply_region.w32 = gf_w32_split_4_32_lazy_altmap_multiply_region_neon;
  else
      gf->multiply_region.w32 = gf_w32_split_4_32_lazy_multiply_region_neon;

}

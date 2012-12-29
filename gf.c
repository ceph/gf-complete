/*
 * gf.c
 *
 * Generic routines for Galois fields
 */

#include "gf_int.h"
#include <stdio.h>
#include <stdlib.h>

int gf_scratch_size(int w, 
                    int mult_type, 
                    int region_type, 
                    int divide_type, 
                    int arg1, 
                    int arg2)
{
  switch(w) {
    case 4: return gf_w4_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    case 8: return gf_w8_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    case 16: return gf_w16_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    case 32: return gf_w32_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    case 64: return gf_w64_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    case 128: return gf_w128_scratch_size(mult_type, region_type, divide_type, arg1, arg2);
    default: return gf_wgen_scratch_size(w, mult_type, region_type, divide_type, arg1, arg2);
  }
}

int gf_dummy_init(gf_t *gf)
{
  return 0;
}

int gf_init_easy(gf_t *gf, int w, int mult_type)
{
  return gf_init_hard(gf, w, mult_type, GF_REGION_DEFAULT, GF_DIVIDE_DEFAULT, 0, 0, 0, NULL, NULL);
}

int gf_init_hard(gf_t *gf, int w, int mult_type, 
                        int region_type,
                        int divide_type,
                        uint64_t prim_poly,
                        int arg1, int arg2,
                        gf_t *base_gf,
                        void *scratch_memory) 
{
  int sz;
  gf_internal_t *h;
  
  sz = gf_scratch_size(w, mult_type, region_type, divide_type, arg1, arg2);

  if (sz <= 0) return 0;

  if (scratch_memory == NULL) {
    h = (gf_internal_t *) malloc(sz);
    h->free_me = 1;
  } else {
    h = scratch_memory;
    h->free_me = 0;
  }
  gf->scratch = (void *) h;
  h->mult_type = mult_type;
  h->region_type = region_type;
  h->divide_type = divide_type;
  h->w = w;
  h->prim_poly = prim_poly;
  h->arg1 = arg1;
  h->arg2 = arg2;
  h->base_gf = base_gf;
  h->private = (void *) gf->scratch;
  h->private += (sizeof(gf_internal_t));
  gf->extract_word.w32 = NULL;

  //printf("Created w=%d, with mult_type=%d and region_type=%d\n", w, mult_type, region_type);

  switch(w) {
    case 4: return gf_w4_init(gf);
    case 8: return gf_w8_init(gf);
    case 16: return gf_w16_init(gf);
    case 32: return gf_w32_init(gf);
    case 64: return gf_w64_init(gf);
    case 128: return gf_w128_init(gf);
    default: return gf_wgen_init(gf);
  }
}

int gf_free(gf_t *gf, int recursive)
{
  gf_internal_t *h;

  h = (gf_internal_t *) gf->scratch;
  if (recursive && h->base_gf != NULL) {
    gf_free(h->base_gf, 1);
    free(h->base_gf);
  }
  if (h->free_me) free(h);
}

void gf_alignment_error(char *s, int a)
{
  fprintf(stderr, "Alignment error in %s:\n", s);
  fprintf(stderr, "   The source and destination buffers must be aligned to each other,\n");
  fprintf(stderr, "   and they must be aligned to a %d-byte address.\n", a);
  exit(1);
}

/* Lifted this code from Jens Gregor -- thanks, Jens */

int gf_is_sse2()
{
  unsigned int cpeinfo;
  unsigned int cpsse;
  asm ( "mov $0x1, %%eax\n\t"
                "cpuid\n\t"
                "mov %%edx, %0\n\t"
            "mov %%ecx, %1\n" : "=m" (cpeinfo), "=m" (cpsse));
  if ((cpeinfo >> 26) & 0x1 ) return 1;
  return 0;
}

static 
void gf_invert_binary_matrix(int *mat, int *inv, int rows) {
  int cols, i, j, k;
  int tmp;

  cols = rows;

  for (i = 0; i < rows; i++) inv[i] = (1 << i);

  /* First -- convert into upper triangular */

  for (i = 0; i < cols; i++) {

    /* Swap rows if we ave a zero i,i element.  If we can't swap, then the
       matrix was not invertible */

    if ((mat[i] & (1 << i)) == 0) {
      for (j = i+1; j < rows && (mat[j] & (1 << i)) == 0; j++) ;
      if (j == rows) {
        fprintf(stderr, "galois_invert_matrix: Matrix not invertible!!\n");
        exit(1);
      }
      tmp = mat[i]; mat[i] = mat[j]; mat[j] = tmp;
      tmp = inv[i]; inv[i] = inv[j]; inv[j] = tmp;
    }

    /* Now for each j>i, add A_ji*Ai to Aj */
    for (j = i+1; j != rows; j++) {
      if ((mat[j] & (1 << i)) != 0) {
        mat[j] ^= mat[i];
        inv[j] ^= inv[i];
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down */

  for (i = rows-1; i >= 0; i--) {
    for (j = 0; j < i; j++) {
      if (mat[j] & (1 << i)) {
        /*  mat[j] ^= mat[i]; */
        inv[j] ^= inv[i];
      }
    }
  }
}

uint32_t gf_bitmatrix_inverse(uint32_t y, int w, uint32_t pp) 
{
  uint32_t mat[32], inv[32], mask;
  int i;

  mask = (w == 32) ? 0xffffffff : (1 << w) - 1;
  for (i = 0; i < w; i++) {
    mat[i] = y;

    if (y & (1 << (w-1))) {
      y = y << 1;
      y = ((y ^ pp) & mask);
    } else {
      y = y << 1;
    }
  }

  gf_invert_binary_matrix(mat, inv, w);
  return inv[0];
}

/*
void gf_two_byte_region_table_multiply(gf_region_data *rd, uint16_t *base)
{
  uint64_t p, ta, shift, tb;
  uint64_t *s64, *d64

  s64 = rd->s_start;
  d64 = rd->d_start;
  
  while (s64 < (uint64_t *) rd->s_top) {
    p = (rd->xor) ? *d64 : 0;
    ta = *s64;

    shift = 0;
    while (ta != 0) {
      tb = base[ta&0xffff];
      p ^= (tb << shift);
      ta >>= 16;
      shift += 16;
    }

    *d64 = p;
    d64++;
    s64++;
  }
}
*/

void gf_two_byte_region_table_multiply(gf_region_data *rd, uint16_t *base)
{
  uint64_t a, prod;
  int j, xor;
  uint64_t *s64, *d64, *top;

  s64 = rd->s_start;
  d64 = rd->d_start;
  top = rd->d_top;
  xor = rd->xor;
  
  if (xor) {
    while (d64 != top) {
      a = *s64;
      prod = base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      prod ^= *d64;
      *d64 = prod;
      *s64++;
      *d64++;
    }
  } else {
    while (d64 != top) {
      a = *s64;
      prod = base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      a <<= 16;
      prod <<= 16;
      prod ^= base[a >> 48];
      *d64 = prod;
      *s64++;
      *d64++;
    }
  }
}

static void gf_slow_multiply_region(gf_region_data *rd, void *src, void *dest, void *s_top)
{
  uint8_t *s8, *d8;
  uint16_t *s16, *d16;
  uint32_t *s32, *d32;
  uint64_t *s64, *d64;
  gf_internal_t *h;
  int wb;
  uint32_t p, a;

  h = rd->gf->scratch;
  wb = (h->w)/8;
  if (wb == 0) wb = 1;
  
  while (src < s_top) {
    switch (h->w) {
    case 8:
      s8 = (uint8_t *) src;
      d8 = (uint8_t *) dest;
      *d8 = (rd->xor) ? (*d8 ^ rd->gf->multiply.w32(rd->gf, rd->val, *s8)) : 
                      rd->gf->multiply.w32(rd->gf, rd->val, *s8);
      break;
    case 4:
      s8 = (uint8_t *) src;
      d8 = (uint8_t *) dest;
      a = *s8;
      p = rd->gf->multiply.w32(rd->gf, rd->val, a&0xf);
      p |= (rd->gf->multiply.w32(rd->gf, rd->val, a >> 4) << 4);
      if (rd->xor) p ^= *d8;
      *d8 = p;
      break;
    case 16:
      s16 = (uint16_t *) src;
      d16 = (uint16_t *) dest;
      *d16 = (rd->xor) ? (*d16 ^ rd->gf->multiply.w32(rd->gf, rd->val, *s16)) : 
                      rd->gf->multiply.w32(rd->gf, rd->val, *s16);
      break;
    case 32:
      s32 = (uint32_t *) src;
      d32 = (uint32_t *) dest;
      *d32 = (rd->xor) ? (*d32 ^ rd->gf->multiply.w32(rd->gf, rd->val, *s32)) : 
                      rd->gf->multiply.w32(rd->gf, rd->val, *s32);
      break;
    case 64:
      s64 = (uint64_t *) src;
      d64 = (uint64_t *) dest;
      *d64 = (rd->xor) ? (*d64 ^ rd->gf->multiply.w64(rd->gf, rd->val, *s64)) : 
                      rd->gf->multiply.w64(rd->gf, rd->val, *s64);
      break;
    default:
      fprintf(stderr, "Error: gf_slow_multiply_region: w=%d not implemented.\n", h->w);
      exit(1);
    }
    src += wb;
    dest += wb;
  }
}

/* If align>16, you align to 16 bytes, but make sure that within the aligned region bytes is a multiple of align.  However, you make sure that the region itself is a multiple of align. 

   If align = -1, then this is cauchy.  You need to make sure that bytes is a multiple of w. */

void gf_set_region_data(gf_region_data *rd,
  gf_t *gf,
  void *src,
  void *dest,
  int bytes,
  uint32_t val,
  int xor,
  int align)
{
  uint8_t *s8, *d8;
  gf_internal_t *h;
  int wb;
  uint32_t a;
  unsigned long uls, uld;

  h = gf->scratch;
  wb = (h->w)/8;
  if (wb == 0) wb = 1;
  
  rd->gf = gf;
  rd->src = src;
  rd->dest = dest;
  rd->bytes = bytes;
  rd->val = val;
  rd->xor = xor;
  rd->align = align;

  uls = (unsigned long) src;
  uld = (unsigned long) dest;

  a = (align <= 16) ? align : 16;

  if (align == -1) { /* This is cauchy.  Error check bytes, then set up the pointers
                        so that there is no alignment regions. */
    if (bytes % h->w != 0) {
      fprintf(stderr, "Error in region multiply operation.\n");
      fprintf(stderr, "The size must be a multiple of %d bytes.\n", h->w);
      exit(1);
    }
  
    rd->s_start = src;
    rd->d_start = dest;
    rd->s_top = src + bytes;
    rd->d_top = src + bytes;
    return;
  }

  if (uls % a != uld % a) {
    fprintf(stderr, "Error in region multiply operation.\n");
    fprintf(stderr, "The source & destination pointers must be aligned with respect\n");
    fprintf(stderr, "to each other along a %d byte boundary.\n", a);
    fprintf(stderr, "Src = 0x%lx.  Dest = 0x%lx\n", (unsigned long) src,
            (unsigned long) dest);
    exit(1);
  }

  if (uls % wb != 0) {
    fprintf(stderr, "Error in region multiply operation.\n");
    fprintf(stderr, "The pointers must be aligned along a %d byte boundary.\n", wb);
    fprintf(stderr, "Src = 0x%lx.  Dest = 0x%lx\n", (unsigned long) src,
            (unsigned long) dest);
    exit(1);
  }

  if (bytes % wb != 0) {
    fprintf(stderr, "Error in region multiply operation.\n");
    fprintf(stderr, "The size must be a multiple of %d bytes.\n", wb);
    exit(1);
  }

  uls %= a;
  if (uls != 0) uls = (align-uls);
  rd->s_start = rd->src + uls;
  rd->d_start = rd->dest + uls;
  bytes -= uls;

  bytes -= (bytes % align);
  rd->s_top = rd->s_start + bytes;
  rd->d_top = rd->d_start + bytes;
}

void gf_do_initial_region_alignment(gf_region_data *rd)
{
  gf_slow_multiply_region(rd, rd->src, rd->dest, rd->s_start);
}

void gf_do_final_region_alignment(gf_region_data *rd)
{
  gf_slow_multiply_region(rd, rd->s_top, rd->d_top, rd->src+rd->bytes);
}

void gf_multby_zero(void *dest, int bytes, int xor) 
{
  if (xor) return;
  bzero(dest, bytes);
  return;
}

void gf_multby_one(gf_t *gf, void *src, void *dest, int bytes, int xor) 
{
#ifdef   INTEL_SSE4
  __m128i ms, md;
#endif
  uint8_t *s8, *d8, *dtop8;
  uint64_t *s64, *d64, *dtop64;
  int abytes;

  gf_region_data rd;
  if (!xor) {
    memcpy(dest, src, bytes);
    return;
  }

#ifdef   INTEL_SSE4
  s8 = (uint8_t *) src;
  d8 = (uint8_t *) dest;
  abytes = bytes & 0xfffffff0;

  while (d8 < (uint8_t *) dest + abytes) {
    ms = _mm_loadu_si128 ((__m128i *)(s8));
    md = _mm_loadu_si128 ((__m128i *)(d8));
    md = _mm_xor_si128(md, ms);
    _mm_storeu_si128((__m128i *)(d8), md);
    s8 += 16;
    d8 += 16;
  }
  while (d8 != (uint8_t *) dest+bytes) {
    *d8 ^= *s8;
    d8++;
    s8++;
  }
  return;
#endif

  /* If you don't have SSE, you'd better be aligned..... */

  gf_set_region_data(&rd, gf, src, dest, bytes, 1, xor, 8);
  s8 = (uint8_t *) src;
  d8 = (uint8_t *) dest;
  while (d8 != rd.d_start) {
    *d8 ^= *s8;
    d8++;
    s8++;
  }
  dtop64 = (uint64_t *) rd.d_top;

  while (d64 < dtop64) {
    *d64 ^= *s64;
    d64++;
    s64++;
  }
  while (d8 != (uint8_t *) dest+bytes) {
    *d8 ^= *s8;
    d8++;
    s8++;
  }
  return;
}

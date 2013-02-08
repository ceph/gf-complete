/*
 * gf_int.h
 *
 * Internal code for Galois field routines.
 */

#pragma once

#include "gf_complete.h"

#include <string.h>

extern void     timer_start (double *t);
extern double   timer_split (const double *t);
extern void     galois_fill_random (void *buf, int len, unsigned int seed);

extern int      gf_is_sse();

typedef struct {
  int mult_type;
  int region_type;
  int divide_type;
  int w;
  uint64_t prim_poly;
  int free_me;
  int arg1;
  int arg2;
  gf_t *base_gf;
  void *private;
} gf_internal_t;

extern int gf_w4_init (gf_t *gf);
extern int gf_w4_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_w8_init (gf_t *gf);
extern int gf_w8_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_w16_init (gf_t *gf);
extern int gf_w16_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_w32_init (gf_t *gf);
extern int gf_w32_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_w64_init (gf_t *gf);
extern int gf_w64_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_w128_init (gf_t *gf);
extern int gf_w128_scratch_size(int mult_type, int region_type, int divide_type, int arg1, int arg2);

extern int gf_wgen_init (gf_t *gf);
extern int gf_wgen_scratch_size(int w, int mult_type, int region_type, int divide_type, int arg1, int arg2);

void gf_wgen_cauchy_region(gf_t *gf, void *src, void *dest, gf_val_32_t val, int bytes, int xor);
gf_val_32_t gf_wgen_extract_word(gf_t *gf, void *start, int bytes, int index);


extern void gf_alignment_error(char *s, int a);

extern uint32_t gf_bitmatrix_inverse(uint32_t y, int w, uint32_t pp);

/* This structure lets you define a region multiply.  It helps because you can handle
   unaligned portions of the data with the procedures below, which really cleans
   up the code. */

typedef struct {
  gf_t *gf;
  void *src;
  void *dest;
  int bytes;
  uint64_t val;
  int xor;
  int align;           /* The number of bytes to which to align. */
  void *s_start;       /* The start and the top of the aligned region. */
  void *d_start;
  void *s_top;
  void *d_top;
} gf_region_data;

/* This lets you set up one of these in one call. It also sets the start/top pointers. */

void gf_set_region_data(gf_region_data *rd,
                        gf_t *gf,
                        void *src,
                        void *dest,
                        int bytes,
                        uint64_t val,
                        int xor,
                        int align);

/* This performs gf->multiply.32() on all of the unaligned bytes in the beginning of the region */

extern void gf_do_initial_region_alignment(gf_region_data *rd);

/* This performs gf->multiply.32() on all of the unaligned bytes in the end of the region */

extern void gf_do_final_region_alignment(gf_region_data *rd);

extern void gf_two_byte_region_table_multiply(gf_region_data *rd, uint16_t *base);

extern void gf_multby_zero(void *dest, int bytes, int xor);

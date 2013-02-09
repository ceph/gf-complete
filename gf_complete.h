/* gf_complete.h
 * External include file for Galois field arithmetic.  */

#pragma once
#include <stdint.h>

#ifdef  INTEL_SSE4
#include <nmmintrin.h>
#include <emmintrin.h>
#endif

/* This does either memcpy or xor, depending on "xor" */

extern void gf_multby_one(void *src, void *dest, int bytes, int xor);

#define GF_W128_IS_ZERO(val) (val[0] == 0 && val[1] == 0)
#define GF_W128_EQUAL(val1, val2) ((val1[0] == val2[0]) && (val1[1] == val2[1]))

/* These are the different ways to perform multiplication.
   Not all are implemented for all values of w.
   See the paper for an explanation of how they work. */

typedef enum {GF_MULT_DEFAULT,   
              GF_MULT_SHIFT,   
              GF_MULT_GROUP,   
              GF_MULT_BYTWO_p,
              GF_MULT_BYTWO_b,
              GF_MULT_TABLE,   
              GF_MULT_LOG_TABLE,   
              GF_MULT_SPLIT_TABLE,   
              GF_MULT_COMPOSITE } gf_mult_type_t;

/* These are the different ways to optimize region 
   operations.  They are bits because you can compose them:
   You can mix SINGLE/DOUBLE/QUAD, LAZY, SSE/NOSSE, STDMAP/ALTMAP/CAUCHY.
   Certain optimizations only apply to certain gf_mult_type_t's.  
   Again, please see documentation for how to use these */
   
#define GF_REGION_DEFAULT      (0x0)
#define GF_REGION_SINGLE_TABLE (0x1)
#define GF_REGION_DOUBLE_TABLE (0x2)
#define GF_REGION_QUAD_TABLE   (0x4)
#define GF_REGION_LAZY         (0x8)
#define GF_REGION_SSE          (0x10)
#define GF_REGION_NOSSE        (0x20)
#define GF_REGION_STDMAP       (0x40)
#define GF_REGION_ALTMAP       (0x80)
#define GF_REGION_CAUCHY       (0x100)

typedef uint32_t gf_region_type_t;

/* These are different ways to implement division.
   Once again, it's best to use "DEFAULT".  However,
   there are times when you may want to experiment
   with the others. */

typedef enum { GF_DIVIDE_DEFAULT,
               GF_DIVIDE_MATRIX,
               GF_DIVIDE_EUCLID } gf_division_type_t;

/* We support w=4,8,16,32,64 and 128 with their own data types and
   operations for multiplication, division, etc.  We also support
   a "gen" type so that you can do general gf arithmetic for any 
   value of w from 1 to 32.  You can perform a "region" operation
   on these if you use "CAUCHY" as the mapping. 
 */

typedef uint32_t    gf_val_32_t;
typedef uint64_t    gf_val_64_t;
typedef uint64_t   *gf_val_128_t;

typedef struct gf *GFP;

typedef union gf_func_a_b {
    gf_val_32_t  (*w32) (GFP gf, gf_val_32_t a,  gf_val_32_t b);
    gf_val_64_t  (*w64) (GFP gf, gf_val_64_t a,  gf_val_64_t b);
    void         (*w128)(GFP gf, gf_val_128_t a, gf_val_128_t b, gf_val_128_t c);
} gf_func_a_b;
  
typedef union {
  gf_val_32_t  (*w32) (GFP gf, gf_val_32_t a);
  gf_val_64_t  (*w64) (GFP gf, gf_val_64_t a);
  void         (*w128)(GFP gf, gf_val_128_t a, gf_val_128_t b);
} gf_func_a;
  
typedef union {
  void  (*w32) (GFP gf, void *src, void *dest, gf_val_32_t val,  int bytes, int add);
  void  (*w64) (GFP gf, void *src, void *dest, gf_val_64_t val,  int bytes, int add);
  void  (*w128)(GFP gf, void *src, void *dest, gf_val_128_t val, int bytes, int add);
} gf_region;

typedef union {
  gf_val_32_t  (*w32) (GFP gf, void *start, int bytes, int index);
  gf_val_64_t  (*w64) (GFP gf, void *start, int bytes, int index);
  void         (*w128)(GFP gf, void *start, int bytes, int index, gf_val_128_t rv);
} gf_extract;

typedef struct gf {
  gf_func_a_b    multiply;
  gf_func_a_b    divide;
  gf_func_a      inverse;
  gf_region      multiply_region;
  gf_extract     extract_word;
  void           *scratch;
} gf_t;
    
extern int gf_init_easy(GFP gf, int w);

extern int gf_init_hard(GFP gf, 
                        int w, 
                        int mult_type, 
                        int region_type, 
                        int divide_type, 
                        uint64_t prim_poly,
                        int arg1, 
                        int arg2,
                        GFP base_gf,
                        void *scratch_memory);

extern int gf_scratch_size(int w, 
                           int mult_type, 
                           int region_type, 
                           int divide_type, 
                           int arg1, 
                           int arg2);

extern int gf_free(GFP gf, int recursive);

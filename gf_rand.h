/* gf_rand.h
 * External include file for random number generation.  */

#pragma once
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* These are all pretty self-explanatory */
uint32_t MOA_Random_32();
uint64_t MOA_Random_64();
void     MOA_Random_128(uint64_t *x);
uint32_t MOA_Random_W(int w, int zero_ok);
void MOA_Fill_Random_Region (void *reg, int size);   /* reg should be aligned to 4 bytes, but
                                                        size can be anything. */
void     MOA_Seed(uint32_t seed);



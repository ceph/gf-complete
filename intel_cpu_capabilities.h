/*
 * Routines to figure out what an Intel CPU's capabilities are.
 */

#pragma once

#include <stdint.h>

/* Words in CPE_INFO */
#define CPU_CPE_INFO            0x1000
#define CPU_CAP_MMX             (CPU_CPE_INFO | 23)
#define CPU_CAP_SSE             (CPU_CPE_INFO | 25)
#define CPU_CAP_SSE2            (CPU_CPE_INFO | 26)

/* Words in CPSSE */
#define CPU_CPSSE               0x2000
#define CPU_CAP_SSE3            (CPU_CPSSE | 0)
#define CPU_CAP_PCLMULQDQ       (CPU_CPSSE | 1)
#define CPU_CAP_SSSE3           (CPU_CPSSE | 10)
#define CPU_CAP_SSE41           (CPU_CPSSE | 19)
#define CPU_CAP_SSE42           (CPU_CPSSE | 20)
#define CPU_CAP_AVX             (CPU_CPSSE | 28)

#define cpuid(func,ax,bx,cx,dx)\
        __asm__ __volatile__ ("cpuid":\
                              "=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func));

inline
int
cpu_has_feature (unsigned which)
{
        uint32_t        cpeinfo;
        uint32_t        cpsse;
        uint32_t        a, b;

        cpuid(1, a, b, cpsse, cpeinfo);
        if (which & CPU_CPE_INFO) {
                return (!! ((cpeinfo >> (which & 0xff)) & 0x1) );
        } else if (which & CPU_CPSSE) {
                return (!! ((cpsse >> (which & 0xff)) & 0x1) );
        } else {
                return (0);
        }
}

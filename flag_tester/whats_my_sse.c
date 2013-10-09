/*
 * whats_my_sse.c - lifted from Jens Gregor -- thanks Jens
 */

#include <stdio.h>
#include <stdlib.h>
#include "intel_cpu_capabilities.h"

struct {
        unsigned        feature_code;
        const char *    feature_name;
} features[] = {
        {CPU_CAP_MMX, "MMX"},
        {CPU_CAP_SSE, "SSE"},
        {CPU_CAP_SSE2, "SSE2"},
        {CPU_CAP_SSE3, "SSE3"},
        {CPU_CAP_SSSE3, "SSSE3"},
        {CPU_CAP_SSE41, "SSE4.1"},
        {CPU_CAP_SSE42, "SSE4.2"},
        {CPU_CAP_PCLMULQDQ, "PCLMULQDQ"},
        {CPU_CAP_AVX, "AVX"},
        {0, NULL},
};


int main()
{
        unsigned        i;

        printf ("CPU has the following instruction set features enabled: " );
        for (i = 0; features[i].feature_name != NULL; ++i) {
                if (cpu_has_feature (features[i].feature_code)) {
                        printf ("%s ", features[i].feature_name);
                }
        }
        printf ("\nCPU is missing the following instruction set features:  ");
        for (i = 0; features[i].feature_name != NULL; ++i) {
                if (! cpu_has_feature (features[i].feature_code)) {
                        printf ("%s ", features[i].feature_name);
                }
        }
        printf ("\n");
}

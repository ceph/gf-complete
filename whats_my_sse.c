#include <stdio.h>
#include <stdlib.h>
int main()
{
        unsigned int cpeinfo;
        unsigned int cpsse;
    asm (
        "mov $0x1, %%eax\n\t"
                "cpuid\n\t"
                "mov %%edx, %0\n\t"
            "mov %%ecx, %1\n" : "=m" (cpeinfo), "=m" (cpsse)
                );

        printf("1 - Instruction set is supported by CPU\n");
        printf("0 - Instruction set not supported\n"); 
        printf("--------------------------------\n");
        printf("MMX:     %d\n", ((cpeinfo >> 23) & 0x1 ));
        printf("SSE:     %d\n", ((cpeinfo >> 25) & 0x1 ));
        printf("SSE2:    %d\n", ((cpeinfo >> 26) & 0x1 ));
        printf("SSE3:    %d\n", ((cpsse       ) & 0x1 ));
        printf("SSE4.1:  %d\n", ((cpsse >> 19) & 0x1 ));
        printf("SSE4.2:  %d\n", ((cpsse >> 20) & 0x1 ));
        printf("SSEAVX:  %d\n", ((cpsse >> 28) & 0x1 ));
}

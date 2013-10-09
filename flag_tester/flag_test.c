/*
 * flag_test.c - copied from whats_my_sse.c to output proper compile
 *  flags for the GNUmakefile
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "intel_cpu_capabilities.h"

void usage()
{
  fprintf(stderr, "usage: flag_test <compiler name>\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  //make sure to extend these buffers if more flags are added to this program
  char cflags[1000], ldflags[1000], buf[1000];
  FILE *file;
  char sse_found = 0;

  if(argc != 2)
    usage();

  sprintf(cflags, "CFLAGS = -O3");
  sprintf(ldflags, "LDFLAGS = -O3");

  if(cpu_has_feature(CPU_CAP_SSE42))
  {
    sprintf(buf, "%s sse_test.c -o sse4 -msse4 -DSSE4 2> /dev/null", argv[1]);
    system(buf);
    if(file = fopen("sse4", "r"))
    {
      fclose(file);

      //run program and compare to the included output
      system("./sse4 > temp.txt 2> /dev/null");
      system("diff sse4_test.txt temp.txt > diff.txt 2> /dev/null");
      file = fopen("diff.txt", "r");
      if(fgetc(file) == EOF)
      {
        strcat(cflags, " -msse4 -DINTEL_SSE4");
        strcat(ldflags, " -msse4");
        sse_found = 1;
      }
      fclose(file);
    }
  }

  if(cpu_has_feature(CPU_CAP_SSSE3) && !sse_found)
  {
    sprintf(buf, "%s sse_test.c -o ssse3 -mssse3 -DSSSE3 2> /dev/null", argv[1]);
    system(buf);
    if(file = fopen("ssse3", "r"))
    {
      fclose(file);

      //run program and compare to the included output
      system("./ssse3 > temp.txt 2> /dev/null");
      system("diff ssse3_test.txt temp.txt > diff.txt 2> /dev/null");
      file = fopen("diff.txt", "r");
      if(fgetc(file) == EOF)
      {
        strcat(cflags, " -mssse3 -DINTEL_SSSE3");
        strcat(ldflags, " -mssse3");
        sse_found = 1;
      }
      fclose(file);
    }
  }

  if(cpu_has_feature(CPU_CAP_SSE2) && !sse_found)
  {
    sprintf(buf, "%s sse_test.c -o sse2 -msse2 -DSSE2 2> /dev/null", argv[1]);
    system(buf);
    if(file = fopen("sse2", "r"))
    {
      fclose(file);

      //run program and compare to the included output
      system("./sse2 > temp.txt 2> /dev/null");
      system("diff sse2_test.txt temp.txt > diff.txt 2> /dev/null");
      file = fopen("diff.txt", "r");
      if(fgetc(file) == EOF)
      {
        strcat(cflags, " -msse2 -DINTEL_SSE2");
        strcat(ldflags, " -msse2");
        sse_found = 1;
      }
      fclose(file);
    }
  }

  if(cpu_has_feature(CPU_CAP_PCLMULQDQ) && sse_found)
  {
    sprintf(buf, "%s pclmul_test.c -o pclmul -maes -mpclmul 2> /dev/null"
      , argv[1]);
    system(buf);
    if(file = fopen("pclmul", "r"))
    {
      fclose(file);

      //run program and compare to the included output
      system("./pclmul > temp.txt 2> /dev/null");
      system("diff pclmul_test.txt temp.txt > diff.txt 2> /dev/null");
      file = fopen("diff.txt", "r");
      if(fgetc(file) == EOF)
      {
        strcat(cflags, " -maes -mpclmul -DINTEL_PCLMUL");
        strcat(ldflags, " -maes -mpclmul");
      }
      fclose(file);
    }
  }

  printf("%s\n%s\n", cflags, ldflags);
}

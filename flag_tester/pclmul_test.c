#include <wmmintrin.h>
#include <stdint.h>
#include <stdio.h>

#define MM_PRINT8(s, r) { uint8_t blah[16], ii; printf("%-20s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 1) printf("%s%02x", (ii%4==0) ? "   " : " ", blah[15-ii]); printf("\n"); }


int main()
{
  uint64_t answer;
  uint32_t pp;
  __m128i a, b, c;

  a = _mm_set1_epi8(0x0D);
  b = _mm_set_epi32(0,0,0,0x0A);
  pp = 0x13;
  MM_PRINT8("a", a);
  MM_PRINT8("b", b);

  c = _mm_clmulepi64_si128(a, b, 0);
  MM_PRINT8("a clm b", c);

  a = _mm_set1_epi8(0xf0);
  MM_PRINT8("a", a);
  b = _mm_and_si128(a, c);
  b = _mm_srli_epi64(b, 4);
  MM_PRINT8("shifted", b);


  a = _mm_set_epi32(0,0,0,pp);
  MM_PRINT8("PP", a);

  b = _mm_clmulepi64_si128(a, b, 0);
  MM_PRINT8("PP clm over", b);

  c = _mm_xor_si128(c,b);
  MM_PRINT8("Answer", c);
  //answer = _mm_extract_epi64(c, 0);
  //printf("%llx\n", answer);
}

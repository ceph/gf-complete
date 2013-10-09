#ifdef SSE4
#define SSSE3
#include <nmmintrin.h>
#endif

#ifdef SSSE3
#define SSE2
#include <tmmintrin.h>
#endif

#ifdef SSE2
#include <emmintrin.h>
#endif

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#define MM_PRINT8(s, r) { uint8_t blah[16], ii; printf("%-20s", s); _mm_storeu_si128((__m128i *)blah, r); for (ii = 0; ii < 16; ii += 1) printf("%s%02x", (ii%4==0) ? "   " : " ", blah[15-ii]); printf("\n"); }

int main()
{
  uint32_t u32;
  uint64_t u64;
  uint8_t *ui8 = malloc(20), i;
  __m128i a, b, c, d;

  for(i=0; i < 20; i++)
    ui8[i] = i;

  a = _mm_load_si128( (__m128i *) ui8 );
  b = _mm_loadu_si128( (__m128i *) (ui8+1));
  c = _mm_loadu_si128( (__m128i *) (ui8+2));
  d = _mm_loadu_si128( (__m128i *) (ui8+3));

  MM_PRINT8("a", a);
  MM_PRINT8("b", b);
  MM_PRINT8("c", c);
  MM_PRINT8("d", d);

  a = _mm_slli_epi16(a, 2);
  b = _mm_slli_epi32(b, 2);
  c = _mm_slli_epi64(c, 2);
  d = _mm_slli_si128(d, 2);

  MM_PRINT8("a sl16", a);
  MM_PRINT8("b sl32", b);
  MM_PRINT8("c sl64", c);
  MM_PRINT8("d sl128", d);

  a = _mm_srli_epi16(a, 2);
  b = _mm_srli_epi32(b, 2);
  c = _mm_srli_epi64(c, 2);
  d = _mm_srli_si128(d, 2);

  MM_PRINT8("a sr16", a);
  MM_PRINT8("b sr32", b);
  MM_PRINT8("c sr64", c);
  MM_PRINT8("d sr128", d);

  d = _mm_xor_si128(a, b);
  MM_PRINT8("d = a^b", d);
  
  d = _mm_sub_epi8(a, b);
  MM_PRINT8("d = a-b epi8", d);
  
  d = _mm_sub_epi16(a, b);
  MM_PRINT8("d = a-b epi16", d);
  
  d = _mm_sub_epi32(a, b);
  MM_PRINT8("d = a-b epi32", d);
  
  d = _mm_sub_epi64(a, b);
  MM_PRINT8("d = a-b epi64", d);
  
  d = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
  MM_PRINT8("d set_epi8", d);
  
  d = _mm_set_epi32(0x12345678, 0x9abcdef0, 0x12345678, 0x9abcdef0);
  MM_PRINT8("d set_epi32", d);
  
  d = _mm_set1_epi64x(0xF0F0F0F0F0F0F0F0ULL);
  MM_PRINT8("d set1_epi64", d);
  
  d = _mm_set1_epi32(0xe2e2e2e2);
  MM_PRINT8("d set1_epi32", d);

  d = _mm_set1_epi16(0xaff3);
  MM_PRINT8("d set1_epi16", d);

  d = _mm_set1_epi8(0xc5);
  MM_PRINT8("d set1_epi8", d);

  d = _mm_packus_epi16(d, d);
  MM_PRINT8("d packus_epi16(d,d)", d);

  c = _mm_unpackhi_epi8(a, d);
  MM_PRINT8("c unpackhi(a,d)", c);

  b = _mm_unpacklo_epi8(c, a);
  MM_PRINT8("b unpacklo(c,a)", b);

  d = _mm_and_si128(d, b);
  MM_PRINT8("d and(d,b)", d);

  _mm_store_si128( (__m128i *) ui8, a);
  printf("a stored to mem: ");
  for(i=0; i < 16; i++)
    printf("%u ", ui8[i]);
  printf("\n");

  d = _mm_setzero_si128();
  MM_PRINT8("d setzero", d);

  u32 = 0xABCD1234;
  u64 = 0xFEDCBA1291827364ULL;
  
  #ifdef SSE4
  d = _mm_insert_epi32(d, u32, 2);
  MM_PRINT8("d insert32 @ 2", d);

  u32 = 0;
  u32 = _mm_extract_epi32(d, 2);
  printf("extract_epi32 @ 2: %x\n", u32);

  d = _mm_insert_epi64(d, u64, 0);
  MM_PRINT8("d insert64 @ 0", d);

  u64 = 0;
  u64 = _mm_extract_epi64(d, 0);
  printf("extract_epi64 @ 0: %" PRIx64 "\n", u64);
  #endif

  c = _mm_set1_epi8(5);
  MM_PRINT8("c", c);

  #ifdef SSSE3
  a = _mm_shuffle_epi8(b, c);
  MM_PRINT8("a shuffle(b, c)", a);
  #endif

}

#include "oblas_lite.h"

u8 GF2_8_MUL[256 * 256];

void obl_fill_mul_tab() {
  for (int i = 0; i < 256; i++) {
    for (int j = 0; j < 256; j++) {
      if (i == 0 || j == 0)
        GF2_8_MUL[256 * i + j] = 0;
      else
        GF2_8_MUL[256 * i + j] = GF2_8_EXP[GF2_8_LOG[i] + GF2_8_LOG[j]];
    }
  }
}

static void obl_axpy_ref(u8 *a, u8 *b, u8 u, int k) {
  register u8 *u_row = &GF2_8_MUL[u << 8];
  register u8 *ap = a, *ae = &a[k], *bp = b;
  for (; ap != ae; ap++, bp++)
    *ap ^= u_row[*bp];
}

static void obl_scal_ref(u8 *a, u8 *b, u8 u, int k) {
  register u8 *u_row = &GF2_8_MUL[u << 8];
  register u8 *ap = a, *ae = &a[k];
  for (; ap != ae; ap++)
    *ap = u_row[*ap];
}

#if defined(__AVX512F__)
#include <immintrin.h>

#define OBL_SHUF(op, a, b, f)                                                  \
  do {                                                                         \
    const u8 *u_lo = GF2_8_SHUF_LO + u * 16;                                   \
    const u8 *u_hi = GF2_8_SHUF_HI + u * 16;                                   \
    const __m512i mask = _mm512_set1_epi8(0x0f);                               \
    const __m128i ulo_128 = _mm_loadu_si128((__m128i *)u_lo);                  \
    const __m128i uhi_128 = _mm_loadu_si128((__m128i *)u_hi);                  \
    const __m512i urow_lo = _mm512_broadcast_i32x4(ulo_128);                   \
    const __m512i urow_hi = _mm512_broadcast_i32x4(uhi_128);                   \
    __m512i *ap = (__m512i *)a,                                                \
            *ae = (__m512i *)(a + k - (k % sizeof(__m256i))),                  \
            *bp = (__m512i *)b;                                                \
    for (; ap < ae; ap++, bp++) {                                              \
      __m512i bx = _mm512_loadu_si512(bp);                                     \
      __m512i lo = _mm512_and_si512(bx, mask);                                 \
      bx = _mm512_srli_epi64(bx, 4);                                           \
      __m512i hi = _mm512_and_si512(bx, mask);                                 \
      lo = _mm512_shuffle_epi8(urow_lo, lo);                                   \
      hi = _mm512_shuffle_epi8(urow_hi, hi);                                   \
      _mm512_storeu_si512(                                                     \
          ap, f(_mm512_loadu_si512(ap), _mm512_xor_si512(lo, hi)));            \
    }                                                                          \
    op##_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(__m512i));                      \
  } while (0)
#define OBL_SHUF_XOR _mm512_xor_si512

#else
#if defined(__AVX2__)
#include <immintrin.h>

#define OBL_SHUF(op, a, b, f)                                                  \
  do {                                                                         \
    const u8 *u_lo = GF2_8_SHUF_LO + u * 16;                                   \
    const u8 *u_hi = GF2_8_SHUF_HI + u * 16;                                   \
    const __m256i mask = _mm256_set1_epi8(0x0f);                               \
    const __m256i urow_lo =                                                    \
        _mm256_loadu2_m128i((__m128i *)u_lo, (__m128i *)u_lo);                 \
    const __m256i urow_hi =                                                    \
        _mm256_loadu2_m128i((__m128i *)u_hi, (__m128i *)u_hi);                 \
    __m256i *ap = (__m256i *)a,                                                \
            *ae = (__m256i *)(a + k - (k % sizeof(__m256i))),                  \
            *bp = (__m256i *)b;                                                \
    for (; ap < ae; ap++, bp++) {                                              \
      __m256i bx = _mm256_loadu_si256(bp);                                     \
      __m256i lo = _mm256_and_si256(bx, mask);                                 \
      bx = _mm256_srli_epi64(bx, 4);                                           \
      __m256i hi = _mm256_and_si256(bx, mask);                                 \
      lo = _mm256_shuffle_epi8(urow_lo, lo);                                   \
      hi = _mm256_shuffle_epi8(urow_hi, hi);                                   \
      _mm256_storeu_si256(                                                     \
          ap, f(_mm256_loadu_si256(ap), _mm256_xor_si256(lo, hi)));            \
    }                                                                          \
    op##_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(__m256i));                      \
  } while (0)
#define OBL_SHUF_XOR _mm256_xor_si256

#else
#if defined(__SSSE3__) || (defined(_MSC_VER) && defined(_M_X64))
#include <emmintrin.h>
#include <tmmintrin.h>

#define OBL_SHUF(op, a, b, f)                                                  \
  do {                                                                         \
    const u8 *u_lo = GF2_8_SHUF_LO + u * 16;                                   \
    const u8 *u_hi = GF2_8_SHUF_HI + u * 16;                                   \
    const __m128i mask = _mm_set1_epi8(0x0f);                                  \
    const __m128i urow_lo = _mm_loadu_si128((__m128i *)u_lo);                  \
    const __m128i urow_hi = _mm_loadu_si128((__m128i *)u_hi);                  \
    __m128i *ap = (__m128i *)a,                                                \
            *ae = (__m128i *)(a + k - (k % sizeof(__m128i))),                  \
            *bp = (__m128i *)b;                                                \
    for (; ap < ae; ap++, bp++) {                                              \
      __m128i bx = _mm_loadu_si128(bp);                                        \
      __m128i lo = _mm_and_si128(bx, mask);                                    \
      bx = _mm_srli_epi64(bx, 4);                                              \
      __m128i hi = _mm_and_si128(bx, mask);                                    \
      lo = _mm_shuffle_epi8(urow_lo, lo);                                      \
      hi = _mm_shuffle_epi8(urow_hi, hi);                                      \
      _mm_storeu_si128(ap, f(_mm_loadu_si128(ap), _mm_xor_si128(lo, hi)));     \
    }                                                                          \
    op##_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(__m128i));                      \
  } while (0)
#define OBL_SHUF_XOR _mm_xor_si128

#else
#if defined(__GNUC__) && defined(__aarch64__)

#include <arm_neon.h>
/* FIXME: support armv7 with mfpu=neon and synthesized vqtbl1q_u8 */
#define OBL_SHUF(op, a, b, f)                                                  \
  do {                                                                         \
    const u8 *u_lo = GF2_8_SHUF_LO + u * 16;                                   \
    const u8 *u_hi = GF2_8_SHUF_HI + u * 16;                                   \
    u8 *ap = a, *ae = a + k - (k % sizeof(uint8x16_t)), *bp = (u8 *)b;         \
    uint8x16_t mask = vdupq_n_u8(0x0f);                                        \
    uint8x16_t urow_lo = vld1q_u8(u_lo);                                       \
    uint8x16_t urow_hi = vld1q_u8(u_hi);                                       \
    for (; ap < ae; ap += sizeof(uint8x16_t), bp += sizeof(uint8x16_t)) {      \
      uint8x16_t bx = vld1q_u8(bp);                                            \
      uint8x16_t lo = vandq_u8(bx, mask);                                      \
      bx = vshrq_n_u8(bx, 4);                                                  \
      uint8x16_t hi = vandq_u8(bx, mask);                                      \
      lo = vqtbl1q_u8(urow_lo, lo);                                            \
      hi = vqtbl1q_u8(urow_hi, hi);                                            \
      vst1q_u8(ap, f(vld1q_u8(ap), veorq_u8(lo, hi)));                         \
    }                                                                          \
    op##_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(uint8x16_t));                   \
  } while (0)
#define OBL_SHUF_XOR veorq_u8

#else
#define OBL_SHUF(op, a, b, f)                                                  \
  do {                                                                         \
    op##_ref(a, b, u, k);                                                      \
  } while (0)
#define OBL_SHUF_XOR
#endif
#endif
#endif
#endif

#define OBL_NOOP(a, b) (b)
void obl_axpy(u8 *a, u8 *b, u8 u, int k) {
  OBL_SHUF(obl_axpy, a, b, OBL_SHUF_XOR);
}

void obl_scal(u8 *a, u8 u, int k) { OBL_SHUF(obl_scal, a, a, OBL_NOOP); }

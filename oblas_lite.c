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

static void axpy_ref(u8 *a, u8 *b, u8 u, int k) {
  register u8 *u_row = &GF2_8_MUL[u << 8];
  register u8 *ap = a, *ae = &a[k], *bp = b;
  for (; ap != ae; ap++, bp++)
    *ap ^= u_row[*bp];
}

static void scal_ref(u8 *a, u8 u, int k) {
  register u8 *u_row = &GF2_8_MUL[u << 8];
  register u8 *ap = a, *ae = &a[k];
  for (; ap != ae; ap++)
    *ap = u_row[*ap];
}

#if defined(__AVX2__)
#include <immintrin.h>

void obl_axpy(u8 *a, u8 *b, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;
  const __m256i mask = _mm256_set1_epi8(0x0f);
  const __m256i urow_lo = _mm256_loadu2_m128i((__m128i *)u_lo, (__m128i *)u_lo);
  const __m256i urow_hi = _mm256_loadu2_m128i((__m128i *)u_hi, (__m128i *)u_hi);
  __m256i *ap = (__m256i *)a, *ae = (__m256i *)(a + k - (k % sizeof(__m256i)));
  __m256i *bp = (__m256i *)b;
  for (; ap < ae; ap++, bp++) {
    __m256i bx = _mm256_loadu_si256(bp);
    __m256i lo = _mm256_and_si256(bx, mask);
    bx = _mm256_srli_epi64(bx, 4);
    __m256i hi = _mm256_and_si256(bx, mask);
    lo = _mm256_shuffle_epi8(urow_lo, lo);
    hi = _mm256_shuffle_epi8(urow_hi, hi);
    _mm256_storeu_si256(
        ap, _mm256_xor_si256(_mm256_loadu_si256(ap), _mm256_xor_si256(lo, hi)));
  }
  axpy_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(__m256i));
}

void obl_scal(u8 *a, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;
  const __m256i mask = _mm256_set1_epi8(0x0f);
  const __m256i urow_lo = _mm256_loadu2_m128i((__m128i *)u_lo, (__m128i *)u_lo);
  const __m256i urow_hi = _mm256_loadu2_m128i((__m128i *)u_hi, (__m128i *)u_hi);
  __m256i *ap = (__m256i *)a, *ae = (__m256i *)(a + k - (k % sizeof(__m256i)));
  for (; ap < ae; ap++) {
    __m256i ax = _mm256_loadu_si256(ap);
    __m256i lo = _mm256_and_si256(ax, mask);
    ax = _mm256_srli_epi64(ax, 4);
    __m256i hi = _mm256_and_si256(ax, mask);
    lo = _mm256_shuffle_epi8(urow_lo, lo);
    hi = _mm256_shuffle_epi8(urow_hi, hi);
    _mm256_storeu_si256(ap, _mm256_xor_si256(lo, hi));
  }
  scal_ref((u8 *)ap, u, k % sizeof(__m256i));
}
#else
#if defined(__SSSE3__) || (defined(_MSC_VER) && defined(_M_X64))
#include <emmintrin.h>
#include <tmmintrin.h>

void obl_axpy(u8 *a, u8 *b, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;
  const __m128i mask = _mm_set1_epi8(0x0f);
  const __m128i urow_lo = _mm_loadu_si128((__m128i *)u_lo);
  const __m128i urow_hi = _mm_loadu_si128((__m128i *)u_hi);
  __m128i *ap = (__m128i *)a, *ae = (__m128i *)(a + k);
  __m128i *bp = (__m128i *)b;
  for (; ap < ae; ap++, bp++) {
    __m128i bx = _mm_loadu_si128(bp);
    __m128i lo = _mm_and_si128(bx, mask);
    bx = _mm_srli_epi64(bx, 4);
    __m128i hi = _mm_and_si128(bx, mask);
    lo = _mm_shuffle_epi8(urow_lo, lo);
    hi = _mm_shuffle_epi8(urow_hi, hi);
    _mm_storeu_si128(ap,
                     _mm_xor_si128(_mm_loadu_si128(ap), _mm_xor_si128(lo, hi)));
  }
  axpy_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(__m128i));
}

void obl_scal(u8 *a, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;

  const __m128i mask = _mm_set1_epi8(0x0f);
  const __m128i urow_lo = _mm_loadu_si128((__m128i *)u_lo);
  const __m128i urow_hi = _mm_loadu_si128((__m128i *)u_hi);
  __m128i *ap = (__m128i *)a, *ae = (__m128i *)(a + k);
  for (; ap < ae; ap++) {
    __m128i ax = _mm_loadu_si128(ap);
    __m128i lo = _mm_and_si128(ax, mask);
    ax = _mm_srli_epi64(ax, 4);
    __m128i hi = _mm_and_si128(ax, mask);
    lo = _mm_shuffle_epi8(urow_lo, lo);
    hi = _mm_shuffle_epi8(urow_hi, hi);
    _mm_storeu_si128(ap, _mm_xor_si128(lo, hi));
  }
  scal_ref((u8 *)ap, u, k % sizeof(__m128i));
}
#else
#if defined(__GNUC__) && defined(__aarch64__)
#include <arm_neon.h>
/* FIXME: support armv7 with mfpu=neon and synthesized vqtbl1q_u8 */
void obl_axpy(u8 *a, u8 *b, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;
  u8 *ap = a, *ae = a + k, *bp = (u8 *)b;
  uint8x16_t mask = vdupq_n_u8(0x0f);
  uint8x16_t urow_lo = vld1q_u8(u_lo);
  uint8x16_t urow_hi = vld1q_u8(u_hi);
  for (; ap < ae; ap += sizeof(uint8x16_t), bp += sizeof(uint8x16_t)) {
    uint8x16_t bx = vld1q_u8(bp);
    uint8x16_t lo = vandq_u8(bx, mask);
    bx = vshrq_n_u8(bx, 4);
    uint8x16_t hi = vandq_u8(bx, mask);
    lo = vqtbl1q_u8(urow_lo, lo);
    hi = vqtbl1q_u8(urow_hi, hi);
    uint8x16_t ux = veorq_u8(lo, hi);
    uint8x16_t ax = vld1q_u8(ap);
    vst1q_u8(ap, veorq_u8(ux, ax));
  }
  axpy_ref((u8 *)ap, (u8 *)bp, u, k % sizeof(uint8x16_t));
}

void obl_scal(u8 *a, u8 u, int k) {
  const u8 *u_lo = GF2_8_SHUF_LO + u * 16;
  const u8 *u_hi = GF2_8_SHUF_HI + u * 16;
  u8 *ap = a, *ae = a + k;
  uint8x16_t mask = vdupq_n_u8(0x0f);
  uint8x16_t urow_lo = vld1q_u8(u_lo);
  uint8x16_t urow_hi = vld1q_u8(u_hi);
  for (; ap < ae; ap += sizeof(uint8x16_t)) {
    uint8x16_t ax = vld1q_u8(ap);
    uint8x16_t lo = vandq_u8(ax, mask);
    ax = vshrq_n_u8(ax, 4);
    uint8x16_t hi = vandq_u8(ax, mask);
    lo = vqtbl1q_u8(urow_lo, lo);
    hi = vqtbl1q_u8(urow_hi, hi);
    vst1q_u8(ap, veorq_u8(lo, hi));
  }
  scal_ref((u8 *)ap, u, k % sizeof(uint8x16_t));
}
#else
void obl_axpy(u8 *a, u8 *b, u8 u, int k) { axpy_ref(a, b, u, k); }
void obl_scal(u8 *a, u8 u, int k) { scal_ref(a, u, k); }
#endif
#endif
#endif

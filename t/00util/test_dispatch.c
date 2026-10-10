#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "oblas_lite.h"

#if defined(OBLAS_ARCH_X86) && defined(__GNUC__)
#include <immintrin.h>

enum { SSSE3 = 1, AVX2 = 2, AVX512F = 4, AVX512BW = 8, AVX512DQ = 16, AVX512VL = 32, GFNI = 64 };
static unsigned features;

static int supports(const char *name)
{
    const char *names[] = {"ssse3", "avx2", "avx512f", "avx512bw", "avx512dq", "avx512vl", "gfni"};
    for (unsigned i = 0; i < sizeof(names) / sizeof(names[0]); i++)
        if (strcmp(name, names[i]) == 0)
            return (features & (1U << i)) != 0;
    return 0;
}

#define __builtin_cpu_supports(name) supports(name)
#define __builtin_cpu_init() ((void)0)
#endif

/* Compile the real initializer and selector with synthetic CPU capabilities.
 * Selected vector functions are compared, never executed. */
#include "../../deps/obl/oblas_lite.c"
#include "../../rs.c"

int main(void)
{
    assert(!gf2_8_mul_initialized);
    reed_solomon_init();
    assert(gf2_8_mul_initialized);
    assert(GF2_8_MUL[0x0203] == 6);
    reed_solomon_init();
    assert(GF2_8_MUL[0x0203] == 6);

#if defined(OBLAS_ARCH_X86) && defined(__GNUC__)
    const unsigned avx512 = AVX512F | AVX512BW | AVX512DQ | AVX512VL;
    struct oblas_impl impl;
    oblas_get_impl(&impl);
    assert(impl.axpy == obl_axpy_ref && impl.align_size == sizeof(void *));

    features = GFNI;
    oblas_get_impl(&impl);
    assert(impl.axpy == obl_axpy_ref);

    for (unsigned gfni = 0; gfni <= GFNI; gfni += GFNI) {
        features = SSSE3 | gfni;
        oblas_get_impl(&impl);
        assert(impl.axpy == (gfni ? obl_axpy_ssse3_gfni : obl_axpy_ssse3_std));
        assert(impl.align_size == 16);

        features |= AVX2;
        oblas_get_impl(&impl);
        assert(impl.axpy == (gfni ? obl_axpy_avx2_gfni : obl_axpy_avx2_std));
        assert(impl.align_size == 32);

        for (unsigned missing = AVX512F; missing <= AVX512VL; missing <<= 1) {
            features = SSSE3 | AVX2 | (avx512 & ~missing) | gfni;
            oblas_get_impl(&impl);
            assert(impl.axpy == (gfni ? obl_axpy_avx2_gfni : obl_axpy_avx2_std));
            assert(impl.align_size == 32);
        }

        features = SSSE3 | AVX2 | avx512 | gfni;
        oblas_get_impl(&impl);
        assert(impl.axpy == (gfni ? obl_axpy_avx512_gfni : obl_axpy_avx512_std));
        assert(impl.align_size == 64);
    }
#endif
    puts("Initialization and dispatch tests passed");
    return 0;
}

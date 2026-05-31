#include <stdint.h>
#include <stddef.h>

#include "gf2_8_tables.h"

#if defined(OBLAS_AVX512)
#define OBLAS_ALIGN 64
#elif defined(OBLAS_AVX2)
#define OBLAS_ALIGN 32
#elif defined(OBLAS_SSE3) || defined(OBLAS_NEON)
#define OBLAS_ALIGN 16
#else
#define OBLAS_ALIGN 1
#endif

typedef uint8_t u8;
typedef uint32_t u32;

void obl_axpy(u8 *a, u8 *b, u8 u, unsigned k);
void obl_scal(u8 *a, u8 u, unsigned k);
void obl_scale_copy(u8 *a, u8 *b, u8 u, unsigned k);
void obl_swap(u8 *a, u8 *b, unsigned k);
void obl_axpyb32(u8 *a, u32 *b, u8 u, unsigned k);

void *obl_alloc(size_t num_rows, size_t row_size);
void obl_free(void *ptr);

#ifndef OBLAS_COMMON_H
#define OBLAS_COMMON_H

#include <stddef.h>
#include <stdint.h>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define OBLAS_ARCH_X86 1
#elif defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__) || defined(_M_ARM)
#define OBLAS_ARCH_ARM 1
#endif

#ifdef __cplusplus
extern "C" {
#endif

void obl_swap(uint8_t *a, uint8_t *b, unsigned k);
void *obl_alloc(size_t num_rows, size_t row_size, size_t alignment);
void obl_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* OBLAS_COMMON_H */

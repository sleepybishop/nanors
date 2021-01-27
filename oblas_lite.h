#include <stdint.h>

#include "gf2_8_tables.h"

typedef uint8_t u8;

void obl_fill_mul_tab(void);
void obl_axpy(u8 *a, u8 *b, u8 u, int k);
void obl_scal(u8 *a, u8 u, int k);

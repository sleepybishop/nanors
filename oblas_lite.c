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

void obl_axpy(u8 *a, u8 *b, u8 u, int k) { axpy_ref(a, b, u, k); }

void obl_scal(u8 *a, u8 u, int k) { scal_ref(a, u, k); }

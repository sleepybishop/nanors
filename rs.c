#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "oblas_lite.h"
#include "rs.h"

static void axpy(u8 *a, u8 *b, u8 u, int k) {
  if (u == 0)
    return;

  if (u == 1) {
    for (int i = 0; i < k; i++)
      a[i] = a[i] ^ b[i];
  } else {
    obl_axpy(a, b, u, k);
  }
}

static void scal(u8 *a, u8 u, int k) {
  if (u < 2)
    return;
  obl_scal(a, u, k);
}

static void gemm(u8 *a, u8 **b, u8 **c, int n, int k, int m) {
  int ci = 0;
  for (int row = 0; row < n; row++, ci++) {
    u8 *ap = a + (row * k);
    memset(c[ci], 0, m);
    for (int idx = 0; idx < k; idx++)
      axpy(c[ci], b[idx], ap[idx], m);
  }
}

static int invert_mat(u8 *src, u8 **dst, int V0, int K, int T, int *c) {
  int V0b = V0, W = K - V0;
  u8 u = 0;
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < W; j++) {
      src[i * W + j] = src[(V0 + i) * K + c[V0 + j]];
    }
  }
  for (; V0 < K; V0++) {
    for (int row = 0; row < V0b; row++) {
      u = src[V0 * K + c[row]];
      axpy(dst[c[V0]], dst[c[row]], u, T);
    }
  }
  for (int x = 0; x < W; x++) {
    u = GF2_8_INV[src[x * W + x]];
    scal(src + x * W + x, u, W);
    scal(dst[c[V0b + x]], u, T);
    for (int row = x + 1; row < W; row++) {
      u = src[row * W + x];
      axpy(src + row * W, src + x * W, u, W);
      axpy(dst[c[V0b + row]], dst[c[V0b + x]], u, T);
    }
  }
  for (int x = W - 1; x >= 0; x--) {
    for (int row = 0; row < x; row++) {
      u = src[row * W + x];
      axpy(dst[c[V0b + row]], dst[c[V0b + x]], u, T);
    }
  }
  return 0;
}

void reed_solomon_init() { obl_fill_mul_tab(); }

reed_solomon *reed_solomon_new(int ds, int ps) {
  reed_solomon *rs = NULL;
  if ((ds + ps) > DATA_SHARDS_MAX || ds <= 0 || ps <= 0 || ps > ds)
    return NULL;
  rs = calloc(1, sizeof(reed_solomon) + 3 * ps * ds);
  if (!rs)
    return NULL;
  rs->ds = ds;
  rs->ps = ps;
  rs->ts = ds + ps;

  for (int j = 0; j < rs->ps; j++) {
    u8 *row = rs->p + j * rs->ds;
    for (int i = 0; i < rs->ds; i++)
      row[i] = GF2_8_INV[(rs->ps + i) ^ j];
  }
  return rs;
}

void reed_solomon_release(reed_solomon *rs) {
  if (rs)
    free(rs);
}

int reed_solomon_decode(reed_solomon *rs, u8 **data, u8 *marks, int nr_shards,
                        int bs) {
  if (nr_shards < rs->ts)
    return -1;

  u8 src[DATA_SHARDS_MAX * DATA_SHARDS_MAX];
  int erased_blocks[DATA_SHARDS_MAX], colperm[DATA_SHARDS_MAX], gaps = 0;
  memset(src, 0, rs->ds * rs->ds);

  for (int i = 0; i < rs->ds; i++) {
    if (marks[i])
      erased_blocks[gaps++] = i;
  }
  for (int i = 0, j = 0; i < rs->ds - gaps; i++, j++) {
    while (marks[j])
      j++;
    colperm[i] = j;
  }
  for (int i = 0, j = rs->ds - gaps; i < gaps; i++, j++)
    colperm[j] = erased_blocks[i];
  for (int i = 0; i < rs->ds - gaps; i++)
    src[i * rs->ds + colperm[i]] = 1;

  int i = rs->ds - gaps;
  for (int j = rs->ds, g = 0, l = 0; i < rs->ds; i++, j++, l++) {
    while (marks[j]) {
      j++;
      l++;
    }
    if (j >= nr_shards)
      break;
    memcpy(src + i * rs->ds, rs->p + l * rs->ds, rs->ds);
    memcpy(data[erased_blocks[g++]], data[j], bs);
  }
  if (i < rs->ds)
    return -1;
  invert_mat(src, data, rs->ds - gaps, rs->ds, bs, colperm);
  return 0;
}

int reed_solomon_encode(reed_solomon *rs, u8 **shards, int nr_shards, int bs) {
  if (nr_shards < rs->ts)
    return -1;
  gemm(rs->p, shards, shards + rs->ds, rs->ps, rs->ds, bs);
  return 0;
}

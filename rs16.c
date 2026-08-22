#include "rs16.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "deps/obl/oblas16.h"

static int invert_mat_gf16(uint16_t *src, uint16_t *dst, int ds)
{
    for (int i = 0; i < ds; i++) {
        for (int j = 0; j < ds; j++) {
            dst[(size_t)i * ds + j] = (i == j) ? 1 : 0;
        }
    }

    for (int i = 0; i < ds; i++) {
        if (src[(size_t)i * ds + i] == 0) {
            int pivot = -1;
            for (int j = i + 1; j < ds; j++) {
                if (src[(size_t)j * ds + i] != 0) {
                    pivot = j;
                    break;
                }
            }
            if (pivot == -1)
                return -1;
            for (int k = 0; k < ds; k++) {
                uint16_t tmp = src[(size_t)i * ds + k];
                src[(size_t)i * ds + k] = src[(size_t)pivot * ds + k];
                src[(size_t)pivot * ds + k] = tmp;

                tmp = dst[(size_t)i * ds + k];
                dst[(size_t)i * ds + k] = dst[(size_t)pivot * ds + k];
                dst[(size_t)pivot * ds + k] = tmp;
            }
        }

        uint16_t pivot_val = src[(size_t)i * ds + i];
        uint16_t inv_p = gf16_inv(pivot_val);
        for (int k = 0; k < ds; k++) {
            src[(size_t)i * ds + k] = gf16_mul(src[(size_t)i * ds + k], inv_p);
            dst[(size_t)i * ds + k] = gf16_mul(dst[(size_t)i * ds + k], inv_p);
        }

        for (int j = 0; j < ds; j++) {
            if (i == j)
                continue;
            uint16_t factor = src[(size_t)j * ds + i];
            if (factor != 0) {
                for (int k = 0; k < ds; k++) {
                    src[(size_t)j * ds + k] ^= gf16_mul(factor, src[(size_t)i * ds + k]);
                    dst[(size_t)j * ds + k] ^= gf16_mul(factor, dst[(size_t)i * ds + k]);
                }
            }
        }
    }
    return 0;
}

void reed_solomon16_init(void)
{
    oblas16_init();
}

reed_solomon16 *reed_solomon16_new(int data_shards, int parity_shards)
{
    if (data_shards <= 0 || data_shards > RS16_DATA_SHARDS_MAX || parity_shards <= 0 || parity_shards > RS16_DATA_SHARDS_MAX ||
        (uint32_t)data_shards + (uint32_t)parity_shards > 65536U)
        return NULL;

    reed_solomon16_init();

    struct oblas16_impl impl;
    oblas16_get_impl(&impl);

    reed_solomon16 *rs = (reed_solomon16 *)calloc(1, sizeof(reed_solomon16));
    if (!rs)
        return NULL;
    rs->ds = data_shards;
    rs->ps = parity_shards;
    rs->ts = data_shards + parity_shards;
    rs->axpy = impl.axpy;
    rs->scal = impl.scal;
    rs->axiy = impl.axiy;
    rs->align_size = impl.align_size;

    size_t matrix_count = (size_t)rs->ps * (size_t)rs->ds;
    if (matrix_count > SIZE_MAX / sizeof(uint16_t)) {
        free(rs);
        return NULL;
    }
    rs->p = (uint16_t *)calloc(matrix_count, sizeof(uint16_t));
    if (!rs->p) {
        free(rs);
        return NULL;
    }
    for (int j = 0; j < rs->ps; j++) {
        for (int i = 0; i < rs->ds; i++) {
            rs->p[(size_t)j * rs->ds + i] = gf16_inv((uint16_t)((rs->ps + i) ^ j));
        }
    }

    return rs;
}

void reed_solomon16_release(reed_solomon16 *rs)
{
    if (rs) {
        if (rs->p)
            free(rs->p);
        free(rs);
    }
}

int reed_solomon16_encode(reed_solomon16 *rs, uint16_t **shards, int nr_shards, int bs)
{
    if (!rs || !shards || nr_shards < rs->ts || bs <= 0)
        return -1;

    for (int i = 0; i < rs->ps; i++) {
        memset(shards[rs->ds + i], 0, bs * sizeof(uint16_t));
        for (int j = 0; j < rs->ds; j++) {
            uint16_t coef = rs->p[(size_t)i * rs->ds + j];
            if (coef == 0)
                continue;
            uint16_t *dst = shards[rs->ds + i];
            uint16_t *src = shards[j];
            rs->axpy(dst, src, coef, bs);
        }
    }
    return 0;
}

int reed_solomon16_decode(reed_solomon16 *rs, uint16_t **shards, uint8_t *marks, int nr_shards, int bs)
{
    if (!rs || !shards || !marks || nr_shards < rs->ts || bs <= 0) {
        return -1;
    }

    int gaps = 0;
    for (int i = 0; i < rs->ds; i++)
        gaps += marks[i] != 0;

    if (gaps > 0) {
        uint16_t *erasures = (uint16_t *)calloc((size_t)gaps, sizeof(uint16_t));
        uint16_t *surviving = (uint16_t *)calloc((size_t)rs->ds, sizeof(uint16_t));
        if (!erasures || !surviving) {
            free(erasures);
            free(surviving);
            return -1;
        }

        int gap_idx = 0;
        for (int i = 0; i < rs->ds; i++) {
            if (marks[i])
                erasures[gap_idx++] = (uint16_t)i;
        }

        int surv_idx = 0;
        for (int i = 0; i < rs->ts && surv_idx < rs->ds; i++) {
            if (!marks[i])
                surviving[surv_idx++] = (uint16_t)i;
        }
        if (surv_idx < rs->ds) {
            free(erasures);
            free(surviving);
            return -1;
        }

        size_t matrix_count = (size_t)rs->ds * (size_t)rs->ds;
        if (matrix_count > SIZE_MAX / sizeof(uint16_t)) {
            free(erasures);
            free(surviving);
            return -1;
        }
        uint16_t *decode_mat = (uint16_t *)calloc(matrix_count, sizeof(uint16_t));
        uint16_t *decode_mat_inv = (uint16_t *)calloc(matrix_count, sizeof(uint16_t));
        if (!decode_mat || !decode_mat_inv) {
            free(decode_mat);
            free(decode_mat_inv);
            free(erasures);
            free(surviving);
            return -1;
        }

        for (int i = 0; i < rs->ds; i++) {
            int r = surviving[i];
            if (r < rs->ds) {
                decode_mat[(size_t)i * rs->ds + r] = 1;
            } else {
                int p_idx = r - rs->ds;
                memcpy(&decode_mat[(size_t)i * rs->ds], &rs->p[(size_t)p_idx * rs->ds], (size_t)rs->ds * sizeof(uint16_t));
            }
        }

        if (invert_mat_gf16(decode_mat, decode_mat_inv, rs->ds) != 0) {
            free(decode_mat);
            free(decode_mat_inv);
            free(erasures);
            free(surviving);
            return -1;
        }

        for (int e = 0; e < gaps; e++) {
            int erased = erasures[e];
            uint16_t *dst = shards[erased];
            memset(dst, 0, (size_t)bs * sizeof(uint16_t));
            for (int j = 0; j < rs->ds; j++) {
                uint16_t coef = decode_mat_inv[(size_t)erased * rs->ds + j];
                if (coef != 0)
                    rs->axpy(dst, shards[surviving[j]], coef, (unsigned)bs);
            }
        }

        free(decode_mat);
        free(decode_mat_inv);
        free(erasures);
        free(surviving);
    }

    for (int i = 0; i < rs->ps; i++) {
        if (marks[rs->ds + i]) {
            memset(shards[rs->ds + i], 0, bs * sizeof(uint16_t));
            for (int j = 0; j < rs->ds; j++) {
                uint16_t coef = rs->p[(size_t)i * rs->ds + j];
                if (coef == 0)
                    continue;
                uint16_t *dst = shards[rs->ds + i];
                uint16_t *src = shards[j];
                rs->axpy(dst, src, coef, bs);
            }
        }
    }

    return 0;
}

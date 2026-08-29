#include "rs_rlrs.h"
#include "rs16_afft.h"
#include "deps/obl/oblas16.h"
#include "deps/obl/oblas16_afft.h"
#include <stdlib.h>
#include <string.h>

struct reed_solomon_rlrs {
    int ds;
    int K_prime;
    int log_K_prime;
    int cached_N;
    void *flat_alloc;
    uint16_t *log_walsh;
    struct afft_workspace ws;

    /* Cache management */
    uint16_t **cached_shards;
    int cached_block_size;
    uint16_t *coeff_buf;
    size_t coeff_buf_size;
    int has_cache;

    /* Erasure locator cache for repeated packet-loss patterns. */
    uint8_t *cached_decode_marks;
    size_t cached_decode_marks_size;
    int *cached_decode_parity_indices;
    int cached_decode_parity_capacity;
    int cached_decode_parity_count;
    int cached_decode_N;
    int has_decode_cache;
};

static int next_power_of_2(int n)
{
    int p = 1;
    while (p < n)
        p <<= 1;
    return p;
}

static int log2_int(int n)
{
    int p = 0;
    while ((1 << p) < n)
        p++;
    return p;
}

static int rlrs_ensure_workspace(reed_solomon_rlrs *rs, int N)
{
    if (N <= rs->cached_N)
        return 0;

    if (rs->flat_alloc) {
        obl_free(rs->flat_alloc);
        rs->flat_alloc = NULL;
    }

    int log_N = log2_int(N);
    if (log_N < 1)
        log_N = 1;
    size_t sz_work_buf = (size_t)N * BATCH_SIZE * sizeof(uint16_t);
    size_t sz_needed_buf = log_N * (1 << (log_N - 1)) * sizeof(uint8_t);
    size_t sz_L_eval = N * sizeof(uint16_t);
    size_t sz_inv_Lp_eval = N * sizeof(uint16_t);
    size_t sz_is_erased = N * sizeof(uint8_t);
    size_t sz_error_locations = N * sizeof(uint16_t);
    size_t sz_log_walsh = N * sizeof(uint16_t);

    sz_work_buf = (sz_work_buf + 63) & ~63;
    sz_needed_buf = (sz_needed_buf + 63) & ~63;
    sz_L_eval = (sz_L_eval + 63) & ~63;
    sz_inv_Lp_eval = (sz_inv_Lp_eval + 63) & ~63;
    sz_is_erased = (sz_is_erased + 63) & ~63;
    sz_error_locations = (sz_error_locations + 63) & ~63;
    sz_log_walsh = (sz_log_walsh + 63) & ~63;

    size_t total_size = sz_work_buf + sz_needed_buf + sz_L_eval + sz_inv_Lp_eval + sz_is_erased + sz_error_locations + sz_log_walsh;

    rs->flat_alloc = obl_alloc(1, total_size, 64);
    if (!rs->flat_alloc) {
        rs->cached_N = 0;
        return -1;
    }

    uint8_t *ptr = (uint8_t *)rs->flat_alloc;
    rs->ws.chunk_buf = (uint16_t *)ptr;
    rs->ws.buf = rs->ws.chunk_buf;
    rs->ws.acc = NULL;
    ptr += sz_work_buf;
    rs->ws.needed_buf = ptr;
    ptr += sz_needed_buf;
    rs->ws.L_eval = (uint16_t *)ptr;
    ptr += sz_L_eval;
    rs->ws.inv_Lp_eval = (uint16_t *)ptr;
    ptr += sz_inv_Lp_eval;
    rs->ws.is_erased = ptr;
    ptr += sz_is_erased;
    rs->ws.error_locations = (uint16_t *)ptr;
    ptr += sz_error_locations;
    rs->log_walsh = (uint16_t *)ptr;
    rs->ws.erasures = NULL;
    rs->ws.E = NULL;

    if (reed_solomon16_afft_build_log_walsh(rs->log_walsh, N) < 0) {
        obl_free(rs->flat_alloc);
        rs->flat_alloc = NULL;
        rs->cached_N = 0;
        return -1;
    }

    rs->cached_N = N;
    rs->has_decode_cache = 0;
    return 0;
}

static int rlrs_cache_matches(reed_solomon_rlrs *rs, uint16_t **data_shards, int block_size)
{
    if (!rs->has_cache || rs->cached_block_size != block_size)
        return 0;
    for (int i = 0; i < rs->ds; i++) {
        if (rs->cached_shards[i] != data_shards[i])
            return 0;
    }
    return 1;
}

static int rlrs_decode_cache_matches(reed_solomon_rlrs *rs, const uint8_t *marks, const int *parity_indices, int parity_count,
                                     int N)
{
    size_t marks_size = (size_t)rs->ds + parity_count;
    return rs->has_decode_cache && rs->cached_decode_N == N && rs->cached_decode_parity_count == parity_count &&
           memcmp(rs->cached_decode_marks, marks, marks_size) == 0 &&
           memcmp(rs->cached_decode_parity_indices, parity_indices, (size_t)parity_count * sizeof(int)) == 0;
}

static void rlrs_update_decode_cache(reed_solomon_rlrs *rs, const uint8_t *marks, const int *parity_indices, int parity_count,
                                     int N)
{
    size_t marks_size = (size_t)rs->ds + parity_count;
    if (marks_size > rs->cached_decode_marks_size) {
        uint8_t *new_marks = (uint8_t *)realloc(rs->cached_decode_marks, marks_size);
        if (!new_marks)
            return;
        rs->cached_decode_marks = new_marks;
        rs->cached_decode_marks_size = marks_size;
    }
    if (parity_count > rs->cached_decode_parity_capacity) {
        int *new_indices = (int *)realloc(rs->cached_decode_parity_indices, (size_t)parity_count * sizeof(int));
        if (!new_indices)
            return;
        rs->cached_decode_parity_indices = new_indices;
        rs->cached_decode_parity_capacity = parity_count;
    }

    memcpy(rs->cached_decode_marks, marks, marks_size);
    memcpy(rs->cached_decode_parity_indices, parity_indices, (size_t)parity_count * sizeof(int));
    rs->cached_decode_parity_count = parity_count;
    rs->cached_decode_N = N;
    rs->has_decode_cache = 1;
}

void reed_solomon_rlrs_invalidate_cache(reed_solomon_rlrs *rs)
{
    if (rs) {
        rs->has_cache = 0;
    }
}

reed_solomon_rlrs *reed_solomon_rlrs_new(int data_shards)
{
    if (data_shards <= 0 || data_shards > RS16_DATA_SHARDS_MAX)
        return NULL;

    int K_prime = next_power_of_2(data_shards);
    if (K_prime >= 65536)
        return NULL;

    reed_solomon16_afft_init();

    reed_solomon_rlrs *rs = (reed_solomon_rlrs *)calloc(1, sizeof(reed_solomon_rlrs));
    if (!rs)
        return NULL;

    rs->ds = data_shards;
    rs->K_prime = K_prime;
    rs->log_K_prime = log2_int(rs->K_prime);
    rs->cached_N = 0;
    rs->flat_alloc = NULL;

    rs->cached_shards = (uint16_t **)calloc(data_shards, sizeof(uint16_t *));
    rs->cached_block_size = 0;
    rs->coeff_buf = NULL;
    rs->coeff_buf_size = 0;
    rs->has_cache = 0;

    oblas16_get_impl(&rs->ws.o16);
    oblas16_afft_get_impl(&rs->ws.afft);
    rs->ws.axpy = rs->ws.o16.axpy;
    rs->ws.axiy = rs->ws.o16.axiy;

    return rs;
}

void reed_solomon_rlrs_release(reed_solomon_rlrs *rs)
{
    if (rs) {
        if (rs->flat_alloc) {
            obl_free(rs->flat_alloc);
        }
        if (rs->cached_shards) {
            free(rs->cached_shards);
        }
        if (rs->coeff_buf) {
            free(rs->coeff_buf);
        }
        free(rs->cached_decode_marks);
        free(rs->cached_decode_parity_indices);
        free(rs);
    }
}

int reed_solomon_rlrs_encode_many(reed_solomon_rlrs *rs, uint16_t **data_shards, const int *parity_indices, uint16_t **out_parity,
                                  int parity_count, int block_size)
{
    if (!rs || !data_shards || !parity_indices || !out_parity || parity_count <= 0 || block_size <= 0)
        return -1;

    int K = rs->ds;
    int K_prime = rs->K_prime;
    int max_parity_index = -1;
    for (int i = 0; i < parity_count; i++) {
        if (parity_indices[i] < 0 || parity_indices[i] >= 65536 - K_prime || !out_parity[i])
            return -1;
        if (parity_indices[i] > max_parity_index)
            max_parity_index = parity_indices[i];
    }

    int N = next_power_of_2(K_prime + max_parity_index + 1);
    if (N > 65536)
        return -1;

    if (rlrs_ensure_workspace(rs, N) < 0)
        return -1;

    int log_N = log2_int(N);
    int log_K_prime = rs->log_K_prime;

    struct afft_workspace *ws = &rs->ws;
    uint16_t *chunk_buf = ws->chunk_buf;
    uint8_t *needed_buf = ws->needed_buf;

    if (log_N > 0) {
        memset(needed_buf, 0, log_N * (1 << (log_N - 1)) * sizeof(uint8_t));
    }

    uint8_t *needed[16];
    int offset = 0;
    for (int k = 0; k < log_N; k++) {
        needed[k] = needed_buf + offset;
        offset += 1 << (log_N - 1 - k);
    }

    for (int i = 0; i < parity_count; i++) {
        int idx = K_prime + parity_indices[i];
        for (int k = 0; k < log_N; k++) {
            needed[k][idx >> (k + 1)] = 1;
        }
    }

    /* Check if cache matches */
    int cache_hit = rlrs_cache_matches(rs, data_shards, block_size);
    if (!cache_hit) {
        size_t required_coeff_size = (size_t)K_prime * block_size;
        if (required_coeff_size > rs->coeff_buf_size) {
            uint16_t *new_buf = (uint16_t *)realloc(rs->coeff_buf, required_coeff_size * sizeof(uint16_t));
            if (!new_buf)
                return -1;
            rs->coeff_buf = new_buf;
            rs->coeff_buf_size = required_coeff_size;
        }
        rs->cached_block_size = block_size;
        for (int i = 0; i < K; i++) {
            rs->cached_shards[i] = data_shards[i];
        }
    }

    for (int start = 0; start < block_size; start += BATCH_SIZE) {
        int current_chunk = (block_size - start < BATCH_SIZE) ? (block_size - start) : BATCH_SIZE;

        memset(chunk_buf, 0, N * BATCH_SIZE * sizeof(uint16_t));

        if (cache_hit) {
            /* Restore precomputed IFFT coefficients directly into chunk_buf */
            if (current_chunk == BATCH_SIZE) {
                memcpy(chunk_buf, &rs->coeff_buf[(size_t)start * K_prime], (size_t)K_prime * BATCH_SIZE * sizeof(uint16_t));
            } else {
                for (int i = 0; i < K_prime; i++) {
                    memcpy(&chunk_buf[i * BATCH_SIZE], &rs->coeff_buf[(size_t)start * K_prime + i * current_chunk],
                           current_chunk * sizeof(uint16_t));
                }
            }
        } else {
            /* Copy data shards and run IFFT */
            for (int i = 0; i < K; i++) {
                memcpy(&chunk_buf[i * BATCH_SIZE], &data_shards[i][start], current_chunk * sizeof(uint16_t));
            }

            oblas16_afft_ifft(chunk_buf, log_K_prime, BATCH_SIZE, K, 0, &ws->o16, &ws->afft);

            /* Save computed IFFT coefficients to cache */
            if (current_chunk == BATCH_SIZE) {
                memcpy(&rs->coeff_buf[(size_t)start * K_prime], chunk_buf, (size_t)K_prime * BATCH_SIZE * sizeof(uint16_t));
            } else {
                for (int i = 0; i < K_prime; i++) {
                    memcpy(&rs->coeff_buf[(size_t)start * K_prime + i * current_chunk], &chunk_buf[i * BATCH_SIZE],
                           current_chunk * sizeof(uint16_t));
                }
            }
        }

        /* Apply one partial FFT for the union of all requested parity blocks. */
        oblas16_afft_fft(chunk_buf, log_N, BATCH_SIZE, needed, 0, &ws->o16, &ws->afft);

        for (int i = 0; i < parity_count; i++) {
            int idx = K_prime + parity_indices[i];
            memcpy(&out_parity[i][start], &chunk_buf[idx * BATCH_SIZE], current_chunk * sizeof(uint16_t));
        }
    }

    if (!cache_hit) {
        rs->has_cache = 1;
    }

    return 0;
}

int reed_solomon_rlrs_encode_block(reed_solomon_rlrs *rs, uint16_t **data_shards, int parity_index, uint16_t *out_parity,
                                   int block_size)
{
    int parity_indices[1] = {parity_index};
    uint16_t *out_parities[1] = {out_parity};
    return reed_solomon_rlrs_encode_many(rs, data_shards, parity_indices, out_parities, 1, block_size);
}

int reed_solomon_rlrs_decode(reed_solomon_rlrs *rs, uint16_t **shards, uint8_t *marks, int *received_parity_indices,
                             int received_parity_count, int block_size)
{
    if (!rs || !shards || !marks || !received_parity_indices || received_parity_count < 0 || block_size <= 0)
        return -1;

    int K = rs->ds;
    int K_prime = rs->K_prime;

    for (int i = 0; i < received_parity_count; i++) {
        if (received_parity_indices[i] < 0 || received_parity_indices[i] >= 65536 - K_prime)
            return -1;
    }

    int num_erasures = 0;
    int erasures[RS16_DATA_SHARDS_MAX];
    for (int i = 0; i < K; i++) {
        if (marks[i]) {
            erasures[num_erasures++] = i;
        }
    }

    if (num_erasures == 0)
        return 0;

    int num_surviving_parities = 0;
    for (int i = 0; i < received_parity_count; i++) {
        if (!marks[K + i]) {
            num_surviving_parities++;
        }
    }

    if (num_surviving_parities < num_erasures)
        return -1;

    int max_parity_index = -1;
    for (int i = 0; i < received_parity_count; i++) {
        if (!marks[K + i] && received_parity_indices[i] > max_parity_index) {
            max_parity_index = received_parity_indices[i];
        }
    }

    int N = next_power_of_2(K_prime + max_parity_index + 1);
    if (N > 65536)
        return -1;

    if (rlrs_ensure_workspace(rs, N) < 0)
        return -1;

    int log_N = log2_int(N);

    struct afft_workspace *ws = &rs->ws;

    if (!rlrs_decode_cache_matches(rs, marks, received_parity_indices, received_parity_count, N)) {
        uint8_t *is_erased = ws->is_erased;
        memset(is_erased, 1, N * sizeof(uint8_t));

        for (int i = 0; i < K; i++) {
            if (!marks[i]) {
                is_erased[i] = 0;
            }
        }
        for (int i = K; i < K_prime; i++) {
            is_erased[i] = 0;
        }
        for (int i = 0; i < received_parity_count; i++) {
            if (!marks[K + i]) {
                int idx = K_prime + received_parity_indices[i];
                if (idx < N) {
                    is_erased[idx] = 0;
                }
            }
        }

        uint16_t *error_locations = ws->error_locations;
        memset(error_locations, 0, N * sizeof(uint16_t));

        for (int i = 0; i < N; i++) {
            if (is_erased[i]) {
                error_locations[i] = 1;
            }
        }

        fwht_mod(error_locations, N);

        for (int i = 0; i < N; i++) {
            uint32_t x = (uint32_t)error_locations[i] * rs->log_walsh[i];
            uint32_t s = (x >> 16) + (x & 65535);
            s = (s >> 16) + (s & 65535);
            error_locations[i] = (s >= 65535) ? s - 65535 : s;
        }

        fwht_mod(error_locations, N);

        uint16_t *L_eval = ws->L_eval;
        uint16_t *inv_Lp_eval = ws->inv_Lp_eval;
        for (int i = 0; i < N; i++) {
            if (!is_erased[i]) {
                L_eval[i] = GF16_EXP[error_locations[i]];
            } else {
                inv_Lp_eval[i] = GF16_EXP[65535 - error_locations[i]];
            }
        }

        rlrs_update_decode_cache(rs, marks, received_parity_indices, received_parity_count, N);
    }

    uint16_t *L_eval = ws->L_eval;
    uint16_t *inv_Lp_eval = ws->inv_Lp_eval;

    uint8_t *needed_buf = ws->needed_buf;
    if (log_N > 0) {
        memset(needed_buf, 0, log_N * (1 << (log_N - 1)) * sizeof(uint8_t));
    }

    uint8_t *needed[16];
    int offset = 0;
    for (int k = 0; k < log_N; k++) {
        needed[k] = needed_buf + offset;
        offset += 1 << (log_N - 1 - k);
    }
    for (int i = 0; i < num_erasures; i++) {
        int idx = erasures[i];
        for (int k = 0; k < log_N; k++) {
            needed[k][idx >> (k + 1)] = 1;
        }
    }

    uint16_t *buf = ws->buf;
    for (int byte_idx = 0; byte_idx < block_size; byte_idx += BATCH_SIZE) {
        int batch = (block_size - byte_idx < BATCH_SIZE) ? (block_size - byte_idx) : BATCH_SIZE;

        memset(buf, 0, N * BATCH_SIZE * sizeof(uint16_t));

        for (int i = 0; i < K; i++) {
            if (!marks[i]) {
                uint16_t L = L_eval[i];
                uint16_t *shard = shards[i];
                ws->axiy(&buf[i * BATCH_SIZE], &shard[byte_idx], L, batch);
            }
        }

        for (int i = 0; i < received_parity_count; i++) {
            if (!marks[K + i]) {
                int idx = K_prime + received_parity_indices[i];
                uint16_t L = L_eval[idx];
                uint16_t *shard = shards[K + i];
                ws->axiy(&buf[idx * BATCH_SIZE], &shard[byte_idx], L, batch);
            }
        }

        oblas16_afft_ifft(buf, log_N, BATCH_SIZE, N, 0, &ws->o16, &ws->afft);

        for (int i = 1; i < N; i++) {
            int width = ((i ^ (i - 1)) + 1) >> 1;
            int dst_idx = (i - width) * BATCH_SIZE;
            int src_idx = i * BATCH_SIZE;
            ws->axpy(&buf[dst_idx], &buf[src_idx], 1, width * BATCH_SIZE);
        }

        oblas16_afft_fft(buf, log_N, BATCH_SIZE, needed, 0, &ws->o16, &ws->afft);

        for (int i = 0; i < num_erasures; i++) {
            int idx = erasures[i];
            uint16_t inv_L = inv_Lp_eval[idx];
            uint16_t *shard = shards[idx];
            memset(&shard[byte_idx], 0, batch * sizeof(uint16_t));
            ws->axiy(&shard[byte_idx], &buf[idx * BATCH_SIZE], inv_L, batch);
        }
    }

    return 0;
}

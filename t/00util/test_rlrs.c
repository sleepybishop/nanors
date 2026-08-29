#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#include "rs_rlrs.h"
#include "rs16_afft.h"

static uint16_t mod65535(uint32_t x)
{
    uint32_t s = (x >> 16) + (x & 65535);
    s = (s >> 16) + (s & 65535);
    return (uint16_t)(s >= 65535 ? s - 65535 : s);
}

static int test_subspace_log_walsh(void)
{
    uint16_t *full = calloc(65536, sizeof(uint16_t));
    uint16_t *subspace = calloc(65536, sizeof(uint16_t));
    uint16_t *kernel = calloc(65536, sizeof(uint16_t));
    assert(full != NULL && subspace != NULL && kernel != NULL);

    int status = 0;
    int sizes[] = {8, 128, 1024, 8192};
    for (size_t nidx = 0; nidx < sizeof(sizes) / sizeof(sizes[0]); nidx++) {
        int N = sizes[nidx];
        memset(full, 0, 65536 * sizeof(uint16_t));
        memset(subspace, 0, N * sizeof(uint16_t));
        for (int i = 0; i < N; i++) {
            uint16_t erased = (i % 3 == 0) || (i >= N / 2 && i % 7 != 0);
            full[i] = erased;
            subspace[i] = erased;
        }

        assert(reed_solomon16_afft_build_log_walsh(kernel, N) == 0);
        fwht_mod(full, 65536);
        for (int i = 0; i < 65536; i++)
            full[i] = mod65535((uint32_t)full[i] * LogWalsh[i]);
        fwht_mod(full, 65536);

        fwht_mod(subspace, N);
        for (int i = 0; i < N; i++)
            subspace[i] = mod65535((uint32_t)subspace[i] * kernel[i]);
        fwht_mod(subspace, N);

        if (memcmp(full, subspace, N * sizeof(uint16_t)) != 0) {
            printf("Subspace locator transform mismatch at N=%d\n", N);
            status = 1;
            break;
        }
    }

    free(full);
    free(subspace);
    free(kernel);
    return status;
}

static int verify_buffers(uint16_t **orig, uint16_t **reconstructed, int K, int T)
{
    for (int i = 0; i < K; i++) {
        assert(orig[i] != NULL);
        assert(reconstructed[i] != NULL);
        for (int j = 0; j < T; j++) {
            if (orig[i][j] != reconstructed[i][j]) {
                printf("Mismatch at shard %d, word %d: expected %04x, got %04x\n", i, j, orig[i][j], reconstructed[i][j]);
                return 0;
            }
        }
    }
    return 1;
}

static int test_bulk_matches_single(void)
{
    enum { K = 17, P = 9, T = 777 };
    int parity_indices[P] = {0, 7, 2, 31, 1, 15, 3, 23, 8};
    uint16_t *data_shards[K];
    uint16_t *single_parity[P];
    uint16_t *bulk_parity[P];

    srand(104);
    for (int i = 0; i < K; i++) {
        data_shards[i] = malloc(T * sizeof(uint16_t));
        assert(data_shards[i] != NULL);
        for (int j = 0; j < T; j++)
            data_shards[i][j] = rand() % 65536;
    }
    for (int i = 0; i < P; i++) {
        single_parity[i] = malloc(T * sizeof(uint16_t));
        bulk_parity[i] = malloc(T * sizeof(uint16_t));
        assert(single_parity[i] != NULL && bulk_parity[i] != NULL);
    }

    reed_solomon_rlrs *single = reed_solomon_rlrs_new(K);
    reed_solomon_rlrs *bulk = reed_solomon_rlrs_new(K);
    assert(single != NULL && bulk != NULL);

    int status = 0;
    for (int i = 0; i < P; i++) {
        if (reed_solomon_rlrs_encode_block(single, data_shards, parity_indices[i], single_parity[i], T) < 0)
            status = -1;
    }
    if (reed_solomon_rlrs_encode_many(bulk, data_shards, parity_indices, bulk_parity, P, T) < 0)
        status = -1;

    int invalid_index = 65536 - 32;
    uint16_t *invalid_output = bulk_parity[0];
    if (reed_solomon_rlrs_encode_many(bulk, data_shards, &invalid_index, &invalid_output, 1, T) == 0)
        status = 1;

    uint16_t *all_shards[K + P];
    uint8_t marks[K + P];
    int invalid_parity_indices[P];
    memset(marks, 0, sizeof(marks));
    for (int i = 0; i < K; i++)
        all_shards[i] = data_shards[i];
    for (int i = 0; i < P; i++) {
        all_shards[K + i] = bulk_parity[i];
        invalid_parity_indices[i] = parity_indices[i];
    }
    invalid_parity_indices[0] = -1;
    if (reed_solomon_rlrs_decode(bulk, all_shards, marks, invalid_parity_indices, P, T) == 0)
        status = 1;

    for (int i = 0; i < P && status == 0; i++) {
        if (memcmp(single_parity[i], bulk_parity[i], T * sizeof(uint16_t)) != 0) {
            printf("Bulk parity mismatch at output %d (index %d)\n", i, parity_indices[i]);
            status = 1;
        }
    }

    /* Exercise explicit invalidation when source buffers change in place. */
    data_shards[3][T - 1] ^= 0x5a5a;
    reed_solomon_rlrs_invalidate_cache(single);
    reed_solomon_rlrs_invalidate_cache(bulk);
    for (int i = 0; i < P; i++) {
        if (reed_solomon_rlrs_encode_block(single, data_shards, parity_indices[i], single_parity[i], T) < 0)
            status = -1;
    }
    if (reed_solomon_rlrs_encode_many(bulk, data_shards, parity_indices, bulk_parity, P, T) < 0)
        status = -1;
    for (int i = 0; i < P && status == 0; i++) {
        if (memcmp(single_parity[i], bulk_parity[i], T * sizeof(uint16_t)) != 0) {
            printf("Bulk parity mismatch after cache invalidation at output %d\n", i);
            status = 1;
        }
    }

    reed_solomon_rlrs_release(single);
    reed_solomon_rlrs_release(bulk);
    for (int i = 0; i < K; i++)
        free(data_shards[i]);
    for (int i = 0; i < P; i++) {
        free(single_parity[i]);
        free(bulk_parity[i]);
    }
    return status;
}

int test_rlrs_run(int K, int P, int T, int erasures_count, int seed)
{
    srand(seed);

    uint16_t **data_shards = calloc(K, sizeof(uint16_t *));
    uint16_t **orig_data = calloc(K, sizeof(uint16_t *));
    uint16_t **parity_shards = calloc(P, sizeof(uint16_t *));
    int *parity_indices = calloc(P, sizeof(int));
    uint8_t *marks = calloc(K + P, sizeof(uint8_t));

    for (int i = 0; i < K; i++) {
        data_shards[i] = malloc(T * sizeof(uint16_t));
        orig_data[i] = malloc(T * sizeof(uint16_t));
        for (int j = 0; j < T; j++) {
            data_shards[i][j] = rand() % 65536;
            orig_data[i][j] = data_shards[i][j];
        }
    }

    for (int i = 0; i < P; i++) {
        parity_shards[i] = malloc(T * sizeof(uint16_t));
        // Pick arbitrary non-contiguous parity indices (e.g. i * 3 + 1)
        parity_indices[i] = i * 3 + 1;
    }

    reed_solomon_rlrs *rs = reed_solomon_rlrs_new(K);
    if (!rs) {
        printf("reed_solomon_rlrs_new failed\n");
        return -1;
    }

    // Encode arbitrary non-contiguous parity blocks in one transform pass.
    int encode_res = reed_solomon_rlrs_encode_many(rs, data_shards, parity_indices, parity_shards, P, T);
    if (encode_res < 0) {
        printf("bulk encode failed\n");
        reed_solomon_rlrs_release(rs);
        return -1;
    }

    // Erase erasures_count data shards
    int erased = 0;
    while (erased < erasures_count) {
        int at = rand() % K;
        if (marks[at] == 0) {
            memset(data_shards[at], 0, T * sizeof(uint16_t));
            marks[at] = 1;
            erased++;
        }
    }

    // Collect all shards: first K are data, next P are parities
    uint16_t **all_shards = calloc(K + P, sizeof(uint16_t *));
    for (int i = 0; i < K; i++) {
        all_shards[i] = data_shards[i];
    }
    for (int i = 0; i < P; i++) {
        all_shards[K + i] = parity_shards[i];
        marks[K + i] = 0; // Parity shards are present
    }

    // Decode
    int decode_res = reed_solomon_rlrs_decode(rs, all_shards, marks, parity_indices, P, T);
    if (decode_res < 0) {
        printf("reed_solomon_rlrs_decode failed\n");
    }

    int status = 0;
    if (decode_res < 0) {
        status = -1;
    } else {
        status = verify_buffers(orig_data, data_shards, K, T) ? 0 : 1;
    }

    /* Decode a changed pattern, then repeat it to exercise both cache paths. */
    for (int pass = 0; pass < 2 && status == 0 && erasures_count > 0; pass++) {
        memset(marks, 0, K + P);
        for (int i = 0; i < erasures_count; i++) {
            int at = (i * 7 + 1) % K;
            marks[at] = 1;
            memset(data_shards[at], 0, T * sizeof(uint16_t));
        }
        decode_res = reed_solomon_rlrs_decode(rs, all_shards, marks, parity_indices, P, T);
        if (decode_res < 0 || !verify_buffers(orig_data, data_shards, K, T))
            status = 1;
    }

    reed_solomon_rlrs_release(rs);

    for (int i = 0; i < K; i++) {
        free(data_shards[i]);
        free(orig_data[i]);
    }
    for (int i = 0; i < P; i++) {
        free(parity_shards[i]);
    }
    free(data_shards);
    free(orig_data);
    free(parity_shards);
    free(parity_indices);
    free(marks);
    free(all_shards);

    return status;
}

int main(void)
{
    printf("=== RUNNING RLRS CODEC TEST BATTERY ===\n");

    reed_solomon_rlrs *oversized = reed_solomon_rlrs_new(32769);
    if (oversized != NULL) {
        reed_solomon_rlrs_release(oversized);
        printf("Oversized data-shard count was accepted\n");
        return 1;
    }

    int status = test_subspace_log_walsh();
    if (status != 0) {
        printf("Subspace locator transform test failed\n");
        return 1;
    }
    printf("Subspace locator transform matches full-field result: OK\n");

    status = test_bulk_matches_single();
    if (status != 0) {
        printf("Bulk encode equivalence test failed\n");
        return 1;
    }
    printf("Bulk encode matches single-block API: OK\n");

    // Test simple K=4, P=2, no erasures
    status = test_rlrs_run(4, 2, 64, 0, 100);
    if (status != 0) {
        printf("Test K=4, P=2, 0 erasures failed\n");
        return 1;
    }
    printf("K=4, P=2, 0 erasures: OK\n");

    // Test K=4, P=2, 1 erasure
    status = test_rlrs_run(4, 2, 64, 1, 101);
    if (status != 0) {
        printf("Test K=4, P=2, 1 erasure failed\n");
        return 1;
    }
    printf("K=4, P=2, 1 erasure: OK\n");

    // Test K=4, P=2, 2 erasures
    status = test_rlrs_run(4, 2, 64, 2, 102);
    if (status != 0) {
        printf("Test K=4, P=2, 2 erasures failed\n");
        return 1;
    }
    printf("K=4, P=2, 2 erasures: OK\n");

    // Test K=50, P=20, 15 erasures
    status = test_rlrs_run(50, 20, 512, 15, 103);
    if (status != 0) {
        printf("Test K=50, P=20, 15 erasures failed\n");
        return 1;
    }
    printf("K=50, P=20, 15 erasures: OK\n");

    printf("=== ALL RLRS TESTS PASSED SUCCESSFULLY ===\n");
    return 0;
}

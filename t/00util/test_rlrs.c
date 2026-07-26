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

    // Encode parity blocks on demand
    for (int i = 0; i < P; i++) {
        int res = reed_solomon_rlrs_encode_block(rs, data_shards, parity_indices[i], parity_shards[i], T);
        if (res < 0) {
            printf("encode block failed for index %d\n", parity_indices[i]);
            reed_solomon_rlrs_release(rs);
            return -1;
        }
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

    reed_solomon_rlrs_release(rs);

    int status = 0;
    if (decode_res < 0) {
        status = -1;
    } else {
        status = verify_buffers(orig_data, data_shards, K, T) ? 0 : 1;
    }

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

    // Test simple K=4, P=2, no erasures
    int status = test_rlrs_run(4, 2, 64, 0, 100);
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

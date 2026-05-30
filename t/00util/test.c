#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "rs.h"

#define MAP(x, max_x, m, n) (m + x / (max_x / (n - m) + 1))

typedef reed_solomon rs_t;

// Helper to check if two buffers are identical
static int verify_buffers(uint8_t **orig, uint8_t **reconstructed, int K, int T)
{
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < T; j++) {
            if (orig[i][j] != reconstructed[i][j]) {
                return 0;
            }
        }
    }
    return 1;
}

// Runs a single RS encode/decode verification
int test_codec_run(int K, int N, int T, int erasures_count, int seed)
{
    uint8_t **buf = calloc(K + N, sizeof(uint8_t *));
    uint8_t **cmp = calloc(K + N, sizeof(uint8_t *));
    uint8_t *marks = calloc(1, K + N);
    int status = 0;

    for (int i = 0; i < K + N; i++) {
        buf[i] = calloc(1, T);
        cmp[i] = calloc(1, T);
    }

    // Populate data shards with deterministic random data
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < T; j++) {
            buf[i][j] = MAP(rand(), RAND_MAX, 0, 256);
            cmp[i][j] = buf[i][j];
        }
    }

    reed_solomon_init();
    rs_t *rs = reed_solomon_new(K, N);
    if (!rs) {
        status = -2; // Initialization failed
        goto cleanup;
    }

    reed_solomon_encode(rs, buf, K + N, T);

    // Erase a specific number of shards
    int erased = 0;
    while (erased < erasures_count) {
        int at = rand() % (K + N);
        if (marks[at] == 0) {
            memset(buf[at], 0, T);
            marks[at] = 1;
            erased++;
        }
    }

    int decode_res = reed_solomon_reconstruct(rs, buf, marks, K + N, T);
    reed_solomon_release(rs);

    if (decode_res < 0) {
        status = -1; // Decode returned error
    } else {
        status = verify_buffers(cmp, buf, K, T) ? 0 : 1; // 0 = SUCCESS, 1 = MISMATCH
    }

cleanup:
    for (int i = 0; i < K + N; i++) {
        free(buf[i]);
        free(cmp[i]);
    }
    free(buf);
    free(cmp);
    free(marks);
    return status;
}

// Runs a single RS encode/decode verification with explicit erasures mask
int test_codec_explicit_erasures(int K, int N, int T, uint32_t erasure_mask)
{
    uint8_t **buf = calloc(K + N, sizeof(uint8_t *));
    uint8_t **cmp = calloc(K + N, sizeof(uint8_t *));
    uint8_t *marks = calloc(1, K + N);
    int status = 0;

    for (int i = 0; i < K + N; i++) {
        buf[i] = calloc(1, T);
        cmp[i] = calloc(1, T);
    }

    // Populate data shards with deterministic random data
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < T; j++) {
            buf[i][j] = MAP(rand(), RAND_MAX, 0, 256);
            cmp[i][j] = buf[i][j];
        }
    }

    reed_solomon_init();
    rs_t *rs = reed_solomon_new(K, N);
    if (!rs) {
        status = -2; // Initialization failed
        goto cleanup;
    }

    reed_solomon_encode(rs, buf, K + N, T);

    for (int i = 0; i < K + N; i++) {
        if ((erasure_mask >> i) & 1) {
            memset(buf[i], 0, T);
            marks[i] = 1;
        }
    }

    int decode_res = reed_solomon_reconstruct(rs, buf, marks, K + N, T);
    reed_solomon_release(rs);

    if (decode_res < 0) {
        status = -1; // Decode returned error
    } else {
        status = verify_buffers(cmp, buf, K, T) ? 0 : 1; // 0 = SUCCESS, 1 = MISMATCH
    }

cleanup:
    for (int i = 0; i < K + N; i++) {
        free(buf[i]);
        free(cmp[i]);
    }
    free(buf);
    free(cmp);
    free(marks);
    return status;
}

// 1. API Parameter Safety Tests
static void run_api_safety_tests(void)
{
    printf("[API SAFETY] Running parameter checks...\n");
    
    // ds + ps > 255
    rs_t *rs1 = reed_solomon_new(250, 10);
    assert(rs1 == NULL);

    // ds <= 0
    rs_t *rs2 = reed_solomon_new(0, 5);
    assert(rs2 == NULL);

    // ps <= 0
    rs_t *rs3 = reed_solomon_new(5, 0);
    assert(rs3 == NULL);

    // Valid instantiation
    rs_t *rs4 = reed_solomon_new(5, 3);
    assert(rs4 != NULL);
    assert(rs4->ds == 5);
    assert(rs4->ps == 3);
    assert(rs4->ts == 8);
    reed_solomon_release(rs4);

    printf("[API SAFETY] OK\n");
}

// 2. Decoder Failure Recovery Tests
static void run_decoder_failure_tests(void)
{
    printf("[DECODER FAILURE] Testing too many erasures...\n");
    // K = 5, N = 3 (Parity = 3, so we can correct at most 3 erasures)
    // Erase 4 shards (should fail)
    int res = test_codec_run(5, 3, 64, 4, 12345);
    assert(res == -1); // Must fail decode
    printf("[DECODER FAILURE] OK\n");
}

// 3. Edge Case Bounds & Alignment Tests
static void run_bounds_and_alignment_tests(void)
{
    printf("[ALIGNMENT & TAILS] Testing non-aligned buffer sizes and SIMD boundary tails...\n");
    
    // Test sizes from 1 to 130 bytes to comprehensively cover all SIMD tail paths (NEON, AVX2, AVX-512)
    for (int T = 1; T <= 130; T++) {
        // Test varying erasure counts
        for (int e = 0; e <= 2; e++) {
            int res = test_codec_run(6, 2, T, e, T * 7 + e);
            if (res != 0) {
                printf("[ALIGNMENT & TAILS] FAILED at T = %d, erasures = %d, res = %d\n", T, e, res);
                exit(1);
            }
        }
    }

    // Extreme configurations
    // K = 1, N = 1 (Minimal dimensions)
    assert(test_codec_run(1, 1, 64, 1, 999) == 0);
    assert(test_codec_run(1, 1, 1, 1, 999) == 0);

    // K = 254, N = 1 (Max data shards)
    assert(test_codec_run(254, 1, 16, 1, 888) == 0);

    // K = 1, N = 254 (Max parity shards)
    assert(test_codec_run(1, 254, 16, 50, 777) == 0);

    printf("[ALIGNMENT & TAILS] OK\n");
}

// 4. Combinatorial Exhaustion for Small Matrices
static void run_combinatorial_tests(void)
{
    printf("[COMBINATORIAL] Testing all combinations for small K, N...\n");
    // Test pairs like (K=3, N=2), (K=4, N=3), (K=5, N=4)
    int test_pairs[][2] = {{3, 2}, {4, 3}, {5, 4}, {6, 2}, {2, 6}, {4, 4}};
    int num_pairs = sizeof(test_pairs) / sizeof(test_pairs[0]);
    int T = 64; // arbitrary size

    for (int p = 0; p < num_pairs; p++) {
        int K = test_pairs[p][0];
        int N = test_pairs[p][1];
        int total = K + N;
        uint32_t max_mask = (1U << total);

        for (uint32_t mask = 0; mask < max_mask; mask++) {
            if (__builtin_popcount(mask) <= N) {
                int res = test_codec_explicit_erasures(K, N, T, mask);
                if (res != 0) {
                    printf("[COMBINATORIAL] FAILED at K=%d, N=%d, mask=0x%x\n", K, N, mask);
                    exit(1);
                }
            }
        }
    }
    printf("[COMBINATORIAL] OK\n");
}

int main(int argc, char *argv[])
{
    int seed = time(NULL);
    if (argc == 2) {
        seed = strtol(argv[1], NULL, 10);
    }
    srand(seed);

    printf("=== RUNNING ROBUST CODEC TEST BATTERY ===\n");
    printf("Seed: %d\n", seed);

    // Execute structured unit tests
    run_api_safety_tests();
    run_decoder_failure_tests();
    run_bounds_and_alignment_tests();
    run_combinatorial_tests();

    // Execute random tests for broader coverage
    printf("[RANDOM SEARCH] Running randomized codec iterations...\n");
    for (int i = 0; i < 200; i++) {
        int K = MAP(rand(), RAND_MAX, 1, 50);
        int N = MAP(rand(), RAND_MAX, 1, 20);
        if ((K + N) > 255) {
            K = 255 - N;
        }
        int T = MAP(rand(), RAND_MAX, 1, 1024);
        int e = rand() % (N + 1); // Random erasures <= N

        int res = test_codec_run(K, N, T, e, seed + i);
        if (res != 0) {
            printf("[RANDOM SEARCH] Mismatch or failure at K = %d, N = %d, T = %d, e = %d\n", K, N, T, e);
            exit(1);
        }
    }
    printf("[RANDOM SEARCH] OK\n");

    printf("===ALL TESTS PASSED SUCCESSFULLY=== \n===OK===\n");
    return 0;
}

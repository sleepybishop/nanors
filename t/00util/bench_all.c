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

#include "rs.h"
#include "rs16_afft.h"
#include "rs_rlrs.h"

double now(time_t epoch)
{
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return ((now.tv_sec - epoch) + now.tv_nsec / 1000000000.0);
}

void bench_rs(int K, int P, int T, int iterations)
{
    if (K + P > 256) {
        printf("%-10s | %4d | %4d | %8s | %8s\n", "RS (8-bit)", K, P, "N/A", "N/A");
        return;
    }

    uint8_t **buf = calloc(K + P, sizeof(uint8_t *));
    uint8_t *marks = calloc(K + P, sizeof(uint8_t));
    for (int i = 0; i < K + P; i++) {
        buf[i] = aligned_alloc(64, T * sizeof(uint16_t));
        memset(buf[i], rand() % 256, T * sizeof(uint16_t));
    }

    reed_solomon *rs = reed_solomon_new(K, P);
    if (!rs) {
        printf("Failed to init RS\n");
        return;
    }

    double t_enc = 0.0;
    for (int it = 0; it < iterations; it++) {
        double t0 = now(0);
        reed_solomon_encode(rs, buf, K + P, T * sizeof(uint16_t));
        t_enc += now(0) - t0;
    }

    // Set erasures
    for (int i = 0; i < P; i++) {
        marks[i] = 1;
        memset(buf[i], 0, T * sizeof(uint16_t));
    }

    double t_dec = 0.0;
    for (int it = 0; it < iterations; it++) {
        // Reset buffers for decoding
        for (int i = 0; i < P; i++) {
            memset(buf[i], 0, T * sizeof(uint16_t));
            marks[i] = 1;
        }
        double t0 = now(0);
        reed_solomon_decode(rs, buf, marks, K + P, T * sizeof(uint16_t));
        t_dec += now(0) - t0;
    }

    reed_solomon_release(rs);

    for (int i = 0; i < K + P; i++) {
        free(buf[i]);
    }
    free(buf);
    free(marks);

    double mb_processed = (double)K * T * sizeof(uint16_t) * iterations / (1024.0 * 1024.0);
    printf("%-10s | %4d | %4d | %6.1f MB/s | %6.1f MB/s\n", "RS (8-bit)", K, P, mb_processed / t_enc, mb_processed / t_dec);
}

void bench_rs16_afft(int K, int P, int T, int iterations)
{
    if (K + P > 65536) {
        printf("%-10s | %4d | %4d | %8s | %8s\n", "RS16_AFFT", K, P, "N/A", "N/A");
        return;
    }

    uint16_t **buf = calloc(K + P, sizeof(uint16_t *));
    uint8_t *marks = calloc(K + P, sizeof(uint8_t));
    for (int i = 0; i < K + P; i++) {
        buf[i] = aligned_alloc(64, T * sizeof(uint16_t));
        for (int j = 0; j < T; j++) {
            buf[i][j] = rand() % 65536;
        }
    }

    reed_solomon16_afft *rs = reed_solomon16_afft_new(K, P);
    if (!rs) {
        printf("Failed to init RS16_AFFT\n");
        return;
    }

    double t_enc = 0.0;
    for (int it = 0; it < iterations; it++) {
        double t0 = now(0);
        reed_solomon16_afft_encode(rs, buf, K + P, T);
        t_enc += now(0) - t0;
    }

    for (int i = 0; i < P; i++) {
        marks[i] = 1;
        memset(buf[i], 0, T * sizeof(uint16_t));
    }

    double t_dec = 0.0;
    for (int it = 0; it < iterations; it++) {
        for (int i = 0; i < P; i++) {
            memset(buf[i], 0, T * sizeof(uint16_t));
            marks[i] = 1;
        }
        double t0 = now(0);
        reed_solomon16_afft_decode(rs, buf, marks, K + P, T);
        t_dec += now(0) - t0;
    }

    reed_solomon16_afft_release(rs);

    for (int i = 0; i < K + P; i++) {
        free(buf[i]);
    }
    free(buf);
    free(marks);

    double mb_processed = (double)K * T * sizeof(uint16_t) * iterations / (1024.0 * 1024.0);
    printf("%-10s | %4d | %4d | %6.1f MB/s | %6.1f MB/s\n", "RS16_AFFT", K, P, mb_processed / t_enc, mb_processed / t_dec);
}

void bench_rlrs(int K, int P, int T, int iterations)
{
    uint16_t **data_shards = calloc(K, sizeof(uint16_t *));
    uint16_t **parity_shards = calloc(P, sizeof(uint16_t *));
    int *parity_indices = calloc(P, sizeof(int));
    uint8_t *marks = calloc(K + P, sizeof(uint8_t));

    for (int i = 0; i < K; i++) {
        data_shards[i] = aligned_alloc(64, T * sizeof(uint16_t));
        for (int j = 0; j < T; j++) {
            data_shards[i][j] = rand() % 65536;
        }
    }
    for (int i = 0; i < P; i++) {
        parity_shards[i] = aligned_alloc(64, T * sizeof(uint16_t));
        parity_indices[i] = i;
    }

    reed_solomon_rlrs *rs = reed_solomon_rlrs_new(K);
    if (!rs) {
        printf("Failed to init RLRS\n");
        return;
    }

    double t_enc = 0.0;
    for (int it = 0; it < iterations; it++) {
        double t0 = now(0);
        for (int i = 0; i < P; i++) {
            reed_solomon_rlrs_encode_block(rs, data_shards, parity_indices[i], parity_shards[i], T);
        }
        t_enc += now(0) - t0;
    }

    // Collect all shards
    uint16_t **all_shards = calloc(K + P, sizeof(uint16_t *));
    for (int i = 0; i < K; i++) {
        all_shards[i] = data_shards[i];
    }
    for (int i = 0; i < P; i++) {
        all_shards[K + i] = parity_shards[i];
    }

    double t_dec = 0.0;
    for (int it = 0; it < iterations; it++) {
        // Erase first P data shards
        memset(marks, 0, K + P);
        for (int i = 0; i < P; i++) {
            marks[i] = 1;
            memset(data_shards[i], 0, T * sizeof(uint16_t));
        }

        double t0 = now(0);
        reed_solomon_rlrs_decode(rs, all_shards, marks, parity_indices, P, T);
        t_dec += now(0) - t0;
    }

    reed_solomon_rlrs_release(rs);

    for (int i = 0; i < K; i++) {
        free(data_shards[i]);
    }
    for (int i = 0; i < P; i++) {
        free(parity_shards[i]);
    }
    free(data_shards);
    free(parity_shards);
    free(parity_indices);
    free(marks);
    free(all_shards);

    double mb_processed = (double)K * T * sizeof(uint16_t) * iterations / (1024.0 * 1024.0);
    printf("%-10s | %4d | %4d | %6.1f MB/s | %6.1f MB/s\n", "RLRS (16b)", K, P, mb_processed / t_enc, mb_processed / t_dec);
}

int main(void)
{
    printf("=========================================================\n");
    printf("%-10s |    K |    P |   Enc Speed  |   Dec Speed  \n", "Codec");
    printf("=========================================================\n");

    int sizes[] = {5, 10, 50, 100, 500, 1000};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_sizes; i++) {
        int K = sizes[i];
        int P = K / 5; // 20% parity
        if (P < 2)
            P = 2; // Minimum 2 parity shards

        int T = 1024; // 1024 words = 2KB per shard
        int iterations = 1000;
        if (K >= 500) {
            iterations = 100; // Reduce iterations for larger sizes to avoid long runtimes
        }

        bench_rs(K, P, T, iterations);
        bench_rs16_afft(K, P, T, iterations);
        bench_rlrs(K, P, T, iterations);
        printf("---------------------------------------------------------\n");
    }

    printf("\n=========================================================\n");
    printf(" LARGE BATCH BENCHMARK (T=1260 payload, Batched 2.5MB/shard)\n");
    printf("=========================================================\n");
    printf("%-10s |    K |    P |   Enc Speed  |   Dec Speed  \n", "Codec");
    printf("=========================================================\n");

    for (int i = 0; i < num_sizes; i++) {
        int K = sizes[i];
        int P = K / 5;
        if (P < 2)
            P = 2;

        int T = 1260 * 1000 / K; // Scaling T so total batch size per call is ~2.5 MB total
        if (T < 1260)
            T = 1260;
        int iterations = 100;

        bench_rs(K, P, T, iterations);
        bench_rs16_afft(K, P, T, iterations);
        bench_rlrs(K, P, T, iterations);
        printf("---------------------------------------------------------\n");
    }

    return 0;
}

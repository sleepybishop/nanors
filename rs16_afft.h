#ifndef RS16_AFFT_H
#define RS16_AFFT_H

#include "rs16.h"
#include "deps/obl/oblas16.h"
#include "deps/obl/oblas16_afft.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BATCH_SIZE 512

struct afft_workspace {
    uint16_t *chunk_buf;
    uint8_t *needed_buf;
    uint16_t *acc;
    uint16_t *gamma_arr;
    uint16_t *L_eval;
    uint16_t *inv_Lp_eval;
    uint8_t *is_erased;
    uint16_t *error_locations;
    uint16_t *buf;
    int *erasures;
    int *E;
    void *flat_alloc;
    int gamma_count;
    volatile int op_lock;
    void (*axpy)(uint16_t *a, const uint16_t *b, uint16_t u, unsigned k);
    void (*axiy)(uint16_t *dst, const uint16_t *src, uint16_t twist, unsigned batch);
    struct oblas16_impl o16;
    struct oblas16_afft_impl afft;
};

extern uint16_t LogWalsh[65536];
void fwht_mod(uint16_t *restrict data, int n);

typedef struct {
    int ds;
    int ps;
    int ts;
    void *p;
} reed_solomon16_afft;

void reed_solomon16_afft_init(void);
/* next_power_of_two(data_shards) + parity_shards must not exceed 65536. */
reed_solomon16_afft *reed_solomon16_afft_new(int data_shards, int parity_shards);
void reed_solomon16_afft_release(reed_solomon16_afft *rs);

/* Encode and decode operations on one codec instance are serialized internally. */
int reed_solomon16_afft_encode(reed_solomon16_afft *rs, uint16_t **shards, int nr_shards, int bs);
/* Reconstructs every data or parity shard whose corresponding marks entry is nonzero. */
int reed_solomon16_afft_decode(reed_solomon16_afft *rs, uint16_t **shards, uint8_t *marks, int nr_shards, int bs);

#ifdef __cplusplus
}
#endif

#endif /* RS16_AFFT_H */

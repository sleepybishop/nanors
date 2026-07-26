#ifndef RS_RLRS_H
#define RS_RLRS_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct reed_solomon_rlrs reed_solomon_rlrs;

/*
 * Initialize a rateless codec instance.
 * @param data_shards Number of source data shards (K).
 * @return A pointer to the initialized codec instance, or NULL on error.
 */
reed_solomon_rlrs *reed_solomon_rlrs_new(int data_shards);

/*
 * Free the codec instance and its resources.
 * @param rs Pointer to the codec instance.
 */
void reed_solomon_rlrs_release(reed_solomon_rlrs *rs);

/*
 * Invalidate the cached IFFT coefficients.
 * Call this if the contents of the data shards change but the pointers stay the same.
 */
void reed_solomon_rlrs_invalidate_cache(reed_solomon_rlrs *rs);

/*
 * Encode a single parity block by index.
 * @param rs Pointer to the codec instance.
 * @param data_shards Array of K pointers to source data buffers.
 * @param parity_index The index of the parity block (0 to 65535 - K).
 * @param out_parity Output buffer for the generated parity block.
 * @param block_size Size of the buffers in bytes.
 * @return 0 on success, or -1 on error.
 */
int reed_solomon_rlrs_encode_block(reed_solomon_rlrs *rs, uint16_t **data_shards, int parity_index, uint16_t *out_parity,
                                   int block_size);

/*
 * Decode the erased data/parity shards.
 * @param rs Pointer to the codec instance.
 * @param shards Array containing data & parity buffers. The first K entries are the source shards.
 * @param marks Boolean array of size K + received_parity_count. 1 if shard is missing/erased, 0 otherwise.
 * @param received_parity_indices Array containing the index of each received parity block (relative to 0).
 * @param received_parity_count Number of parity shards provided. Must be >= the number of erasures.
 * @param block_size Size of the buffers in bytes.
 * @return 0 on success, or -1 on error.
 */
int reed_solomon_rlrs_decode(reed_solomon_rlrs *rs, uint16_t **shards, uint8_t *marks, int *received_parity_indices,
                             int received_parity_count, int block_size);

#ifdef __cplusplus
}
#endif

#endif /* RS_RLRS_H */

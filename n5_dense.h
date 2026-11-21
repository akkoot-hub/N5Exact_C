#ifndef N5_DENSE_H
#define N5_DENSE_H

#include <stdint.h>

#include "n5_sparse.h"
#include "n5_slicing.h"

#define N5_DENSE_MAX_MODULI 8

/* Result of a shared-structure multi-modulus dense contraction. */
typedef struct
{
    uint8_t modulus_count;
    uint32_t values[N5_DENSE_MAX_MODULI];
    uint64_t max_table_entries;
    unsigned int eliminated_variables;
} N5DenseMultiModularContractionResult;

/*
 * Contract one sparse network for several independent CRT moduli.
 * Structural planning and dense-table traversal are shared; each residue
 * lane remains independent in its own Z/pZ arithmetic.
 */
int n5_dense_modular_contract_sparse_network_multi(
    const N5SparseNetwork* network,
    const N5SparseContractionPlan* plan,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5DenseMultiModularContractionResult* result
);

/* Result of a sliced, shared-structure multi-modulus dense contraction. */
typedef struct
{
    uint8_t modulus_count;
    uint32_t values[N5_DENSE_MAX_MODULI];
    uint64_t raw_assignment_count;
    uint64_t zero_assignment_count;
    uint64_t max_table_entries;
} N5DenseSlicedMultiResult;

/*
 * Enumerate a slice plan and contract all modulus lanes over one shared
 * conditioning and dense-structure traversal.
 */
int n5_dense_modular_contract_sliced_multi(
    const N5SparseNetwork* network,
    const N5SlicePlan* slice_plan,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5DenseSlicedMultiResult* result
);

#endif

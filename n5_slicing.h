#ifndef N5_SLICING_H
#define N5_SLICING_H

#include <stdint.h>

#include "n5_sparse.h"

#define N5_SLICING_MAX_VARIABLES 8

/* A structural slice plan over concrete sparse-variable domains. */
typedef struct
{
    uint8_t variable_count;
    uint8_t variables[N5_SLICING_MAX_VARIABLES];
    uint64_t estimated_peak;
    uint64_t raw_assignment_count;
    uint64_t target_peak;
    int target_reached;
} N5SlicePlan;

/*
 * Compiled conditioning layout for a fixed slice-variable scope.
 * Fixed values change between representatives; variable positions do not.
 */
typedef struct
{
    const N5SparseNetwork* source;
    unsigned int fixed_count;
    uint8_t fixed_variables[N5_SLICING_MAX_VARIABLES];
    struct
    {
        uint8_t fixed_count;
        uint8_t fixed_positions[N5_SPARSE_MAX_ARITY];
        uint8_t fixed_value_indices[N5_SPARSE_MAX_ARITY];
        uint8_t keep_count;
        uint8_t keep_positions[N5_SPARSE_MAX_ARITY];
        uint8_t output_vars[N5_SPARSE_MAX_ARITY];
        uint8_t canonical_input;
    } factor[N5_SPARSE_NETWORK_MAX_FACTORS];
} N5PreparedCondition;

/* Estimate the dense structural peak after fixing the selected variables. */
int n5_slicing_estimate_peak(
    const N5SparseNetwork* network,
    const uint8_t* variables,
    unsigned int variable_count,
    uint64_t* estimated_peak
);

/* Greedily choose slice variables until target_peak or max_slices is reached. */
int n5_slicing_choose_greedy(
    const N5SparseNetwork* network,
    uint64_t target_peak,
    unsigned int max_slices,
    N5SlicePlan* plan
);

/* Compile a reusable conditioning layout for a fixed variable scope. */
int n5_slicing_prepare_condition(
    const N5SparseNetwork* network,
    const uint8_t* fixed_variables,
    unsigned int fixed_count,
    N5PreparedCondition* compiled
);

/* Restrict one factor using the compiled layout; caller owns its storage. */
int n5_slicing_condition_prepared_factor(
    const N5PreparedCondition* compiled,
    unsigned int factor_index,
    const uint8_t* fixed_values,
    N5SparseFactor* conditioned
);

/* Condition a complete network using a previously compiled layout. */
int n5_slicing_condition_prepared(
    const N5PreparedCondition* compiled,
    const uint8_t* fixed_values,
    N5SparseNetwork* output
);

/* Fix selected variables and create an independently owned output network. */
int n5_slicing_condition_sparse_network(
    const N5SparseNetwork* network,
    const uint8_t* fixed_variables,
    const uint8_t* fixed_values,
    unsigned int fixed_count,
    N5SparseNetwork* output
);

#endif

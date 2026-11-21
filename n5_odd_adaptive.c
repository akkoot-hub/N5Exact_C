#include "n5_odd_adaptive.h"

#include <stdlib.h>
#include <string.h>

#include "n5_contraction.h"
#include "n5_dense.h"
#include "n5_kernel.h"
#include "n5_odd_network.h"
#ifdef N5_STEP11_PROFILE_TIMING
#include <stdio.h>
#include "n5_parallel.h"
#endif

/*
 * 一个 fixed-r k=2 second-Walsh decomposition。
 * channel IDs 使用 canonical k=2 even basis 的原始 ID。
 */
typedef struct
{
    N5V3 r;

    uint8_t channel_count;
    uint8_t channel_ids[7];
    int8_t gauge[7];

    uint8_t lambda_count;
    N5V3 lambda_values[8];

    int scale;
} N5OddK2WalshDecomposition;

/* fixed-K cardinality helper。 */
static unsigned int n5_odd_adaptive_cell_size(
    const N5SizeMatrix* size_matrix,
    unsigned int cell
)
{
    return size_matrix->matrix[cell >> 2u][cell & 3u];
}

/* 在 decomposition 中找 canonical q ID 的 local position。 */
static int n5_odd_adaptive_q_position(
    const N5OddK2WalshDecomposition* decomposition,
    unsigned int q
)
{
    unsigned int index;

    for (index = 0u; index < decomposition->channel_count; ++index)
    {
        if (decomposition->channel_ids[index] == q)
        {
            return (int)index;
        }
    }

    return -1;
}

/* definition-derived Walsh leg sign。 */
static int n5_odd_adaptive_leg_value(
    const N5OddK2WalshDecomposition* decomposition,
    const N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT],
    unsigned int q,
    N5V3 lambda
)
{
    int position;

    position = n5_odd_adaptive_q_position(decomposition, q);
    if (position < 0)
    {
        return 0;
    }

    return
        (int)decomposition->gauge[position] *
        n5_walsh_character(lambda, channels[q].s);
}

/*
 * 从当前 canonical all-even k=2 basis 的 transformed cell tensor
 * 自动恢复 fixed-r gauge；不保存/读取历史 sign table。
 */
static int n5_odd_adaptive_build_k2_decomposition(
    N5V3 r,
    N5OddK2WalshDecomposition* decomposition
)
{
    N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT];
    unsigned int channel_count;
    unsigned int q;
    unsigned int candidate_mask;
    unsigned int candidate_count;

#ifndef N5_STEP8_DISABLE_K2_GAUGE_CACHE
    /* The definition-derived k=2 decomposition depends on r alone.
       Each worker owns eight immutable, lazily computed values; no locks,
       input-dependent shortcut, or cross-modulus arithmetic is involved. */
#if defined(_MSC_VER)
    static __declspec(thread) uint8_t cached[8];
    static __declspec(thread) N5OddK2WalshDecomposition saved[8];
#else
    static _Thread_local uint8_t cached[8];
    static _Thread_local N5OddK2WalshDecomposition saved[8];
#endif
#endif
    if (decomposition == NULL || r >= N5_V3_SIZE)
    {
        return 0;
    }
#ifndef N5_STEP8_DISABLE_K2_GAUGE_CACHE
    if (cached[(unsigned int)r] != 0u)
    {
        *decomposition = saved[(unsigned int)r];
        return 1;
    }
#endif

    memset(decomposition, 0, sizeof(*decomposition));

    channel_count = n5_generate_even_fourier_channels(2u, channels, 28u);
    if (channel_count != 28u)
    {
        return 0;
    }

    decomposition->r = r;

    for (q = 0u; q < channel_count; ++q)
    {
        if (channels[q].r == r)
        {
            if (decomposition->channel_count >= 7u)
            {
                return 0;
            }

            decomposition->channel_ids[decomposition->channel_count] =
                (uint8_t)q;
            ++decomposition->channel_count;
        }
    }

    if (decomposition->channel_count == 0u)
    {
        return 0;
    }

    candidate_count = 1u << (decomposition->channel_count - 1u);

    for (candidate_mask = 0u;
         candidate_mask < candidate_count;
         ++candidate_mask)
    {
        unsigned int index;
        unsigned int ia;
        unsigned int ib;
        unsigned int ic;
        unsigned int nonzero_checked;
        int ok;

        decomposition->gauge[0] = 1;
        for (index = 1u; index < decomposition->channel_count; ++index)
        {
            decomposition->gauge[index] =
                ((candidate_mask >> (index - 1u)) & 1u) ? -1 : 1;
        }

        ok = 1;
        nonzero_checked = 0u;

        for (ia = 0u; ia < decomposition->channel_count && ok; ++ia)
        {
            for (ib = 0u; ib < decomposition->channel_count && ok; ++ib)
            {
                for (ic = 0u; ic < decomposition->channel_count; ++ic)
                {
                    unsigned int qa;
                    unsigned int qb;
                    unsigned int qc;
                    int actual;
                    int predicted_sign;

                    qa = decomposition->channel_ids[ia];
                    qb = decomposition->channel_ids[ib];
                    qc = decomposition->channel_ids[ic];

                    actual = n5_even_cell_tensor_value(
                        &channels[qa],
                        &channels[qb],
                        &channels[qc]
                    );

                    if (actual == 0)
                    {
                        continue;
                    }

                    ++nonzero_checked;

                    if (actual != 32 && actual != -32)
                    {
                        return 0;
                    }

                    predicted_sign =
                        (int)decomposition->gauge[ia] *
                        (int)decomposition->gauge[ib] *
                        (int)decomposition->gauge[ic];

                    if ((actual > 0 ? 1 : -1) != predicted_sign)
                    {
                        ok = 0;
                        break;
                    }
                }
            }
        }

        if (ok && nonzero_checked > 0u)
        {
            break;
        }
    }

    if (candidate_mask == candidate_count)
    {
        return 0;
    }

    if (r == 0u)
    {
        unsigned int lambda;

        for (lambda = 0u; lambda < 8u; ++lambda)
        {
            decomposition->lambda_values[decomposition->lambda_count] =
                (N5V3)lambda;
            ++decomposition->lambda_count;
        }

        decomposition->scale = 4;
    }
    else
    {
        uint8_t seen[8];
        unsigned int lambda;

        memset(seen, 0, sizeof(seen));

        for (lambda = 0u; lambda < 8u; ++lambda)
        {
            unsigned int partner;
            unsigned int representative;

            if (seen[lambda] != 0u)
            {
                continue;
            }

            partner = lambda ^ r;
            representative = lambda < partner ? lambda : partner;

            decomposition->lambda_values[decomposition->lambda_count] =
                (N5V3)representative;
            ++decomposition->lambda_count;

            seen[lambda] = 1u;
            seen[partner] = 1u;
        }

        if (decomposition->lambda_count != 4u)
        {
            return 0;
        }

        decomposition->scale = 8;
    }

#ifndef N5_STEP8_DISABLE_K2_GAUGE_CACHE
    saved[(unsigned int)r] = *decomposition;
    cached[(unsigned int)r] = 1u;
#endif
    return 1;
}

/* 构造 q-incidence -> shared lambda_cell 的 exact leg/norm factor。 */
static int n5_odd_adaptive_build_incidence_leg(
    unsigned int cell,
    unsigned int side,
    const N5OddK2WalshDecomposition* decomposition,
    N5SparseFactor* factor
)
{
    N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT];
    unsigned int channel_count;
    unsigned int q_index;
    unsigned int lambda_index;
    uint8_t vars[2];

    if (decomposition == NULL || factor == NULL || cell >= 16u || side >= 3u)
    {
        return 0;
    }

    channel_count = n5_generate_even_fourier_channels(2u, channels, 28u);
    if (channel_count != 28u)
    {
        return 0;
    }

    vars[0] = (uint8_t)n5_contraction_incidence_var(cell, side);
    vars[1] = (uint8_t)(48u + cell);

    if (!n5_sparse_factor_init(
            factor,
            vars,
            2u,
            (unsigned int)decomposition->channel_count *
            (unsigned int)decomposition->lambda_count
        ))
    {
        return 0;
    }

    for (q_index = 0u;
         q_index < decomposition->channel_count;
         ++q_index)
    {
        unsigned int q;

        q = decomposition->channel_ids[q_index];

        for (lambda_index = 0u;
             lambda_index < decomposition->lambda_count;
             ++lambda_index)
        {
            N5V3 lambda;
            int leg;
            uint8_t assignment[2];

            lambda = decomposition->lambda_values[lambda_index];
            leg = n5_odd_adaptive_leg_value(
                decomposition,
                channels,
                q,
                lambda
            );

            assignment[0] = (uint8_t)q;
            assignment[1] = lambda;

            if (!n5_sparse_factor_append(
                    factor,
                    assignment,
                    n5_small_rational_make(leg, channels[q].norm)
                ))
            {
                n5_sparse_factor_free(factor);
                return 0;
            }
        }
    }

    return n5_sparse_factor_compress(factor);
}

/* deep clone 一个 factor，并保持原 vars。 */
static int n5_odd_adaptive_clone_factor(
    const N5SparseFactor* source,
    N5SparseFactor* output
)
{
    return n5_sparse_factor_clone_with_vars(source, source->vars, output);
}

/* 判断 network 是否已经因为某个 factor 无非零 entries 而为零。 */
static int n5_odd_adaptive_network_is_zero(const N5SparseNetwork* network)
{
    unsigned int factor_index;

    if (network == NULL)
    {
        return 1;
    }

    for (factor_index = 0u;
         factor_index < network->factor_count;
         ++factor_index)
    {
        if (network->factors[factor_index].nnz == 0u)
        {
            return 1;
        }
    }

    return 0;
}

int n5_odd_adaptive_build_second_walsh_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[16],
    const N5OddKernelLibrary* library,
    const N5SparseNetwork* ordinary_network,
    N5SparseNetwork* output
)
{
    unsigned int cell;
    unsigned int line_id;
    int64_t total_scale;

    (void)library;

    if (size_matrix == NULL ||
        r_by_cell == NULL ||
        ordinary_network == NULL ||
        output == NULL ||
        ordinary_network->factor_count != 28u)
    {
        return 0;
    }

    n5_sparse_network_init(output);
    total_scale = 1;

    /* ordinary factor[0..15] 是 normalized cell factors。 */
    for (cell = 0u; cell < 16u; ++cell)
    {
        unsigned int k;

        k = n5_odd_adaptive_cell_size(size_matrix, cell);

        if (k == 2u)
        {
            N5OddK2WalshDecomposition decomposition;

            if (!n5_odd_adaptive_build_k2_decomposition(
                    r_by_cell[cell],
                    &decomposition
                ))
            {
                n5_sparse_network_free(output);
                return 0;
            }

            if (total_scale > INT64_MAX / decomposition.scale)
            {
                n5_sparse_network_free(output);
                return 0;
            }

            total_scale *= decomposition.scale;
            continue;
        }

        {
            N5SparseFactor factor;

            if (!n5_odd_adaptive_clone_factor(
                    &ordinary_network->factors[cell],
                    &factor
                ) ||
                !n5_sparse_network_add_move(output, &factor))
            {
                n5_sparse_factor_free(&factor);
                n5_sparse_network_free(output);
                return 0;
            }
        }
    }

    /* ordinary factor[16..27] 是 ordered partition line factors。 */
    for (line_id = 0u; line_id < 12u; ++line_id)
    {
        unsigned int cells[4];
        unsigned int side;
        unsigned int position;
        N5SparseFactor current;

        if (!n5_contraction_line_cells(line_id, cells))
        {
            n5_sparse_network_free(output);
            return 0;
        }

        side = n5_contraction_line_side(line_id);
        if (side > 2u)
        {
            n5_sparse_network_free(output);
            return 0;
        }

        if (!n5_odd_adaptive_clone_factor(
                &ordinary_network->factors[16u + line_id],
                &current
            ))
        {
            n5_sparse_network_free(output);
            return 0;
        }

        for (position = 0u; position < 4u; ++position)
        {
            unsigned int incident_cell;

            incident_cell = cells[position];

            if (n5_odd_adaptive_cell_size(size_matrix, incident_cell) == 2u)
            {
                N5OddK2WalshDecomposition decomposition;
                N5SparseFactor leg;
                N5SparseFactor multiplied;
                N5SparseFactor summed;
                unsigned int q_var;

                if (!n5_odd_adaptive_build_k2_decomposition(
                        r_by_cell[incident_cell],
                        &decomposition
                    ) ||
                    !n5_odd_adaptive_build_incidence_leg(
                        incident_cell,
                        side,
                        &decomposition,
                        &leg
                    ))
                {
                    n5_sparse_factor_free(&current);
                    n5_sparse_network_free(output);
                    return 0;
                }

                if (!n5_sparse_factor_multiply(&current, &leg, &multiplied))
                {
                    n5_sparse_factor_free(&leg);
                    n5_sparse_factor_free(&current);
                    n5_sparse_network_free(output);
                    return 0;
                }

                q_var = n5_contraction_incidence_var(incident_cell, side);

                if (!n5_sparse_factor_sum_out(&multiplied, q_var, &summed))
                {
                    n5_sparse_factor_free(&multiplied);
                    n5_sparse_factor_free(&leg);
                    n5_sparse_factor_free(&current);
                    n5_sparse_network_free(output);
                    return 0;
                }

                n5_sparse_factor_free(&multiplied);
                n5_sparse_factor_free(&leg);
                n5_sparse_factor_free(&current);
                current = summed;
            }
        }

        if (!n5_sparse_network_add_move(output, &current))
        {
            n5_sparse_factor_free(&current);
            n5_sparse_network_free(output);
            return 0;
        }
    }

    if (total_scale != 1)
    {
        N5SparseFactor scalar;

        if (!n5_sparse_factor_make_scalar(
                n5_small_rational_make(total_scale, 1),
                &scalar
            ) ||
            !n5_sparse_network_add_move(output, &scalar))
        {
            n5_sparse_factor_free(&scalar);
            n5_sparse_network_free(output);
            return 0;
        }
    }

    return 1;
}

/* 收集一个 variable 在 network 中真实出现过的 distinct values。 */
static int n5_odd_adaptive_collect_variable_domain(
    const N5SparseNetwork* network,
    unsigned int variable,
    uint8_t values[N5_SPARSE_MAX_VALUE + 1],
    unsigned int* value_count
)
{
    uint8_t seen[N5_SPARSE_MAX_VALUE + 1];
    unsigned int factor_index;
    unsigned int count;
    unsigned int value;

    if (network == NULL || values == NULL || value_count == NULL ||
        variable >= N5_SPARSE_MAX_VARIABLES)
    {
        return 0;
    }

    memset(seen, 0, sizeof(seen));

    for (factor_index = 0u;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        int position;
        unsigned int p;

        factor = &network->factors[factor_index];
        position = -1;

        for (p = 0u; p < factor->arity; ++p)
        {
            if (factor->vars[p] == variable)
            {
                position = (int)p;
                break;
            }
        }

        if (position >= 0)
        {
            unsigned int entry_index;

            for (entry_index = 0u; entry_index < factor->nnz; ++entry_index)
            {
                seen[factor->entries[entry_index].assignment[position]] = 1u;
            }
        }
    }

    count = 0u;
    for (value = 0u; value <= N5_SPARSE_MAX_VALUE; ++value)
    {
        if (seen[value] != 0u)
        {
            values[count] = (uint8_t)value;
            ++count;
        }
    }

    *value_count = count;
    return count > 0u;
}

#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
/*
 * Cross-slice exact factor-restriction memo, NOT a cache of entire networks.
 * A source factor depends only on the selected slice variables in its scope.
 * All other slice coordinates can therefore share the same exact restriction.
 * The source and cached entries remain immutable until the class finishes.
 * The dense six-prime contraction itself is never approximated or memoized.
 */
#define N5_STEP11_FACTOR_CACHE_MAX_VARIANTS 128u
#define N5_STEP11_FACTOR_CACHE_BUDGET_BYTES UINT64_C(4194304)

typedef struct
{
    N5SparseFactor* variants;
    unsigned int variant_count;
} N5Step11FactorCacheEntry;
#endif

/*
 * 递归枚举 slice assignments；每个 nonzero slice 只构造一次 conditioned network，
 * 然后 shared dense traversal 同时计算全部六个 primes。
 */
typedef struct
{
    const N5SparseNetwork* network;
    const N5SlicePlan* slice_plan;

    uint8_t domain_values[N5_SLICING_MAX_VARIABLES]
                         [N5_SPARSE_MAX_VALUE + 1];
    unsigned int domain_count[N5_SLICING_MAX_VARIABLES];

    uint8_t current_values[N5_SLICING_MAX_VARIABLES];
#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
    unsigned int current_indices[N5_SLICING_MAX_VARIABLES];
    N5Step11FactorCacheEntry factor_cache[N5_SPARSE_NETWORK_MAX_FACTORS];
#endif

#ifndef N5_STEP11_DISABLE_COMPILED_CONDITION
    /* Slice scope remains fixed for the entire charge class; values change. */
    N5PreparedCondition prepared_condition;
#endif

    uint32_t residues[N5_CRT_PRIME_COUNT];
    uint64_t visited_assignments;
    uint64_t zero_assignments;
    uint64_t max_dense_table_entries;
#ifdef N5_STEP11_PROFILE_TIMING
    double condition_seconds;
    double plan_seconds;
    double dense_seconds;
    double precompute_seconds;
#endif
    int success;
} N5OddAdaptiveSliceContext;

#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
static void n5_step11_cache_release(N5OddAdaptiveSliceContext* ctx)
{
    unsigned int f, j;
    for (f = 0; f < ctx->network->factor_count; ++f) {
        N5Step11FactorCacheEntry* c = &ctx->factor_cache[f];
        for (j = 0; j < c->variant_count; ++j)
            n5_sparse_factor_free(&c->variants[j]);
        free(c->variants);
        c->variants = NULL;
        c->variant_count = 0;
    }
}

static void n5_step11_cache_prepare(N5OddAdaptiveSliceContext* ctx)
{
    unsigned int f, k;
    uint64_t remaining = N5_STEP11_FACTOR_CACHE_BUDGET_BYTES;
    for (f = 0; f < ctx->network->factor_count; ++f) {
        const N5SparseFactor* src = &ctx->network->factors[f];
        unsigned int fixed = ctx->prepared_condition.factor[f].fixed_count;
        unsigned int combinations = 1u;
        uint64_t budget;
        N5Step11FactorCacheEntry* cache = &ctx->factor_cache[f];
        if (fixed == 0u) continue; /* source factor is directly reusable */
        for (k = 0; k < fixed; ++k) {
            unsigned int v = ctx->prepared_condition.factor[f].fixed_value_indices[k];
            if (ctx->domain_count[v] == 0u ||
                combinations > N5_STEP11_FACTOR_CACHE_MAX_VARIANTS / ctx->domain_count[v])
                break;
            combinations *= ctx->domain_count[v];
        }
        if (k != fixed || ctx->slice_plan->raw_assignment_count <= combinations)
            continue;
        /* Bound simultaneous per-worker memory even in pathological networks. */
        budget = (uint64_t)combinations * (uint64_t)src->nnz * sizeof(N5SparseEntry)
            + (uint64_t)combinations * sizeof(N5SparseFactor);
        if (budget > remaining) continue;
        cache->variants = (N5SparseFactor*)calloc(combinations, sizeof(N5SparseFactor));
        if (!cache->variants) continue;
        for (k = 0; k < combinations; ++k) {
            uint8_t values[N5_SLICING_MAX_VARIABLES] = {0};
            unsigned int code = k, j;
            for (j = 0; j < fixed; ++j) {
                unsigned int v = ctx->prepared_condition.factor[f].fixed_value_indices[j];
                unsigned int digit = code % ctx->domain_count[v];
                code /= ctx->domain_count[v];
                values[v] = ctx->domain_values[v][digit];
            }
            if (!n5_slicing_condition_prepared_factor(
                    &ctx->prepared_condition, f, values, &cache->variants[k]))
                break;
            ++cache->variant_count;
        }
        if (k != combinations) {
            unsigned int j;
            for (j = 0; j < cache->variant_count; ++j)
                n5_sparse_factor_free(&cache->variants[j]);
            free(cache->variants);
            cache->variants = NULL;
            cache->variant_count = 0u;
            continue; /* Revert to ordinary exact conditioning. */
        }
        remaining -= budget;
    }
}

/* Returns an immutable/owned mixture.  The caller must use free_view, NEVER
 * n5_sparse_network_free, because cached factors and source aliases are shared.
 */
static int n5_step11_condition_view(
    N5OddAdaptiveSliceContext* ctx,
    N5SparseNetwork* view,
    uint8_t owned[N5_SPARSE_NETWORK_MAX_FACTORS])
{
    unsigned int f;
    n5_sparse_network_init(view);
    memset(owned, 0, N5_SPARSE_NETWORK_MAX_FACTORS);
    for (f = 0; f < ctx->network->factor_count; ++f) {
        unsigned int fixed = ctx->prepared_condition.factor[f].fixed_count;
        N5SparseFactor* dst = &view->factors[f];
        if (fixed == 0u) {
            *dst = ctx->network->factors[f];
        } else if (ctx->factor_cache[f].variants != NULL) {
            unsigned int j, idx = 0u, stride = 1u;
            for (j = 0; j < fixed; ++j) {
                unsigned int v = ctx->prepared_condition.factor[f].fixed_value_indices[j];
                idx += ctx->current_indices[v] * stride;
                stride *= ctx->domain_count[v];
            }
            if (idx >= ctx->factor_cache[f].variant_count) return 0;
            *dst = ctx->factor_cache[f].variants[idx];
        } else {
            if (!n5_slicing_condition_prepared_factor(
                    &ctx->prepared_condition, f, ctx->current_values, dst))
                return 0;
            owned[f] = 1u;
        }
        view->factor_count = f + 1u;
        if (dst->nnz == 0u) break; /* exact absorbing-zero shortcut */
    }
    return 1;
}

static void n5_step11_free_view(
    N5SparseNetwork* view,
    const uint8_t owned[N5_SPARSE_NETWORK_MAX_FACTORS])
{
    unsigned int f;
    for (f = 0; f < view->factor_count; ++f) {
        if (owned[f]) n5_sparse_factor_free(&view->factors[f]);
    }
    view->factor_count = 0u;
}
#endif

static void n5_odd_adaptive_slice_recurse(
    N5OddAdaptiveSliceContext* context,
    unsigned int depth
)
{
    if (context == NULL || !context->success)
    {
        return;
    }

    if (depth == context->slice_plan->variable_count)
    {
        N5SparseNetwork conditioned;
        N5SparseContractionPlan plan;
        N5DenseMultiModularContractionResult dense_result;
        unsigned int prime_index;

        ++context->visited_assignments;
        n5_sparse_network_init(&conditioned);
#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
        {
            uint8_t owned[N5_SPARSE_NETWORK_MAX_FACTORS];
#ifdef N5_STEP11_PROFILE_TIMING
            double a = n5_parallel_wall_time_seconds();
#endif
            if (!n5_step11_condition_view(context, &conditioned, owned)) {
                n5_step11_free_view(&conditioned, owned);
                context->success = 0;
                return;
            }
            /* The borrowed view must remain live until dense contraction.
               Ownership flags are needed through the entire leaf scope. */
            if (n5_odd_adaptive_network_is_zero(&conditioned)) {
                ++context->zero_assignments;
                n5_step11_free_view(&conditioned, owned);
                return;
            }
#ifdef N5_STEP11_PROFILE_TIMING
            context->condition_seconds += n5_parallel_wall_time_seconds() - a;
            a = n5_parallel_wall_time_seconds();
#endif
            if (!n5_sparse_choose_plan(&conditioned, &plan)) {
                context->success = 0;
                n5_step11_free_view(&conditioned, owned);
                return;
            }
#ifdef N5_STEP11_PROFILE_TIMING
            context->plan_seconds += n5_parallel_wall_time_seconds() - a;
            a = n5_parallel_wall_time_seconds();
#endif
            if (!n5_dense_modular_contract_sparse_network_multi(
                    &conditioned, &plan, n5_crt_prime_basis(),
                    N5_CRT_PRIME_COUNT, &dense_result)) {
                context->success = 0;
                n5_step11_free_view(&conditioned, owned);
                return;
            }
#ifdef N5_STEP11_PROFILE_TIMING
            context->dense_seconds += n5_parallel_wall_time_seconds() - a;
#endif
            n5_step11_free_view(&conditioned, owned);
        }
#else
        n5_sparse_network_init(&conditioned);

#ifdef N5_STEP11_DISABLE_COMPILED_CONDITION
        if (!n5_slicing_condition_sparse_network(
                context->network,
                context->slice_plan->variables,
                context->current_values,
                context->slice_plan->variable_count,
                &conditioned
            ))
#else
        if (!n5_slicing_condition_prepared(
                &context->prepared_condition,
                context->current_values,
                &conditioned
            ))
#endif
        {
            context->success = 0;
            return;
        }

        if (n5_odd_adaptive_network_is_zero(&conditioned))
        {
            ++context->zero_assignments;
            n5_sparse_network_free(&conditioned);
            return;
        }

        if (!n5_sparse_choose_plan(&conditioned, &plan) ||
            !n5_dense_modular_contract_sparse_network_multi(
                &conditioned,
                &plan,
                n5_crt_prime_basis(),
                N5_CRT_PRIME_COUNT,
                &dense_result
            ))
        {
            n5_sparse_network_free(&conditioned);
            context->success = 0;
            return;
        }

#endif

        if (dense_result.max_table_entries > context->max_dense_table_entries)
        {
            context->max_dense_table_entries = dense_result.max_table_entries;
        }

        for (prime_index = 0u;
             prime_index < N5_CRT_PRIME_COUNT;
             ++prime_index)
        {
            uint32_t modulus;

            modulus = n5_crt_prime_basis()[prime_index];
            context->residues[prime_index] = (uint32_t)(
                ((uint64_t)context->residues[prime_index] +
                 dense_result.values[prime_index]) % modulus
            );
        }

#if defined(N5_STEP11_DISABLE_COMPILED_CONDITION) || defined(N5_STEP11_DISABLE_FACTOR_CACHE)
        n5_sparse_network_free(&conditioned);
#endif
        return;
    }

    {
        unsigned int value_index;

        for (value_index = 0u;
             value_index < context->domain_count[depth];
             ++value_index)
        {
            context->current_values[depth] =
                context->domain_values[depth][value_index];
#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
            context->current_indices[depth] = value_index;
#endif
            n5_odd_adaptive_slice_recurse(context, depth + 1u);
            if (!context->success)
            {
                return;
            }
        }
    }
}

/* sliced multi-prime contraction helper。 */
static int n5_odd_adaptive_contract_sliced_multi(
    const N5SparseNetwork* network,
    const N5SlicePlan* slice_plan,
    uint32_t residues[N5_CRT_PRIME_COUNT],
    uint64_t* raw_assignments,
    uint64_t* zero_assignments,
    uint64_t* max_dense_table_entries
)
{
    N5OddAdaptiveSliceContext context;
    unsigned int index;

    if (network == NULL || slice_plan == NULL || residues == NULL ||
        raw_assignments == NULL || zero_assignments == NULL ||
        max_dense_table_entries == NULL)
    {
        return 0;
    }

    memset(&context, 0, sizeof(context));
    context.network = network;
    context.slice_plan = slice_plan;
    context.success = 1;

    for (index = 0u; index < slice_plan->variable_count; ++index)
    {
        if (!n5_odd_adaptive_collect_variable_domain(
                network,
                slice_plan->variables[index],
                context.domain_values[index],
                &context.domain_count[index]
            ))
        {
            return 0;
        }

    }
#ifndef N5_STEP11_DISABLE_COMPILED_CONDITION
    if (!n5_slicing_prepare_condition(
            network,
            slice_plan->variables,
            slice_plan->variable_count,
            &context.prepared_condition
        ))
    {
        return 0;
    }
#endif

#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
#ifdef N5_STEP11_PROFILE_TIMING
    {
        double t = n5_parallel_wall_time_seconds();
        n5_step11_cache_prepare(&context);
        context.precompute_seconds += n5_parallel_wall_time_seconds() - t;
    }
#else
    n5_step11_cache_prepare(&context);
#endif
#endif
    n5_odd_adaptive_slice_recurse(&context, 0u);
#if !defined(N5_STEP11_DISABLE_COMPILED_CONDITION) && !defined(N5_STEP11_DISABLE_FACTOR_CACHE)
    n5_step11_cache_release(&context);
#endif

#ifdef N5_STEP11_PROFILE_TIMING
    fprintf(stderr, "STEP11_PROFILE sliced=%llu zero=%llu precompute=%.6f condition=%.6f plan=%.6f dense=%.6f\n",
        (unsigned long long)context.visited_assignments,
        (unsigned long long)context.zero_assignments,
        context.precompute_seconds, context.condition_seconds,
        context.plan_seconds, context.dense_seconds);
#endif
    if (!context.success)
    {
        return 0;
    }

    memcpy(
        residues,
        context.residues,
        sizeof(context.residues)
    );

    *raw_assignments = context.visited_assignments;
    *zero_assignments = context.zero_assignments;
    *max_dense_table_entries = context.max_dense_table_entries;

    return 1;
}

int n5_odd_adaptive_contract_fixed_charge(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[16],
    const N5OddKernelLibrary* library,
    uint64_t target_peak,
    unsigned int max_slices,
    N5OddAdaptiveClassResult* result
)
{
    N5OddRepresentationNetwork ordinary;
    N5SparseContractionPlan ordinary_plan;

    N5SparseNetwork rewritten;
    N5SparseContractionPlan rewritten_plan;
    int rewritten_built;
    int rewrite_selected;

    const N5SparseNetwork* selected_network;
    const N5SparseContractionPlan* selected_plan;

    unsigned int cell;
    int has_k2;

    if (size_matrix == NULL || r_by_cell == NULL || library == NULL ||
        result == NULL || target_peak == 0u ||
        max_slices > N5_SLICING_MAX_VARIABLES)
    {
        return 0;
    }

    memset(result, 0, sizeof(*result));
    n5_odd_network_init(&ordinary);
    n5_sparse_network_init(&rewritten);
    rewritten_built = 0;
    rewrite_selected = 0;

    if (!n5_odd_network_build(size_matrix, r_by_cell, library, &ordinary))
    {
        n5_odd_network_free(&ordinary);
        return 0;
    }

    result->original_factor_count = ordinary.sparse.factor_count;
    result->original_primitive_nnz = ordinary.primitive_nnz;

    if (n5_odd_adaptive_network_is_zero(&ordinary.sparse))
    {
        result->success = 1;
        result->mode = N5_ODD_ADAPTIVE_MODE_ZERO;
        n5_odd_network_free(&ordinary);
        return 1;
    }

    if (!n5_sparse_choose_plan(&ordinary.sparse, &ordinary_plan))
    {
        n5_odd_network_free(&ordinary);
        return 0;
    }

    result->original_peak = ordinary_plan.estimated_peak;

    selected_network = &ordinary.sparse;
    selected_plan = &ordinary_plan;

    has_k2 = 0;
    for (cell = 0u; cell < 16u; ++cell)
    {
        if (n5_odd_adaptive_cell_size(size_matrix, cell) == 2u)
        {
            has_k2 = 1;
            break;
        }
    }

    if (ordinary_plan.estimated_peak > target_peak && has_k2)
    {
        if (!n5_odd_adaptive_build_second_walsh_network(
                size_matrix,
                r_by_cell,
                library,
                &ordinary.sparse,
                &rewritten
            ))
        {
            n5_odd_network_free(&ordinary);
            return 0;
        }

        rewritten_built = 1;

        if (n5_odd_adaptive_network_is_zero(&rewritten))
        {
            result->success = 1;
            result->mode = N5_ODD_ADAPTIVE_MODE_ZERO;
            result->rewrite_applied = 1;
            result->selected_factor_count = rewritten.factor_count;
            result->selected_primitive_nnz =
                n5_sparse_network_primitive_nnz(&rewritten);
            n5_sparse_network_free(&rewritten);
            n5_odd_network_free(&ordinary);
            return 1;
        }

        if (!n5_sparse_choose_plan(&rewritten, &rewritten_plan))
        {
            n5_sparse_network_free(&rewritten);
            n5_odd_network_free(&ordinary);
            return 0;
        }

        if (rewritten_plan.estimated_peak < ordinary_plan.estimated_peak)
        {
            selected_network = &rewritten;
            selected_plan = &rewritten_plan;
            rewrite_selected = 1;
        }
    }

    result->rewrite_applied = rewrite_selected;
    result->selected_factor_count = selected_network->factor_count;
    result->selected_primitive_nnz =
        n5_sparse_network_primitive_nnz(selected_network);
    result->selected_peak = selected_plan->estimated_peak;

    if (selected_plan->estimated_peak <= target_peak)
    {
        N5DenseMultiModularContractionResult dense_result;

        if (!n5_dense_modular_contract_sparse_network_multi(
                selected_network,
                selected_plan,
                n5_crt_prime_basis(),
                N5_CRT_PRIME_COUNT,
                &dense_result
            ))
        {
            if (rewritten_built)
            {
                n5_sparse_network_free(&rewritten);
            }
            n5_odd_network_free(&ordinary);
            return 0;
        }

        memcpy(
            result->residues,
            dense_result.values,
            sizeof(result->residues)
        );

        result->raw_slice_assignments = 1u;
        result->sliced_peak = selected_plan->estimated_peak;
        result->max_dense_table_entries = dense_result.max_table_entries;
        result->mode = rewrite_selected
            ? N5_ODD_ADAPTIVE_MODE_DIRECT_REWRITTEN
            : N5_ODD_ADAPTIVE_MODE_DIRECT;
        result->success = 1;
    }
    else
    {
        N5SlicePlan slice_plan;
        uint64_t raw_assignments;
        uint64_t zero_assignments;
        uint64_t max_dense_entries;

        if (!n5_slicing_choose_greedy(
                selected_network,
                target_peak,
                max_slices,
                &slice_plan
            ) ||
            !slice_plan.target_reached)
        {
            if (rewritten_built)
            {
                n5_sparse_network_free(&rewritten);
            }
            n5_odd_network_free(&ordinary);
            return 0;
        }

        if (!n5_odd_adaptive_contract_sliced_multi(
                selected_network,
                &slice_plan,
                result->residues,
                &raw_assignments,
                &zero_assignments,
                &max_dense_entries
            ))
        {
            if (rewritten_built)
            {
                n5_sparse_network_free(&rewritten);
            }
            n5_odd_network_free(&ordinary);
            return 0;
        }

        result->slice_variable_count = slice_plan.variable_count;
        memcpy(
            result->slice_variables,
            slice_plan.variables,
            sizeof(result->slice_variables)
        );
        result->raw_slice_assignments = raw_assignments;
        result->zero_slice_assignments = zero_assignments;
        result->sliced_peak = slice_plan.estimated_peak;
        result->max_dense_table_entries = max_dense_entries;
        result->mode = rewrite_selected
            ? N5_ODD_ADAPTIVE_MODE_SLICED_REWRITTEN
            : N5_ODD_ADAPTIVE_MODE_SLICED;
        result->success = 1;
    }

    if (rewritten_built)
    {
        n5_sparse_network_free(&rewritten);
    }
    n5_odd_network_free(&ordinary);

    return result->success;
}

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "n5_dense.h"
#include "n5_modular.h"

/* Step11: exact reciprocal reduction, x mod p, via the high half of 64x64.
 * Always retain a portable modulo fallback.  This does not assume prime p.
 * __umulh is an x64 MSVC intrinsic; __uint128_t is available on GCC/Clang.
 */
#if !defined(N5_STEP11_DISABLE_RECIPROCAL_MOD) && \
   ((defined(_MSC_VER) && defined(_M_X64)) || defined(__SIZEOF_INT128__))
#define N5_STEP11_FASTMOD_SUPPORTED 1
#endif
#if defined(N5_STEP11_FASTMOD_SUPPORTED) && defined(_MSC_VER)
#include <intrin.h>
#endif

#ifdef N5_STEP11_FASTMOD_SUPPORTED
static inline uint32_t n5_step11_reduce_u64(
    uint64_t x, uint32_t p, uint64_t reciprocal)
{
    uint64_t q, r;
#if defined(_MSC_VER) && defined(_M_X64)
    q = __umulh(x, reciprocal);
#else
    q = (uint64_t)(((__uint128_t)x * (__uint128_t)reciprocal) >> 64);
#endif
    r = x - q * (uint64_t)p;
    /* reciprocal=floor((2^64-1)/p), hence floor(x/p)-q in {0,1}
       for all 0 <= x < 2^64, p>=2.  One subtraction is sufficient. */
    if (r >= (uint64_t)p) r -= (uint64_t)p;
    return (uint32_t)r;
}
#define N5_STEP11_REDUCE(x,p,r) n5_step11_reduce_u64((x),(p),(r))
#else
#define N5_STEP11_REDUCE(x,p,r) ((uint32_t)((x) % (p)))
#endif

#if !defined(N5_STEP10_DISABLE_TRIPLE_PAIR_CACHE) && \
    !defined(N5_STEP6_DISABLE_TRIPLE_KERNEL)
#define N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE 1
#endif

/*
 * 单变量 compact domain。
 */
typedef struct
{
    uint16_t count;
    uint8_t values[256];
    int16_t index_of[256];
} N5DenseVariableDomain;

/*
 * 一个 dense modular factor。
 *
 * table index 采用 factor.vars 顺序的 mixed-radix 编码。
 */
typedef struct
{
    uint8_t arity;
    uint8_t vars[N5_SPARSE_MAX_ARITY];
    uint16_t domain_sizes[N5_SPARSE_MAX_ARITY];
    uint64_t strides[N5_SPARSE_MAX_ARITY];
    uint64_t table_size;
    uint32_t* values;
} N5DenseModularFactor;

/*
 * 模加法。
 */
static uint32_t n5_dense_mod_add(
    uint32_t a,
    uint32_t b,
    uint32_t modulus
)
{
    uint64_t total;

    /* The operands are canonical residues: 0 <= a,b < modulus.
       Their sum is below 2*modulus; subtract once instead of division.
       Use uint64_t to retain correctness even near UINT32_MAX. */
    total = (uint64_t)a + (uint64_t)b;
    if (total >= (uint64_t)modulus)
    {
        total -= (uint64_t)modulus;
    }

    return (uint32_t)total;
}

/*
 * 模乘法。
 */
static uint32_t n5_dense_mod_mul(
    uint32_t a,
    uint32_t b,
    uint32_t modulus
)
{
    return (uint32_t)(
        ((uint64_t)a * (uint64_t)b) % modulus
    );
}

/*
 * 从 sparse network 恢复所有 variable 的实际 value domains。
 */
static int n5_dense_build_domains(
    const N5SparseNetwork* network,
    N5DenseVariableDomain domains[N5_SPARSE_MAX_VARIABLES]
)
{
    unsigned int var;
    unsigned int factor_index;

    if (network == NULL || domains == NULL)
    {
        return 0;
    }

    for (var = 0; var < N5_SPARSE_MAX_VARIABLES; ++var)
    {
        unsigned int value;

        domains[var].count = 0;

        for (value = 0; value < 256; ++value)
        {
            domains[var].index_of[value] = -1;
        }
    }

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        unsigned int entry_index;

        factor = &network->factors[factor_index];

        for (entry_index = 0;
             entry_index < factor->nnz;
             ++entry_index)
        {
            unsigned int position;

            for (position = 0;
                 position < factor->arity;
                 ++position)
            {
                unsigned int variable;
                unsigned int value;

                variable = factor->vars[position];
                value = factor->entries[entry_index].assignment[position];

                if (variable >= N5_SPARSE_MAX_VARIABLES)
                {
                    return 0;
                }

                if (domains[variable].index_of[value] < 0)
                {
                    uint16_t slot;

                    slot = domains[variable].count;

                    if (slot >= 256)
                    {
                        return 0;
                    }

                    domains[variable].values[slot] = (uint8_t)value;
                    domains[variable].index_of[value] = (int16_t)slot;
                    domains[variable].count = (uint16_t)(slot + 1u);
                }
            }
        }
    }

    /*
     * canonical domain order = ascending concrete value。
     */
    for (var = 0; var < N5_SPARSE_MAX_VARIABLES; ++var)
    {
        uint16_t count;
        uint16_t i;

        count = domains[var].count;

        for (i = 1; i < count; ++i)
        {
            uint8_t current;
            uint16_t j;

            current = domains[var].values[i];
            j = i;

            while (j > 0 &&
                   current < domains[var].values[j - 1])
            {
                domains[var].values[j] = domains[var].values[j - 1];
                --j;
            }

            domains[var].values[j] = current;
        }

        for (i = 0; i < 256; ++i)
        {
            domains[var].index_of[i] = -1;
        }

        for (i = 0; i < count; ++i)
        {
            domains[var].index_of[
                domains[var].values[i]
            ] = (int16_t)i;
        }
    }

    return 1;
}

/*
 * 查找 variable 在 scope 中的位置。
 */
static int n5_dense_scope_position(
    const uint8_t* vars,
    unsigned int arity,
    unsigned int variable
)
{
    unsigned int position;

    for (position = 0; position < arity; ++position)
    {
        if (vars[position] == variable)
        {
            return (int)position;
        }
    }

    return -1;
}

/*
 * 判断 network 是否含显式 zero factor。
 */
static int n5_dense_network_is_zero(
    const N5SparseNetwork* network
)
{
    unsigned int factor_index;

    if (network == NULL)
    {
        return 1;
    }

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        if (network->factors[factor_index].nnz == 0)
        {
            return 1;
        }
    }

    return 0;
}

/*
 * 收集一个 sparse variable 的实际 domain。
 */
static int n5_dense_collect_variable_domain(
    const N5SparseNetwork* network,
    unsigned int variable,
    uint8_t values[256],
    unsigned int* value_count
)
{
    int seen[256];
    unsigned int factor_index;
    unsigned int count;
    unsigned int value;

    if (network == NULL || values == NULL || value_count == NULL)
    {
        return 0;
    }

    for (value = 0; value < 256; ++value)
    {
        seen[value] = 0;
    }

    count = 0;

    for (factor_index = 0;
         factor_index < network->factor_count;
         ++factor_index)
    {
        const N5SparseFactor* factor;
        int position;
        unsigned int entry_index;

        factor = &network->factors[factor_index];
        position = n5_dense_scope_position(
            factor->vars,
            factor->arity,
            variable
        );

        if (position < 0)
        {
            continue;
        }

        for (entry_index = 0;
             entry_index < factor->nnz;
             ++entry_index)
        {
            value = factor->entries[entry_index].assignment[position];

            if (!seen[value])
            {
                seen[value] = 1;
            }
        }
    }

    for (value = 0; value < 256; ++value)
    {
        if (seen[value])
        {
            values[count] = (uint8_t)value;
            ++count;
        }
    }

    *value_count = count;

    return count > 0;
}

/* ===================================================================== */
/* Shared-structure multi-modulus dense backend. */
/* ===================================================================== */

/*
 * 多模数 dense factor。
 *
 * values 的布局：
 *
 *   values[table_index * modulus_count + lane]
 *
 * 每个 lane 使用自己独立的 modulus。
 */
typedef struct
{
    uint8_t arity;
    uint8_t vars[N5_SPARSE_MAX_ARITY];
    uint16_t domain_sizes[N5_SPARSE_MAX_ARITY];
    uint64_t strides[N5_SPARSE_MAX_ARITY];
    uint64_t table_size;

    uint8_t modulus_count;
    uint32_t* values;
} N5DenseMultiFactor;

/*
 * 释放 multi-modulus dense factor。
 */
static void n5_dense_multi_factor_free(
    N5DenseMultiFactor* factor
)
{
    if (factor == NULL)
    {
        return;
    }

    free(factor->values);
    factor->values = NULL;
    factor->arity = 0;
    factor->table_size = 0;
    factor->modulus_count = 0;
}

/*
 * 初始化 multi-modulus dense factor storage。
 */
/*
 * zero_initialize=1 is required when expanding a sparse factor, because
 * absent assignments represent zero. For an elimination output, pass 0:
 * the contraction loop writes every (table_index, lane) before any read.
 */
static int n5_dense_multi_factor_allocate(
    const uint8_t* vars,
    unsigned int arity,
    const N5DenseVariableDomain domains[N5_SPARSE_MAX_VARIABLES],
    unsigned int modulus_count,
    int zero_initialize,
    N5DenseMultiFactor* factor
)
{
    uint64_t size;
    uint64_t scalar_count;
    unsigned int position;

    if (factor == NULL ||
        modulus_count == 0 ||
        modulus_count > N5_DENSE_MAX_MODULI ||
        arity > N5_SPARSE_MAX_ARITY ||
        (arity > 0 && vars == NULL))
    {
        return 0;
    }

    memset(factor, 0, sizeof(*factor));

    factor->arity = (uint8_t)arity;
    factor->modulus_count = (uint8_t)modulus_count;
    size = 1;

    for (position = 0; position < arity; ++position)
    {
        factor->vars[position] = vars[position];

        if (vars[position] >= N5_SPARSE_MAX_VARIABLES ||
            domains[vars[position]].count == 0)
        {
            return 0;
        }

        factor->domain_sizes[position] =
            domains[vars[position]].count;
    }

    if (arity > 0)
    {
        int pos;

        factor->strides[arity - 1u] = 1;

        for (pos = (int)arity - 2; pos >= 0; --pos)
        {
            factor->strides[pos] =
                factor->strides[pos + 1] *
                (uint64_t)factor->domain_sizes[pos + 1];
        }

        for (position = 0; position < arity; ++position)
        {
            if (size > UINT64_MAX /
                (uint64_t)factor->domain_sizes[position])
            {
                return 0;
            }

            size *= (uint64_t)factor->domain_sizes[position];
        }
    }

    factor->table_size = size;

    if (size > UINT64_MAX / (uint64_t)modulus_count)
    {
        return 0;
    }

    scalar_count = size * (uint64_t)modulus_count;

    if (scalar_count > (uint64_t)SIZE_MAX / sizeof(uint32_t))
    {
        return 0;
    }

    factor->values = zero_initialize
        ? (uint32_t*)calloc((size_t)scalar_count, sizeof(uint32_t))
        : (uint32_t*)malloc((size_t)scalar_count * sizeof(uint32_t));

    return factor->values != NULL;
}

/*
 * exact sparse factor -> multi-modulus dense factor。
 *
 * 每个 exact rational coefficient 独立映射到每个 Z/pZ lane。
 */
static int n5_dense_multi_factor_from_sparse(
    const N5SparseFactor* sparse,
    const N5DenseVariableDomain domains[N5_SPARSE_MAX_VARIABLES],
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5DenseMultiFactor* dense
)
{
    unsigned int entry_index;

    if (sparse == NULL || dense == NULL ||
        moduli == NULL ||
        modulus_count == 0 ||
        modulus_count > N5_DENSE_MAX_MODULI)
    {
        return 0;
    }

    if (!n5_dense_multi_factor_allocate(
            sparse->vars,
            sparse->arity,
            domains,
            modulus_count,
            1,
            dense
        ))
    {
        return 0;
    }

    for (entry_index = 0;
         entry_index < sparse->nnz;
         ++entry_index)
    {
        uint64_t table_index;
        unsigned int position;
        unsigned int lane;

        table_index = 0;

        for (position = 0;
             position < sparse->arity;
             ++position)
        {
            unsigned int variable;
            unsigned int value;
            int local_index;

            variable = sparse->vars[position];
            value = sparse->entries[entry_index].assignment[position];
            local_index = domains[variable].index_of[value];

            if (local_index < 0)
            {
                n5_dense_multi_factor_free(dense);
                return 0;
            }

            table_index +=
                (uint64_t)local_index * dense->strides[position];
        }

        for (lane = 0; lane < modulus_count; ++lane)
        {
            uint32_t residue;
            uint64_t flat_index;

            if (moduli[lane] <= 1 ||
                !n5_modular_rational_to_residue(
                    sparse->entries[entry_index].value,
                    moduli[lane],
                    &residue
                ))
            {
                n5_dense_multi_factor_free(dense);
                return 0;
            }

            flat_index =
                table_index * (uint64_t)modulus_count +
                (uint64_t)lane;

            dense->values[flat_index] =
                n5_dense_mod_add(
                    dense->values[flat_index],
                    residue,
                    moduli[lane]
                );
        }
    }

    return 1;
}

/*
 * multi-modulus dense factor scope lookup。
 */
static int n5_dense_multi_scope_position(
    const uint8_t* vars,
    unsigned int arity,
    unsigned int variable
)
{
    return n5_dense_scope_position(
        vars,
        arity,
        variable
    );
}

/*
 * 在所有 modulus lanes 中同时消去一个 variable。
 *
 * 纯结构工作：
 *
 *   - incident-factor discovery；
 *   - output scope；
 *   - mixed-radix index mapping；
 *   - output assignment decoding；
 *
 * 只执行一次。
 *
 * arithmetic work 在每个 lane 内独立执行。
 */
static int n5_dense_multi_eliminate_variable(
    N5DenseMultiFactor* factors,
    unsigned int* factor_count,
    unsigned int variable,
    const N5DenseVariableDomain domains[N5_SPARSE_MAX_VARIABLES],
    const uint32_t* moduli,
    unsigned int modulus_count,
    uint64_t* max_table_entries
)
{
    unsigned int incident_indices[N5_SPARSE_NETWORK_MAX_FACTORS];
    unsigned int incident_count;

    uint8_t output_vars[N5_SPARSE_MAX_ARITY];
    unsigned int output_arity;

    uint64_t factor_output_stride[
        N5_SPARSE_NETWORK_MAX_FACTORS
    ][N5_SPARSE_MAX_ARITY];

    unsigned int factor_index;
    N5DenseMultiFactor output;
#ifdef N5_STEP11_FASTMOD_SUPPORTED
    uint64_t step11_reciprocal[N5_DENSE_MAX_MODULI];
#endif
#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
    uint64_t *step10_pair_products = NULL;
    uint64_t step10_pair_output_stride[N5_SPARSE_MAX_ARITY] = {0};
    uint64_t step10_pair_output_index = 0;
    unsigned int step10_pair_a = 0, step10_pair_b = 0, step10_pair_c = 0;
#endif
#ifndef N5_STEP6_DISABLE_TRIPLE_KERNEL
    int triple_product_safe = 1;
    int triple_sum_safe = 1;
#endif
#ifndef N5_STEP10_DISABLE_WIDE_LOWORDER
    int loworder_sum_safe = 1;
#endif

    if (factors == NULL || factor_count == NULL ||
        moduli == NULL ||
        modulus_count == 0 ||
        modulus_count > N5_DENSE_MAX_MODULI ||
        variable >= N5_SPARSE_MAX_VARIABLES ||
        domains[variable].count == 0)
    {
        return 0;
    }

#ifndef N5_STEP6_DISABLE_TRIPLE_KERNEL
    /* Check three-factor multiplication and entire eliminated-domain sum
       against 64-bit overflow. The six production primes are ~10^6. */
    {
        unsigned int lane;
        for (lane = 0; lane < modulus_count; ++lane)
        {
            const uint64_t bound = (uint64_t)moduli[lane] - 1u;
            if (moduli[lane] <= 1u ||
                (bound != 0u && bound * bound > UINT64_MAX / bound))
            {
                triple_product_safe = 0;
                triple_sum_safe = 0;
                break;
            }
            else if (bound != 0u &&
                     bound * bound * bound >
                     UINT64_MAX / (uint64_t)domains[variable].count)
            {
                triple_sum_safe = 0;
            }
        }
    }

#endif /* N5_STEP6_DISABLE_TRIPLE_KERNEL */

#ifndef N5_STEP10_DISABLE_WIDE_LOWORDER
    /* For 1- or 2-incident eliminations, exact accumulation of the full
       eliminated domain fits uint64_t if (p-1)^2 * domain_size <= UINT64_MAX.
       Production CRT primes (~1e6) satisfy this for every compact domain.
       For larger user-supplied moduli, retain the generic modular path. */
    {
        unsigned int q;
        for (q = 0; q < modulus_count; ++q) {
            uint64_t bound = (uint64_t)moduli[q] - 1u;
            if (moduli[q] <= 1u ||
                bound * bound > UINT64_MAX / (uint64_t)domains[variable].count) {
                loworder_sum_safe = 0;
                break;
            }
        }
    }
#endif

#ifdef N5_STEP11_FASTMOD_SUPPORTED
    {
        unsigned int i;
        for (i = 0; i < modulus_count; ++i) {
            if (moduli[i] <= 1u) return 0;
            step11_reciprocal[i] = UINT64_MAX / (uint64_t)moduli[i];
        }
    }
#endif

    incident_count = 0;
    output_arity = 0;

    for (factor_index = 0;
         factor_index < *factor_count;
         ++factor_index)
    {
        N5DenseMultiFactor* factor;
        int position;

        factor = &factors[factor_index];

        position = n5_dense_multi_scope_position(
            factor->vars,
            factor->arity,
            variable
        );

        if (position >= 0)
        {
            unsigned int p;

            incident_indices[incident_count] = factor_index;
            ++incident_count;

            for (p = 0; p < factor->arity; ++p)
            {
                unsigned int v;

                v = factor->vars[p];

                if (v == variable)
                {
                    continue;
                }

                if (n5_dense_multi_scope_position(
                        output_vars,
                        output_arity,
                        v
                    ) < 0)
                {
                    if (output_arity >= N5_SPARSE_MAX_ARITY)
                    {
                        return 0;
                    }

                    output_vars[output_arity] = (uint8_t)v;
                    ++output_arity;
                }
            }
        }
    }

    if (incident_count == 0)
    {
        return 1;
    }

    /* canonical output scope order */
    {
        unsigned int i;

        for (i = 1; i < output_arity; ++i)
        {
            uint8_t current;
            unsigned int j;

            current = output_vars[i];
            j = i;

            while (j > 0 && current < output_vars[j - 1])
            {
                output_vars[j] = output_vars[j - 1];
                --j;
            }

            output_vars[j] = current;
        }
    }

    if (!n5_dense_multi_factor_allocate(
            output_vars,
            output_arity,
            domains,
            modulus_count,
            0,
            &output
        ))
    {
        return 0;
    }

    if (output.table_size > *max_table_entries)
    {
        *max_table_entries = output.table_size;
    }

    /*
     * Precompute each incident factor's offset increment for every output
     * coordinate. The output table is traversed in mixed-radix order, so
     * adjacent output indices differ in exactly one coordinate (plus
     * occasional carry). We can update offsets incrementally rather than
     * decoding every output index with division/modulo and rebuilding each
     * factor offset from scratch.
     */
    {
        unsigned int incident;

        memset(factor_output_stride, 0, sizeof(factor_output_stride));

        for (incident = 0; incident < incident_count; ++incident)
        {
            const N5DenseMultiFactor* factor = &factors[incident_indices[incident]];
            unsigned int fp;

            for (fp = 0; fp < factor->arity; ++fp)
            {
                int pos;
                if (factor->vars[fp] == variable)
                    continue;

                pos = n5_dense_multi_scope_position(
                    output.vars, output.arity, factor->vars[fp]);
                if (pos < 0)
                {
                    n5_dense_multi_factor_free(&output);
                    return 0;
                }
                factor_output_stride[incident][(unsigned int)pos] =
                    factor->strides[fp];
            }
        }
    }

    {
        uint64_t out_index;
        uint16_t assignment_indices[N5_SPARSE_MAX_ARITY];
        uint64_t incident_offsets[N5_SPARSE_NETWORK_MAX_FACTORS];
        uint64_t incident_elimination_strides[N5_SPARSE_NETWORK_MAX_FACTORS];
        unsigned int incident;

        memset(assignment_indices, 0, sizeof(assignment_indices));
        memset(incident_offsets, 0, sizeof(incident_offsets));

        /* Fixed for the entire elimination: x advances by its factor stride. */
        for (incident = 0; incident < incident_count; ++incident)
        {
            const N5DenseMultiFactor* factor = &factors[incident_indices[incident]];
            unsigned int fp;
            incident_elimination_strides[incident] = 0;
            for (fp = 0; fp < factor->arity; ++fp)
            {
                if (factor->vars[fp] == variable)
                {
                    incident_elimination_strides[incident] = factor->strides[fp];
                    break;
                }
            }
        }

#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
    /* Tensor-network pair-product memoization: when two of the three
     * factors do not depend on every output coordinate, their exact product
     * can be reused for all assignments to the omitted coordinates.
     * Only use this path when the ENTIRE 3-factor sum fits uint64_t.
     * Memory is capped at 16 MiB per active elimination. */
#ifndef N5_STEP6_DISABLE_TRIPLE_KERNEL
    if (incident_count == 3u && triple_sum_safe && output.table_size >= 256u)
    {
        unsigned int a, b, p;
        uint64_t best_reuse = 0u, best_pair_size = 0u;
        for (a = 0; a < 3u; ++a) for (b = a + 1u; b < 3u; ++b)
        {
            uint64_t pair_size = 1u;
            for (p = 0; p < output.arity; ++p)
                if (factor_output_stride[a][p] || factor_output_stride[b][p])
                    pair_size *= output.domain_sizes[p];
            if (pair_size && output.table_size / pair_size > best_reuse)
            {
                best_reuse = output.table_size / pair_size;
                best_pair_size = pair_size;
                step10_pair_a = a;
                step10_pair_b = b;
            }
        }
        if (best_reuse >= 4u &&
            best_pair_size <= UINT64_MAX / domains[variable].count &&
            best_pair_size * domains[variable].count <=
              (UINT64_C(16) * 1024u * 1024u) /
              (sizeof(uint64_t) * modulus_count))
        {
            uint64_t compact_stride = 1u, j;
            unsigned int q;
            step10_pair_c = 3u - step10_pair_a - step10_pair_b;
            for (p = output.arity; p > 0; --p)
            {
                unsigned int coordinate = p - 1u;
                if (factor_output_stride[step10_pair_a][coordinate] ||
                    factor_output_stride[step10_pair_b][coordinate])
                {
                    step10_pair_output_stride[coordinate] = compact_stride;
                    compact_stride *= output.domain_sizes[coordinate];
                }
            }
            step10_pair_products = (uint64_t*)malloc(
                (size_t)best_pair_size * domains[variable].count *
                modulus_count * sizeof(uint64_t));
            if (step10_pair_products != NULL)
            {
                const uint32_t *ta = factors[incident_indices[step10_pair_a]].values;
                const uint32_t *tb = factors[incident_indices[step10_pair_b]].values;
                for (j = 0; j < best_pair_size; ++j)
                {
                    uint64_t rem = j, base_a = 0u, base_b = 0u;
                    for (p = 0; p < output.arity; ++p)
                    {
                        if (step10_pair_output_stride[p])
                        {
                            const uint64_t coord = rem / step10_pair_output_stride[p];
                            rem %= step10_pair_output_stride[p];
                            base_a += coord * factor_output_stride[step10_pair_a][p];
                            base_b += coord * factor_output_stride[step10_pair_b][p];
                        }
                    }
                    for (q = 0; q < domains[variable].count; ++q)
                    {
                        uint64_t sa = (base_a + q *
                            incident_elimination_strides[step10_pair_a]) * modulus_count;
                        uint64_t sb = (base_b + q *
                            incident_elimination_strides[step10_pair_b]) * modulus_count;
                        uint64_t dest = (j * domains[variable].count + q) * modulus_count;
                        unsigned int lane;
                        for (lane = 0; lane < modulus_count; ++lane)
                            step10_pair_products[dest + lane] =
                                (uint64_t)ta[sa + lane] * (uint64_t)tb[sb + lane];
                    }
                }
            }
        }
    }
#endif
#endif

        for (out_index = 0;
             out_index < output.table_size;
             ++out_index)
        {
            uint16_t eliminated_index;
            uint32_t sums[N5_DENSE_MAX_MODULI];
            unsigned int lane;

            for (lane = 0; lane < modulus_count; ++lane)
            {
                sums[lane] = 0;
            }

#ifndef N5_STEP6_DISABLE_TRIPLE_KERNEL
            /*
             * Specialized three-incident-factor kernel. The generic kernel
             * handles any incident count; the hot all-even path often joins
             * precisely three factors. All indices are bounded by the same
             * precomputed output offsets and eliminated-variable strides.
             * Zero and unit identities are applied independently in each
             * CRT lane, so a zero in one lane never suppresses another.
             */
#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
            if (step10_pair_products != NULL)
            {
                const uint32_t * const tc = factors[incident_indices[step10_pair_c]].values;
                uint64_t pc = incident_offsets[step10_pair_c] * modulus_count;
                const uint64_t stepc =
                    incident_elimination_strides[step10_pair_c] * modulus_count;
                const uint64_t *pairs = step10_pair_products +
                    step10_pair_output_index * domains[variable].count * modulus_count;
                uint64_t acc[N5_DENSE_MAX_MODULI] = {0};
                for (eliminated_index = 0;
                     eliminated_index < domains[variable].count;
                     ++eliminated_index)
                {
                    unsigned int lane;
                    for (lane = 0; lane < modulus_count; ++lane)
                        acc[lane] += pairs[(uint64_t)eliminated_index * modulus_count + lane] *
                                     tc[pc + lane];
                    pc += stepc;
                }
                for (lane = 0; lane < modulus_count; ++lane)
                    sums[lane] = N5_STEP11_REDUCE(acc[lane], moduli[lane], step11_reciprocal[lane]);
            }
            else
#endif
            if (incident_count == 3u && triple_product_safe)
            {
                const uint32_t * table0 =
                    factors[incident_indices[0]].values;
                const uint32_t * table1 =
                    factors[incident_indices[1]].values;
                const uint32_t * table2 =
                    factors[incident_indices[2]].values;
                const uint64_t base0 =
                    incident_offsets[0] * (uint64_t)modulus_count;
                const uint64_t base1 =
                    incident_offsets[1] * (uint64_t)modulus_count;
                const uint64_t base2 =
                    incident_offsets[2] * (uint64_t)modulus_count;
                const uint64_t step0 =
                    incident_elimination_strides[0] * (uint64_t)modulus_count;
                const uint64_t step1 =
                    incident_elimination_strides[1] * (uint64_t)modulus_count;
                const uint64_t step2 =
                    incident_elimination_strides[2] * (uint64_t)modulus_count;
                uint64_t pos0 = base0, pos1 = base1, pos2 = base2;
                uint64_t accumulated[N5_DENSE_MAX_MODULI] = {0};

                for (eliminated_index = 0;
                     eliminated_index < domains[variable].count;
                     ++eliminated_index)
                {
                    for (lane = 0; lane < modulus_count; ++lane)
                    {
                        const uint32_t a = table0[pos0 + lane];
                        const uint32_t b = table1[pos1 + lane];
                        const uint32_t c = table2[pos2 + lane];
                        uint32_t product;

                        if (a == 0u || b == 0u || c == 0u)
                            continue;

                        /* If the entire eliminated-domain sum fits in
                           uint64_t, postpone all modular reductions until
                           after accumulation. Otherwise use one reduction
                           per triple as in the previous safe fast path. */
                        const uint64_t wide_product =
                            ((uint64_t)a * (uint64_t)b) * (uint64_t)c;
                        if (triple_sum_safe)
                            accumulated[lane] += wide_product;
                        else
                        {
                            product = N5_STEP11_REDUCE(wide_product, moduli[lane], step11_reciprocal[lane]);
                            sums[lane] = n5_dense_mod_add(
                                sums[lane], product, moduli[lane]);
                        }
                    }
                    pos0 += step0;
                    pos1 += step1;
                    pos2 += step2;
                }
                if (triple_sum_safe)
                {
                    for (lane = 0; lane < modulus_count; ++lane)
                        sums[lane] = N5_STEP11_REDUCE(accumulated[lane], moduli[lane], step11_reciprocal[lane]);
                }
            }
            else
#endif /* N5_STEP6_DISABLE_TRIPLE_KERNEL */
#ifndef N5_STEP10_DISABLE_WIDE_LOWORDER
            if (loworder_sum_safe && (incident_count == 1u || incident_count == 2u))
            {
                const uint32_t * const table0 =
                    factors[incident_indices[0]].values;
                const uint32_t * const table1 = incident_count == 2u
                    ? factors[incident_indices[1]].values : NULL;
                const uint64_t base0 = incident_offsets[0] * (uint64_t)modulus_count;
                const uint64_t base1 = incident_count == 2u
                    ? incident_offsets[1] * (uint64_t)modulus_count : 0u;
                const uint64_t step0 =
                    incident_elimination_strides[0] * (uint64_t)modulus_count;
                const uint64_t step1 = incident_count == 2u
                    ? incident_elimination_strides[1] * (uint64_t)modulus_count : 0u;
                uint64_t accumulated[N5_DENSE_MAX_MODULI] = {0};
                uint64_t pos0 = base0, pos1 = base1;
                if (incident_count == 1u) {
                    for (eliminated_index = 0;
                         eliminated_index < domains[variable].count;
                         ++eliminated_index) {
                        for (lane = 0; lane < modulus_count; ++lane)
                            accumulated[lane] += (uint64_t)table0[pos0 + lane];
                        pos0 += step0;
                    }
                } else {
                    for (eliminated_index = 0;
                         eliminated_index < domains[variable].count;
                         ++eliminated_index) {
                        for (lane = 0; lane < modulus_count; ++lane)
                            accumulated[lane] += (uint64_t)table0[pos0 + lane] *
                                                  (uint64_t)table1[pos1 + lane];
                        pos0 += step0;
                        pos1 += step1;
                    }
                }
                for (lane = 0; lane < modulus_count; ++lane)
                    sums[lane] = N5_STEP11_REDUCE(accumulated[lane], moduli[lane], step11_reciprocal[lane]);
            }
            else
#endif /* N5_STEP10_DISABLE_WIDE_LOWORDER */
            for (eliminated_index = 0;
                 eliminated_index < domains[variable].count;
                 ++eliminated_index)
            {
                uint32_t products[N5_DENSE_MAX_MODULI];
                unsigned int incident;

                for (lane = 0; lane < modulus_count; ++lane)
                {
                    products[lane] = 1u;
                }

                for (incident = 0;
                     incident < incident_count;
                     ++incident)
                {
                    const N5DenseMultiFactor* factor;
                    uint64_t local_index;
                    uint64_t base;
                    int all_zero;

                    factor = &factors[incident_indices[incident]];
                    local_index = incident_offsets[incident] +
                        (uint64_t)eliminated_index * incident_elimination_strides[incident];
                    base = local_index * (uint64_t)modulus_count;

                    all_zero = 1;

                    for (lane = 0; lane < modulus_count; ++lane)
                    {
                        uint32_t local_value;

                        local_value =
                            factor->values[base + lane];

                        if (local_value != 0)
                        {
                            all_zero = 0;
                        }

                        /* Exact identities in Z/pZ: 1*x=x, x*1=x,
                         * and 0*x=0. The first incident factor never
                         * requires a modular multiplication. */
                        if (incident == 0u) {
                            products[lane] = local_value;
                        } else if (products[lane] != 0u) {
                            if (local_value == 0u)
                                products[lane] = 0u;
                            else if (local_value != 1u && products[lane] != 1u)
                                products[lane] = n5_dense_mod_mul(
                                    products[lane], local_value, moduli[lane]);
                            else if (products[lane] == 1u)
                                products[lane] = local_value;
                        }
                    }

                    if (all_zero)
                    {
                        for (lane = 0;
                             lane < modulus_count;
                             ++lane)
                        {
                            products[lane] = 0;
                        }

                        break;
                    }
                }

                for (lane = 0; lane < modulus_count; ++lane)
                {
                    if (products[lane] != 0)
                    {
                        sums[lane] =
                            n5_dense_mod_add(
                                sums[lane],
                                products[lane],
                                moduli[lane]
                            );
                    }
                }
            }

            {
                uint64_t base;

                base =
                    out_index * (uint64_t)modulus_count;

                for (lane = 0; lane < modulus_count; ++lane)
                {
                    output.values[base + lane] = sums[lane];
                }
            }

            /* Advance the mixed-radix output assignment and all factor offsets. */
            {
                unsigned int pos = output.arity;
                while (pos > 0u)
                {
                    unsigned int coordinate = pos - 1u;
                    uint16_t radix = output.domain_sizes[coordinate];
                    uint16_t old_value = assignment_indices[coordinate];

                    if ((uint32_t)old_value + 1u < (uint32_t)radix)
                    {
                        assignment_indices[coordinate] = (uint16_t)(old_value + 1u);
                        for (incident = 0; incident < incident_count; ++incident)
                        {
                            incident_offsets[incident] +=
                                factor_output_stride[incident][coordinate];
                        }
#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
                        step10_pair_output_index +=
                            step10_pair_output_stride[coordinate];
#endif
                        break;
                    }
                    else
                    {
                        assignment_indices[coordinate] = 0u;
                        for (incident = 0; incident < incident_count; ++incident)
                        {
                            incident_offsets[incident] -=
                                factor_output_stride[incident][coordinate] *
                                (uint64_t)(radix - 1u);
                        }
#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
                        step10_pair_output_index -=
                            step10_pair_output_stride[coordinate] *
                            (uint64_t)(radix - 1u);
#endif
                    }
                    --pos;
                }
            }
        }
    }

#ifdef N5_STEP10_ENABLE_TRIPLE_PAIR_CACHE
    free(step10_pair_products);
#endif

    /*
     * 只有当所有 lanes 都为零时，network 才能提前折叠为零。
     *
     * CRT primes 下 exact integer/rational network 若为真正零，
     * 所有 lanes 必然为零。
     *
     * 若仅某一 lane 偶然为零，不能据此删掉其它 lanes。
     */
    {
        uint64_t index;
        int any_nonzero;

        any_nonzero = 0;

        for (index = 0;
             index < output.table_size && !any_nonzero;
             ++index)
        {
            uint64_t base;
            unsigned int lane;

            base =
                index * (uint64_t)modulus_count;

            for (lane = 0; lane < modulus_count; ++lane)
            {
                if (output.values[base + lane] != 0)
                {
                    any_nonzero = 1;
                    break;
                }
            }
        }

        if (!any_nonzero)
        {
            unsigned int i;

            for (i = 0; i < *factor_count; ++i)
            {
                n5_dense_multi_factor_free(&factors[i]);
            }

            *factor_count = 0;

            output.arity = 0;
            output.table_size = 1;

            factors[0] = output;
            *factor_count = 1;

            return 1;
        }
    }

    /* remove incident factors, append output */
    {
        N5DenseMultiFactor survivors[N5_SPARSE_NETWORK_MAX_FACTORS];
        unsigned int survivor_count;
        unsigned int i;

        survivor_count = 0;

        for (i = 0; i < *factor_count; ++i)
        {
            int incident;
            unsigned int j;

            incident = 0;

            for (j = 0; j < incident_count; ++j)
            {
                if (incident_indices[j] == i)
                {
                    incident = 1;
                    break;
                }
            }

            if (incident)
            {
                n5_dense_multi_factor_free(&factors[i]);
            }
            else
            {
                survivors[survivor_count] = factors[i];
                ++survivor_count;
            }
        }

        for (i = 0; i < survivor_count; ++i)
        {
            factors[i] = survivors[i];
        }

        factors[survivor_count] = output;
        *factor_count = survivor_count + 1u;
    }

    return 1;
}

/*
 * multi-modulus dense contraction public API。
 */
int n5_dense_modular_contract_sparse_network_multi(
    const N5SparseNetwork* network,
    const N5SparseContractionPlan* plan,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5DenseMultiModularContractionResult* result
)
{
    N5DenseVariableDomain domains[N5_SPARSE_MAX_VARIABLES];
    N5DenseMultiFactor factors[N5_SPARSE_NETWORK_MAX_FACTORS];

    unsigned int factor_count;
    unsigned int factor_index;
    unsigned int step;
    unsigned int lane;

    uint64_t max_table_entries;

    if (network == NULL || plan == NULL ||
        moduli == NULL || result == NULL ||
        modulus_count == 0 ||
        modulus_count > N5_DENSE_MAX_MODULI ||
        !n5_sparse_plan_is_valid(plan) ||
        !n5_dense_build_domains(network, domains))
    {
        return 0;
    }

    for (lane = 0; lane < modulus_count; ++lane)
    {
        if (moduli[lane] <= 1)
        {
            return 0;
        }
    }

    memset(factors, 0, sizeof(factors));
    memset(result, 0, sizeof(*result));

    result->modulus_count = (uint8_t)modulus_count;

    factor_count = network->factor_count;
    max_table_entries = 1;

    for (factor_index = 0;
         factor_index < factor_count;
         ++factor_index)
    {
        if (!n5_dense_multi_factor_from_sparse(
                &network->factors[factor_index],
                domains,
                moduli,
                modulus_count,
                &factors[factor_index]
            ))
        {
            unsigned int cleanup;

            for (cleanup = 0; cleanup < factor_index; ++cleanup)
            {
                n5_dense_multi_factor_free(&factors[cleanup]);
            }

            return 0;
        }

        if (factors[factor_index].table_size > max_table_entries)
        {
            max_table_entries = factors[factor_index].table_size;
        }
    }

    for (step = 0;
         step < plan->variable_count;
         ++step)
    {
        if (!n5_dense_multi_eliminate_variable(
                factors,
                &factor_count,
                plan->elimination_order[step],
                domains,
                moduli,
                modulus_count,
                &max_table_entries
            ))
        {
            unsigned int cleanup;

            for (cleanup = 0; cleanup < factor_count; ++cleanup)
            {
                n5_dense_multi_factor_free(&factors[cleanup]);
            }

            return 0;
        }
    }

    for (lane = 0; lane < modulus_count; ++lane)
    {
        result->values[lane] = 1u;
    }

    for (factor_index = 0;
         factor_index < factor_count;
         ++factor_index)
    {
        uint64_t base;

        if (factors[factor_index].arity != 0 ||
            factors[factor_index].table_size != 1)
        {
            unsigned int cleanup;

            for (cleanup = 0; cleanup < factor_count; ++cleanup)
            {
                n5_dense_multi_factor_free(&factors[cleanup]);
            }

            return 0;
        }

        base = 0;

        for (lane = 0; lane < modulus_count; ++lane)
        {
            result->values[lane] =
                n5_dense_mod_mul(
                    result->values[lane],
                    factors[factor_index].values[base + lane],
                    moduli[lane]
                );
        }
    }

    for (factor_index = 0;
         factor_index < factor_count;
         ++factor_index)
    {
        n5_dense_multi_factor_free(&factors[factor_index]);
    }

    result->max_table_entries = max_table_entries;
    result->eliminated_variables = plan->variable_count;

    return 1;
}

/* =============================================================== */
/* Generic sliced shared multi-prime traversal (experimental)       */
/* =============================================================== */

typedef struct
{
    const N5SparseNetwork* network;
    const N5SlicePlan* plan;
    const uint32_t* moduli;
    unsigned int modulus_count;
    uint8_t domains[N5_SLICING_MAX_VARIABLES][256];
    unsigned int domain_counts[N5_SLICING_MAX_VARIABLES];
    uint8_t values[N5_SLICING_MAX_VARIABLES];
    uint32_t residues[N5_DENSE_MAX_MODULI];
    uint64_t visited;
    uint64_t zero_assignments;
    uint64_t max_entries;
    int success;
} N5DenseSharedSliceContext;

static void n5_dense_shared_slice_recursive(N5DenseSharedSliceContext* ctx,
                                             unsigned int depth)
{
    if (!ctx->success) return;
    if (depth != ctx->plan->variable_count)
    {
        unsigned int j;
        for (j=0;j<ctx->domain_counts[depth];++j)
        {
            ctx->values[depth]=ctx->domains[depth][j];
            n5_dense_shared_slice_recursive(ctx,depth+1u);
            if (!ctx->success) return;
        }
        return;
    }
    {
        N5SparseNetwork conditioned;
        N5SparseContractionPlan local_plan;
        N5DenseMultiModularContractionResult contracted;
        unsigned int lane;
        ++ctx->visited;
        if (!n5_slicing_condition_sparse_network(
                 ctx->network,ctx->plan->variables,ctx->values,
                 ctx->plan->variable_count,&conditioned))
        {
            ctx->success=0;
            return;
        }
        if (n5_dense_network_is_zero(&conditioned))
        {
            ++ctx->zero_assignments;
            n5_sparse_network_free(&conditioned);
            return;
        }
        if (!n5_sparse_choose_plan(&conditioned,&local_plan) ||
            !n5_dense_modular_contract_sparse_network_multi(
                 &conditioned,&local_plan,ctx->moduli,
                 ctx->modulus_count,&contracted))
        {
            ctx->success=0;
            n5_sparse_network_free(&conditioned);
            return;
        }
        for (lane=0;lane<ctx->modulus_count;++lane)
        {
            ctx->residues[lane]=n5_dense_mod_add(
                ctx->residues[lane],contracted.values[lane],ctx->moduli[lane]);
        }
        if (contracted.max_table_entries>ctx->max_entries)
            ctx->max_entries=contracted.max_table_entries;
        n5_sparse_network_free(&conditioned);
    }
}

int n5_dense_modular_contract_sliced_multi(
    const N5SparseNetwork* network,
    const N5SlicePlan* slice_plan,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5DenseSlicedMultiResult* result
)
{
    N5DenseSharedSliceContext ctx;
    unsigned int i;
    if (!network || !slice_plan || !moduli || !result ||
        modulus_count==0u || modulus_count>N5_DENSE_MAX_MODULI ||
        slice_plan->variable_count>N5_SLICING_MAX_VARIABLES)
        return 0;
    memset(&ctx,0,sizeof(ctx));
    ctx.network=network;ctx.plan=slice_plan;
    ctx.moduli=moduli;ctx.modulus_count=modulus_count;ctx.success=1;
    for (i=0;i<modulus_count;++i)
        if (moduli[i]<=1u) return 0;
    for (i=0;i<slice_plan->variable_count;++i)
    {
        if (!n5_dense_collect_variable_domain(network,
            slice_plan->variables[i],ctx.domains[i],&ctx.domain_counts[i]) ||
            ctx.domain_counts[i]==0u)
            return 0;
    }
    n5_dense_shared_slice_recursive(&ctx,0u);
    if (!ctx.success) return 0;
    memset(result,0,sizeof(*result));
    result->modulus_count=(uint8_t)modulus_count;
    memcpy(result->values,ctx.residues,sizeof(ctx.residues));
    result->raw_assignment_count=ctx.visited;
    result->zero_assignment_count=ctx.zero_assignments;
    result->max_table_entries=ctx.max_entries;
    return 1;
}

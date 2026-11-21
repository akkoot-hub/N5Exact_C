#include "n5_odd_network.h"

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "n5_contraction.h"
#include "n5_kernel.h"

/*
 * Step 9: exact memoization of the even cell T3 tensor.
 * A value depends ONLY on the cardinality k and the three canonical channel
 * indices, not on the current K orbit, fixed r charge, or contraction order.
 * The channel ordering comes from n5_generate_even_fourier_channels, exactly
 * as in the unoptimized construction below.
 *
 * Every worker owns its immutable basis and lazy value table. No shared
 * mutable state or initialization races; at k=4 the largest table has
 * the fourth-cardinality channel cube of signed 32-bit entries per worker.
 */
#ifndef N5_STEP9_DISABLE_EVEN_CELL_CACHE
#if defined(_MSC_VER)
#define N5_STEP9_TLS __declspec(thread)
#else
#define N5_STEP9_TLS _Thread_local
#endif

typedef struct {
    unsigned int initialized;
    unsigned int count;
    N5EvenFourierChannel channels[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];
    uint32_t* values;
} N5Step9EvenTensorCache;

static N5_STEP9_TLS N5Step9EvenTensorCache n5_step9_even_cache[9];
/* The OS releases these thread-local tables on worker exit. No heap allocation
   is required for these per-thread lookup tables. */
static N5_STEP9_TLS uint32_t n5_step9_values_k0[1];
static N5_STEP9_TLS uint32_t n5_step9_values_k2[28u*28u*28u];
static N5_STEP9_TLS uint32_t n5_step9_values_k4[70u*70u*70u];
static N5_STEP9_TLS uint32_t n5_step9_values_k6[28u*28u*28u];
static N5_STEP9_TLS uint32_t n5_step9_values_k8[1];

static int n5_step9_even_cell_value(
    unsigned int k,
    unsigned int qa,
    unsigned int qb,
    unsigned int qc,
    unsigned int expected_count,
    int* value)
{
    N5Step9EvenTensorCache* cache;
    size_t offset;
    uint32_t packed;
    int computed;

    if (k > 8u || (k & 1u) != 0u || value == NULL)
        return 0;
    cache = &n5_step9_even_cache[k];
    if (!cache->initialized) {
        cache->count = n5_generate_even_fourier_channels(
            k, cache->channels, 70u);
        if (cache->count == 0u || cache->count > 70u)
            return 0;
        cache->initialized = 1u;
    }
    if (cache->count != expected_count ||
        qa >= cache->count || qb >= cache->count || qc >= cache->count)
        return 0;

    if (cache->values == NULL) {
        switch (k) {
            case 0u: cache->values = n5_step9_values_k0; break;
            case 2u: cache->values = n5_step9_values_k2; break;
            case 4u: cache->values = n5_step9_values_k4; break;
            case 6u: cache->values = n5_step9_values_k6; break;
            case 8u: cache->values = n5_step9_values_k8; break;
            default: return 0;
        }
    }
    offset = ((size_t)qa * (size_t)cache->count + (size_t)qb)
             * (size_t)cache->count + (size_t)qc;
    /* A stored zero means 'not computed'. Encode the signed exact result
       by a bias of 2^31; INT32_MIN alone would encode as zero and simply
       remains uncached. Conversion uses int64_t, avoiding signed overflow. */
    packed = cache->values[offset];
    if (packed != 0u) {
        *value = (int)((int64_t)packed - INT64_C(2147483648));
        return 1;
    }
    computed = n5_even_cell_tensor_value(
        &cache->channels[qa], &cache->channels[qb], &cache->channels[qc]);
    if (computed != INT32_MIN)
        cache->values[offset] = (uint32_t)((int64_t)computed + INT64_C(2147483648));
    *value = computed;
    return 1;
}
#endif /* N5_STEP9_DISABLE_EVEN_CELL_CACHE */

/* fixed-K cell cardinality helper。 */
static unsigned int n5_odd_network_cell_size(
    const N5SizeMatrix* size_matrix,
    unsigned int cell
)
{
    return size_matrix->matrix[cell >> 2u][cell & 3u];
}

/* unified incidence variable ID。 */
static unsigned int n5_odd_network_var(
    unsigned int cell,
    unsigned int side
)
{
    return 3u * cell + side;
}

/* 判断 cardinality 是否 odd。 */
static int n5_odd_network_is_odd_k(unsigned int k)
{
    return (k & 1u) != 0u;
}

/*
 * 从 KernelLibrary 按 fixed r 构造一个 cell 的 channel domain；
 * 三条 incidence 使用同一个 domain。
 */
static int n5_odd_network_fill_cell_domain(
    const N5SizeMatrix* size_matrix,
    unsigned int cell,
    N5V3 r,
    const N5OddKernelLibrary* library,
    N5OddFixedChargeDomains* domains
)
{
    unsigned int k;
    unsigned int q;
    unsigned int count;
    unsigned int side;

    k = n5_odd_network_cell_size(size_matrix, cell);
    if (k > 8u)
    {
        return 0;
    }

    count = 0u;
    for (q = 0u; q < library->channel_count[k]; ++q)
    {
        if (library->channels[k][q].r == r)
        {
            if (count >= N5_ODD_NETWORK_MAX_DOMAIN)
            {
                return 0;
            }
            for (side = 0u; side < 3u; ++side)
            {
                domains->channel_ids[n5_odd_network_var(cell, side)][count] =
                    (uint8_t)q;
            }
            ++count;
        }
    }

    if (count == 0u)
    {
        return 0;
    }

    for (side = 0u; side < 3u; ++side)
    {
        domains->domain_size[n5_odd_network_var(cell, side)] = (uint8_t)count;
    }

    return 1;
}

/*
 * odd cell 的 transformed T3 exact value。
 * fixed charge 保证三个 legs 的 r 相同；domain construction 已经施加该条件。
 *
 * k=1: 64
 * k=3: 64*M3(a,b,c)
 * k=5: 64*M5(a,b,c)
 * k=7: 64*48
 */
static int n5_odd_network_odd_cell_value(
    unsigned int k,
    const N5OddBasisChannel* qa,
    const N5OddBasisChannel* qb,
    const N5OddBasisChannel* qc,
    const N5OddKernelLibrary* library
)
{
    int multiplicity;

    if (qa == NULL || qb == NULL || qc == NULL || library == NULL)
    {
        return 0;
    }

    if (qa->r != qb->r || qa->r != qc->r)
    {
        return 0;
    }

    if (k == 1u)
    {
        return 64 * (int)library->odd_scalar_k1;
    }
    if (k == 7u)
    {
        return 64 * (int)library->odd_scalar_k7;
    }

    if (qa->fano_label < 1u || qa->fano_label > 7u ||
        qb->fano_label < 1u || qb->fano_label > 7u ||
        qc->fano_label < 1u || qc->fano_label > 7u)
    {
        return 0;
    }

    if (k == 3u)
    {
        multiplicity = n5_odd_kernel_m3(
            library,
            qa->fano_label,
            qb->fano_label,
            qc->fano_label
        );
    }
    else if (k == 5u)
    {
        multiplicity = n5_odd_kernel_m5(
            library,
            qa->fano_label,
            qb->fano_label,
            qc->fano_label
        );
    }
    else
    {
        return 0;
    }

    return 64 * multiplicity;
}

/* 构造一个 cell factor。 */
static int n5_odd_network_build_cell_factor(
    const N5SizeMatrix* size_matrix,
    unsigned int cell,
    const N5OddKernelLibrary* library,
    const N5OddFixedChargeDomains* domains,
    N5SparseFactor* factor
)
{
    unsigned int k;
    unsigned int var_a;
    unsigned int var_b;
    unsigned int var_c;
    uint8_t vars[3];
    unsigned int ia;
    unsigned int ib;
    unsigned int ic;
    unsigned int reserve;

#ifdef N5_STEP9_DISABLE_EVEN_CELL_CACHE
    N5EvenFourierChannel even_channels[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];
#endif
    unsigned int even_count;

    k = n5_odd_network_cell_size(size_matrix, cell);
    var_a = n5_odd_network_var(cell, N5_CONTRACTION_SIDE_A);
    var_b = n5_odd_network_var(cell, N5_CONTRACTION_SIDE_B);
    var_c = n5_odd_network_var(cell, N5_CONTRACTION_SIDE_C);

    vars[0] = (uint8_t)var_a;
    vars[1] = (uint8_t)var_b;
    vars[2] = (uint8_t)var_c;

    reserve =
        (unsigned int)domains->domain_size[var_a] *
        (unsigned int)domains->domain_size[var_b] *
        (unsigned int)domains->domain_size[var_c];

    if (!n5_sparse_factor_init(factor, vars, 3u, reserve))
    {
        return 0;
    }

    even_count = 0u;
    if (!n5_odd_network_is_odd_k(k))
    {
#ifdef N5_STEP9_DISABLE_EVEN_CELL_CACHE
        even_count = n5_generate_even_fourier_channels(k, even_channels, N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY);
#else
        even_count = library->channel_count[k];
#endif
        if (even_count == 0u || even_count > N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY ||
            even_count != library->channel_count[k])
        {
            n5_sparse_factor_free(factor);
            return 0;
        }
    }

    for (ia = 0u; ia < domains->domain_size[var_a]; ++ia)
    {
        unsigned int qa_id = domains->channel_ids[var_a][ia];
        for (ib = 0u; ib < domains->domain_size[var_b]; ++ib)
        {
            unsigned int qb_id = domains->channel_ids[var_b][ib];
            for (ic = 0u; ic < domains->domain_size[var_c]; ++ic)
            {
                unsigned int qc_id = domains->channel_ids[var_c][ic];
                int value;
                uint8_t assignment[3];

                if (n5_odd_network_is_odd_k(k))
                {
                    value = n5_odd_network_odd_cell_value(
                        k,
                        &library->channels[k][qa_id],
                        &library->channels[k][qb_id],
                        &library->channels[k][qc_id],
                        library
                    );
                }
                else
                {
                    if (qa_id >= even_count || qb_id >= even_count || qc_id >= even_count)
                    {
                        n5_sparse_factor_free(factor);
                        return 0;
                    }
#ifdef N5_STEP9_DISABLE_EVEN_CELL_CACHE
                    value = n5_even_cell_tensor_value(
                        &even_channels[qa_id],
                        &even_channels[qb_id],
                        &even_channels[qc_id]
                    );
#else
                    if (!n5_step9_even_cell_value(k, qa_id, qb_id, qc_id,
                                                  even_count, &value)) {
                        n5_sparse_factor_free(factor);
                        return 0;
                    }
#endif
                }

                if (value == 0)
                {
                    continue;
                }

                assignment[0] = (uint8_t)qa_id;
                assignment[1] = (uint8_t)qb_id;
                assignment[2] = (uint8_t)qc_id;

                {
                    int64_t denominator;
                    denominator =
                        (int64_t)library->channels[k][qa_id].norm *
                        (int64_t)library->channels[k][qb_id].norm *
                        (int64_t)library->channels[k][qc_id].norm;

                    if (!n5_sparse_factor_append(
                            factor,
                            assignment,
                            n5_small_rational_make(value, denominator)
                        ))
                {
                        n5_sparse_factor_free(factor);
                        return 0;
                    }
                }
            }
        }
    }

    return n5_sparse_factor_compress(factor);
}

/* 构造一条 ordered partition line factor。 */
static int n5_odd_network_build_line_factor(
    const N5SizeMatrix* size_matrix,
    unsigned int line_id,
    const N5OddKernelLibrary* library,
    const N5OddFixedChargeDomains* domains,
    N5SparseFactor* factor,
    unsigned int* odd_degree_out
)
{
    unsigned int cells[4];
    unsigned int side;
    uint8_t vars[4];
    unsigned int sizes[4];
    unsigned int leg;
    unsigned int reserve;
    unsigned int i0, i1, i2, i3;
    unsigned int odd_degree;

    if (!n5_contraction_line_cells(line_id, cells))
    {
        return 0;
    }

    side = n5_contraction_line_side(line_id);
    if (side > N5_CONTRACTION_SIDE_C)
    {
        return 0;
    }

    reserve = 1u;
    odd_degree = 0u;
    for (leg = 0u; leg < 4u; ++leg)
    {
        unsigned int var;
        sizes[leg] = n5_odd_network_cell_size(size_matrix, cells[leg]);
        if (sizes[leg] & 1u)
        {
            ++odd_degree;
        }
        var = n5_odd_network_var(cells[leg], side);
        vars[leg] = (uint8_t)var;
        reserve *= domains->domain_size[var];
    }

    if (odd_degree != 0u && odd_degree != 2u && odd_degree != 4u)
    {
        return 0;
    }

    if (!n5_sparse_factor_init(factor, vars, 4u, reserve))
    {
        return 0;
    }

    for (i0 = 0u; i0 < domains->domain_size[vars[0]]; ++i0)
    {
        unsigned int q0 = domains->channel_ids[vars[0]][i0];
        for (i1 = 0u; i1 < domains->domain_size[vars[1]]; ++i1)
        {
            unsigned int q1 = domains->channel_ids[vars[1]][i1];
            for (i2 = 0u; i2 < domains->domain_size[vars[2]]; ++i2)
            {
                unsigned int q2 = domains->channel_ids[vars[2]][i2];
                for (i3 = 0u; i3 < domains->domain_size[vars[3]]; ++i3)
                {
                    unsigned int q3 = domains->channel_ids[vars[3]][i3];
                    int value;
                    uint8_t assignment[4];

                    value = n5_odd_kernel_line4(
                        &library->channels[sizes[0]][q0],
                        &library->channels[sizes[1]][q1],
                        &library->channels[sizes[2]][q2],
                        &library->channels[sizes[3]][q3]
                    );

                    if (value == 0)
                    {
                        continue;
                    }

                    assignment[0] = (uint8_t)q0;
                    assignment[1] = (uint8_t)q1;
                    assignment[2] = (uint8_t)q2;
                    assignment[3] = (uint8_t)q3;

                    if (!n5_sparse_factor_append(
                            factor,
                            assignment,
                            n5_small_rational_make(value, 1)
                        ))
                    {
                        n5_sparse_factor_free(factor);
                        return 0;
                    }
                }
            }
        }
    }

    if (!n5_sparse_factor_compress(factor))
    {
        n5_sparse_factor_free(factor);
        return 0;
    }

    if (odd_degree_out != NULL)
    {
        *odd_degree_out = odd_degree;
    }

    return 1;
}

void n5_odd_network_init(N5OddRepresentationNetwork* network)
{
    if (network == NULL)
    {
        return;
    }

    memset(network, 0, sizeof(*network));
    n5_sparse_network_init(&network->sparse);
}

void n5_odd_network_free(N5OddRepresentationNetwork* network)
{
    if (network == NULL)
    {
        return;
    }

    n5_sparse_network_free(&network->sparse);
    memset(network, 0, sizeof(*network));
}

int n5_odd_network_build(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[N5_ODD_NETWORK_CELL_COUNT],
    const N5OddKernelLibrary* library,
    N5OddRepresentationNetwork* network
)
{
    unsigned int cell;
    unsigned int line_id;

    if (size_matrix == NULL || r_by_cell == NULL || library == NULL || network == NULL)
    {
        return 0;
    }

    n5_odd_network_init(network);

    for (cell = 0u; cell < 16u; ++cell)
    {
        unsigned int k = n5_odd_network_cell_size(size_matrix, cell);

        if ((k == 0u || k == 8u) && r_by_cell[cell] != 0u)
        {
            n5_odd_network_free(network);
            return 0;
        }

        if (!n5_odd_network_fill_cell_domain(
                size_matrix,
                cell,
                r_by_cell[cell],
                library,
                &network->domains
            ))
        {
            n5_odd_network_free(network);
            return 0;
        }

        if (k & 1u)
        {
            ++network->odd_cell_count;
        }
        else
        {
            ++network->even_cell_count;
        }
    }

    /* 16 cell factors */
    for (cell = 0u; cell < 16u; ++cell)
    {
        N5SparseFactor factor;
        memset(&factor, 0, sizeof(factor));

        if (!n5_odd_network_build_cell_factor(
                size_matrix,
                cell,
                library,
                &network->domains,
                &factor
            ) ||
            !n5_sparse_network_add_move(&network->sparse, &factor))
        {
            n5_sparse_factor_free(&factor);
            n5_odd_network_free(network);
            return 0;
        }
    }

    /* 12 line factors */
    for (line_id = 0u; line_id < 12u; ++line_id)
    {
        N5SparseFactor factor;
        unsigned int odd_degree;
        memset(&factor, 0, sizeof(factor));

        if (!n5_odd_network_build_line_factor(
                size_matrix,
                line_id,
                library,
                &network->domains,
                &factor,
                &odd_degree
            ) ||
            !n5_sparse_network_add_move(&network->sparse, &factor))
        {
            n5_sparse_factor_free(&factor);
            n5_odd_network_free(network);
            return 0;
        }

        if (odd_degree == 0u)
        {
            ++network->line_odd_degree_0;
        }
        else if (odd_degree == 2u)
        {
            ++network->line_odd_degree_2;
        }
        else
        {
            ++network->line_odd_degree_4;
        }
    }

    /* inverse-Gram metrics have been absorbed into each cell factor. */

    if (network->sparse.factor_count != 28u)
    {
        n5_odd_network_free(network);
        return 0;
    }

    network->primitive_nnz = n5_sparse_network_primitive_nnz(&network->sparse);
    return 1;
}


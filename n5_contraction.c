#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "n5_contraction.h"
#include "n5_sparse.h"

#include "n5_kernel.h"
#include "n5_orbit.h"

/*
 * 将一个 cell-side incidence 编码为：
 *
 *   var = 3*cell + side。
 */
unsigned int n5_contraction_incidence_var(
    unsigned int cell,
    unsigned int side
)
{
    if (cell >= N5_CONTRACTION_CELL_COUNT)
    {
        return N5_CONTRACTION_INVALID_ID;
    }

    if (side > N5_CONTRACTION_SIDE_C)
    {
        return N5_CONTRACTION_INVALID_ID;
    }

    return
        3u * cell + side;
}

/*
 * 返回 line 所属 side。
 *
 *   0..3  -> A
 *   4..7  -> B
 *   8..11 -> C
 */
unsigned int n5_contraction_line_side(
    unsigned int line_id
)
{
    if (line_id >= N5_CONTRACTION_LINE_COUNT)
    {
        return N5_CONTRACTION_INVALID_ID;
    }

    return line_id >> 2u;
}

/*
 * 生成一条 top-level line 上的四个 cell IDs。
 *
 * 顺序固定为：
 *
 *   A-row：
 *     beta = 0,1,2,3；
 *
 *   B-column：
 *     alpha = 0,1,2,3；
 *
 *   C-xor diagonal：
 *     alpha = 0,1,2,3，
 *     beta = alpha XOR gamma。
 */
int n5_contraction_line_cells(
    unsigned int line_id,
    unsigned int cells[4]
)
{
    unsigned int index;

    if (cells == NULL ||
        line_id >= N5_CONTRACTION_LINE_COUNT)
    {
        return 0;
    }

    if (line_id < 4)
    {
        unsigned int alpha;

        alpha = line_id;

        for (index = 0; index < 4; ++index)
        {
            cells[index] =
                4u * alpha + index;
        }

        return 1;
    }

    if (line_id < 8)
    {
        unsigned int beta;

        beta = line_id - 4u;

        for (index = 0; index < 4; ++index)
        {
            cells[index] =
                4u * index + beta;
        }

        return 1;
    }

    {
        unsigned int gamma;

        gamma = line_id - 8u;

        for (index = 0; index < 4; ++index)
        {
            unsigned int alpha;
            unsigned int beta;

            alpha = index;
            beta = alpha ^ gamma;

            cells[index] =
                4u * alpha + beta;
        }
    }

    return 1;
}

/*
 * 判断固定 K 是否为合法 all-even size matrix。
 */
int n5_contraction_k_is_all_even_legal(
    const N5SizeMatrix* size_matrix
)
{
    unsigned int row;
    unsigned int column;

    if (size_matrix == NULL)
    {
        return 0;
    }

    if (!n5_size_matrix_is_legal(size_matrix))
    {
        return 0;
    }

    for (row = 0;
         row < N5_V2_SIZE;
         ++row)
    {
        for (column = 0;
             column < N5_V2_SIZE;
             ++column)
        {
            if ((size_matrix->matrix[row][column] & 1u) != 0)
            {
                return 0;
            }
        }
    }

    return 1;
}

/*
 * 构造固定 all-even K 的 scalar GF(2) r-flow basis。
 *
 * 12 条 line equations 先表示成 16-bit row masks，
 * 然后执行 GF(2) reduced row echelon elimination。
 */
int n5_contraction_build_scalar_flow_basis(
    const N5SizeMatrix* size_matrix,
    N5ScalarFlowBasis* basis
)
{
    uint16_t active_mask;

    uint16_t rows[N5_CONTRACTION_LINE_COUNT];
    unsigned int row_count;

    unsigned int pivot_columns[N5_CONTRACTION_CELL_COUNT];
    unsigned int rank;

    unsigned int cell;
    unsigned int line_id;

    if (basis == NULL ||
        !n5_contraction_k_is_all_even_legal(
            size_matrix
        ))
    {
        return 0;
    }

    basis->active_count = 0;
    basis->rank = 0;
    basis->nullity = 0;

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        basis->active_cells[cell] = 0;
        basis->basis_masks[cell] = 0;
    }

    active_mask = 0;

    /*
     * k=0 和 k=8 只有 r=0 channel。
     *
     * 其余 even cardinalities 可携带非零 r，
     * 因此作为 active scalar variables。
     */
    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int row;
        unsigned int column;
        unsigned int k;

        row =
            cell >> 2u;

        column =
            cell & 3u;

        k =
            size_matrix->matrix[row][column];

        if (k != 0 &&
            k != 8)
        {
            basis->active_cells[
                basis->active_count
            ] =
                (uint8_t)cell;

            ++basis->active_count;

            active_mask |=
                (uint16_t)(1u << cell);
        }
    }

    /*
     * 建立 12 条 line constraint rows。
     *
     * inactive cells 被自动移除。
     */
    row_count = 0;

    for (line_id = 0;
         line_id < N5_CONTRACTION_LINE_COUNT;
         ++line_id)
    {
        unsigned int cells[4];
        unsigned int index;

        uint16_t row_mask;

        if (!n5_contraction_line_cells(
                line_id,
                cells
            ))
        {
            return 0;
        }

        row_mask = 0;

        for (index = 0; index < 4; ++index)
        {
            row_mask |=
                (uint16_t)(
                    1u << cells[index]
                );
        }

        row_mask &=
            active_mask;

        if (row_mask != 0)
        {
            rows[row_count] =
                row_mask;

            ++row_count;
        }
    }

    /*
     * GF(2) RREF。
     */
    rank = 0;

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int pivot_row;
        unsigned int search_row;

        if ((active_mask &
             (uint16_t)(1u << cell)) == 0)
        {
            continue;
        }

        pivot_row =
            row_count;

        for (search_row = rank;
             search_row < row_count;
             ++search_row)
        {
            if ((rows[search_row] &
                 (uint16_t)(1u << cell)) != 0)
            {
                pivot_row =
                    search_row;

                break;
            }
        }

        if (pivot_row == row_count)
        {
            continue;
        }

        {
            uint16_t temporary;

            temporary =
                rows[rank];

            rows[rank] =
                rows[pivot_row];

            rows[pivot_row] =
                temporary;
        }

        /*
         * 消掉其他所有 rows 中的当前 pivot bit，
         * 直接得到 reduced row echelon form。
         */
        for (search_row = 0;
             search_row < row_count;
             ++search_row)
        {
            if (search_row == rank)
            {
                continue;
            }

            if ((rows[search_row] &
                 (uint16_t)(1u << cell)) != 0)
            {
                rows[search_row] ^=
                    rows[rank];
            }
        }

        pivot_columns[rank] =
            cell;

        ++rank;

        if (rank == row_count)
        {
            break;
        }
    }

    basis->rank =
        (uint8_t)rank;

    basis->nullity =
        (uint8_t)(
            basis->active_count - rank
        );

    /*
     * 为每个 free active column 构造一个
     * global 16-bit nullspace basis vector。
     */
    {
        unsigned int basis_index;

        basis_index = 0;

        for (cell = 0;
             cell < N5_CONTRACTION_CELL_COUNT;
             ++cell)
        {
            int is_active;
            int is_pivot;

            unsigned int pivot_index;

            uint16_t vector;

            is_active =
                (active_mask &
                 (uint16_t)(1u << cell)) != 0;

            if (!is_active)
            {
                continue;
            }

            is_pivot = 0;

            for (pivot_index = 0;
                 pivot_index < rank;
                 ++pivot_index)
            {
                if (pivot_columns[pivot_index] ==
                    cell)
                {
                    is_pivot = 1;
                    break;
                }
            }

            if (is_pivot)
            {
                continue;
            }

            vector =
                (uint16_t)(1u << cell);

            /*
             * RREF 中：
             *
             *   x_pivot
             *   =
             *   sum(coeff * x_free)
             *
             * over GF(2)。
             *
             * 当前只有一个 free variable = 1。
             */
            for (pivot_index = 0;
                 pivot_index < rank;
                 ++pivot_index)
            {
                if ((rows[pivot_index] &
                     (uint16_t)(1u << cell)) != 0)
                {
                    vector |=
                        (uint16_t)(
                            1u <<
                            pivot_columns[pivot_index]
                        );
                }
            }

            basis->basis_masks[basis_index] =
                vector;

            ++basis_index;
        }

        if (basis_index !=
            basis->nullity)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 根据 scalar nullity 计算 V3-valued flow 数量：
 *
 *   8^nullity。
 */
uint64_t n5_contraction_count_v3_flows(
    const N5ScalarFlowBasis* basis
)
{
    uint64_t count;
    unsigned int index;

    if (basis == NULL)
    {
        return 0;
    }

    count = 1;

    for (index = 0;
         index < basis->nullity;
         ++index)
    {
        count <<= 3u;
    }

    return count;
}

/*
 * 检查一个 V3-valued cell flow。
 */
int n5_contraction_v3_flow_is_valid(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT]
)
{
    unsigned int cell;
    unsigned int line_id;

    if (!n5_contraction_k_is_all_even_legal(
            size_matrix
        ) ||
        flow == NULL)
    {
        return 0;
    }

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int row;
        unsigned int column;
        unsigned int k;

        if (!n5_v3_is_valid(flow[cell]))
        {
            return 0;
        }

        row =
            cell >> 2u;

        column =
            cell & 3u;

        k =
            size_matrix->matrix[row][column];

        if ((k == 0 || k == 8) &&
            flow[cell] != 0)
        {
            return 0;
        }
    }

    for (line_id = 0;
         line_id < N5_CONTRACTION_LINE_COUNT;
         ++line_id)
    {
        unsigned int cells[4];
        unsigned int index;

        N5V3 total;

        if (!n5_contraction_line_cells(
                line_id,
                cells
            ))
        {
            return 0;
        }

        total = 0;

        for (index = 0; index < 4; ++index)
        {
            total ^=
                flow[cells[index]];
        }

        if (total != 0)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 构造 48 个 fixed-flow restricted channel domains。
 */
int n5_contraction_build_fixed_flow_domains(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5FixedFlowDomains* domains
)
{
    unsigned int cell;

    if (domains == NULL ||
        !n5_contraction_v3_flow_is_valid(
            size_matrix,
            flow
        ))
    {
        return 0;
    }

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        N5EvenFourierChannel channels[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];

        unsigned int row;
        unsigned int column;
        unsigned int k;

        unsigned int channel_count;
        unsigned int channel_index;

        uint8_t surviving_ids[
            N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN
        ];

        unsigned int surviving_count;

        unsigned int side;

        row =
            cell >> 2u;

        column =
            cell & 3u;

        k =
            size_matrix->matrix[row][column];

        channel_count =
            n5_generate_even_fourier_channels(
                k,
                channels,
                N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
            );

        if (channel_count == 0)
        {
            return 0;
        }

        surviving_count = 0;

        for (channel_index = 0;
             channel_index < channel_count;
             ++channel_index)
        {
            if (channels[channel_index].r ==
                flow[cell])
            {
                if (surviving_count >=
                    N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN)
                {
                    return 0;
                }

                surviving_ids[surviving_count] =
                    (uint8_t)channel_index;

                ++surviving_count;
            }
        }

        if (surviving_count == 0)
        {
            return 0;
        }

        domains->cell_r[cell] =
            flow[cell];

        for (side = 0;
             side <= N5_CONTRACTION_SIDE_C;
             ++side)
        {
            unsigned int var;
            unsigned int slot;

            var =
                n5_contraction_incidence_var(
                    cell,
                    side
                );

            if (var >=
                N5_CONTRACTION_INCIDENCE_COUNT)
            {
                return 0;
            }

            domains->domain_size[var] =
                (uint8_t)surviving_count;

            for (slot = 0;
                 slot < surviving_count;
                 ++slot)
            {
                domains->channel_ids[var][slot] =
                    surviving_ids[slot];
            }

            for (slot = surviving_count;
                 slot <
                    N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN;
                 ++slot)
            {
                domains->channel_ids[var][slot] = 0;
            }
        }
    }

    return 1;
}

/* ================================================================
 * Stage 17: concrete sparse factors + exact k=2 second-Walsh rewrite
 * ================================================================ */

/*
 * fixed-K cell cardinality helper。
 */
static unsigned int n5_contraction_cell_size(
    const N5SizeMatrix* size_matrix,
    unsigned int cell
)
{
    unsigned int row;
    unsigned int column;

    row = cell >> 2u;
    column = cell & 3u;

    return size_matrix->matrix[row][column];
}

/*
 * 构造 fixed-flow restricted transformed cell factor。
 */
static int n5_contraction_build_cell_sparse_factor(
    const N5SizeMatrix* size_matrix,
    unsigned int cell,
    const N5FixedFlowDomains* domains,
    N5SparseFactor* factor
)
{
    N5EvenFourierChannel channels[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];

    unsigned int k;
    unsigned int channel_count;

    unsigned int var_a;
    unsigned int var_b;
    unsigned int var_c;

    uint8_t vars[3];

    unsigned int ia;
    unsigned int ib;
    unsigned int ic;

    unsigned int reserve;

    if (size_matrix == NULL ||
        domains == NULL ||
        factor == NULL ||
        cell >= N5_CONTRACTION_CELL_COUNT)
    {
        return 0;
    }

    k =
        n5_contraction_cell_size(
            size_matrix,
            cell
        );

    channel_count =
        n5_generate_even_fourier_channels(
            k,
            channels,
            N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
        );

    if (channel_count == 0)
    {
        return 0;
    }

    var_a =
        n5_contraction_incidence_var(
            cell,
            N5_CONTRACTION_SIDE_A
        );

    var_b =
        n5_contraction_incidence_var(
            cell,
            N5_CONTRACTION_SIDE_B
        );

    var_c =
        n5_contraction_incidence_var(
            cell,
            N5_CONTRACTION_SIDE_C
        );

    vars[0] = (uint8_t)var_a;
    vars[1] = (uint8_t)var_b;
    vars[2] = (uint8_t)var_c;

    reserve =
        (unsigned int)domains->domain_size[var_a] *
        (unsigned int)domains->domain_size[var_b] *
        (unsigned int)domains->domain_size[var_c];

    if (!n5_sparse_factor_init(
            factor,
            vars,
            3,
            reserve
        ))
    {
        return 0;
    }

    for (ia = 0;
         ia < domains->domain_size[var_a];
         ++ia)
    {
        unsigned int qa;

        qa =
            domains->channel_ids[var_a][ia];

        if (qa >= channel_count)
        {
            n5_sparse_factor_free(factor);
            return 0;
        }

        for (ib = 0;
             ib < domains->domain_size[var_b];
             ++ib)
        {
            unsigned int qb;

            qb =
                domains->channel_ids[var_b][ib];

            if (qb >= channel_count)
            {
                n5_sparse_factor_free(factor);
                return 0;
            }

            for (ic = 0;
                 ic < domains->domain_size[var_c];
                 ++ic)
            {
                unsigned int qc;
                int value;

                uint8_t assignment[3];

                qc =
                    domains->channel_ids[var_c][ic];

                if (qc >= channel_count)
                {
                    n5_sparse_factor_free(factor);
                    return 0;
                }

                value =
                    n5_even_cell_tensor_value(
                        &channels[qa],
                        &channels[qb],
                        &channels[qc]
                    );

                if (value == 0)
                {
                    continue;
                }

                assignment[0] = (uint8_t)qa;
                assignment[1] = (uint8_t)qb;
                assignment[2] = (uint8_t)qc;

                if (!n5_sparse_factor_append(
                        factor,
                        assignment,
                        n5_small_rational_make(
                            value,
                            1
                        )
                    ))
                {
                    n5_sparse_factor_free(factor);
                    return 0;
                }
            }
        }
    }

    return 1;
}

/*
 * 构造 fixed-flow restricted ordered-partition line factor。
 *
 * ordered sizes 使用 line_cells() 的 canonical leg 顺序，
 * 因而不会把 motif 当成无序对象来生成 tensor。
 */
static int n5_contraction_build_line_sparse_factor(
    const N5SizeMatrix* size_matrix,
    unsigned int line_id,
    const N5FixedFlowDomains* domains,
    N5SparseFactor* factor
)
{
    unsigned int cells[4];
    unsigned int sizes[4];
    unsigned int side;

    N5EvenPartitionTensor tensor;

    uint8_t vars[4];

    unsigned int i0;
    unsigned int i1;
    unsigned int i2;
    unsigned int i3;

    unsigned int reserve;
    unsigned int leg;

    if (size_matrix == NULL ||
        domains == NULL ||
        factor == NULL ||
        line_id >= N5_CONTRACTION_LINE_COUNT)
    {
        return 0;
    }

    if (!n5_contraction_line_cells(
            line_id,
            cells
        ))
    {
        return 0;
    }

    side =
        n5_contraction_line_side(
            line_id
        );

    for (leg = 0; leg < 4; ++leg)
    {
        unsigned int var;

        sizes[leg] =
            n5_contraction_cell_size(
                size_matrix,
                cells[leg]
            );

        var =
            n5_contraction_incidence_var(
                cells[leg],
                side
            );

        vars[leg] = (uint8_t)var;
    }

    if (!n5_build_even_partition_tensor(
            sizes,
            &tensor
        ))
    {
        return 0;
    }

    reserve = 1;

    for (leg = 0; leg < 4; ++leg)
    {
        reserve *=
            domains->domain_size[vars[leg]];
    }

    if (!n5_sparse_factor_init(
            factor,
            vars,
            4,
            reserve
        ))
    {
        n5_free_even_partition_tensor(
            &tensor
        );

        return 0;
    }

    for (i0 = 0;
         i0 < domains->domain_size[vars[0]];
         ++i0)
    {
        unsigned int q0;

        q0 =
            domains->channel_ids[vars[0]][i0];

        for (i1 = 0;
             i1 < domains->domain_size[vars[1]];
             ++i1)
        {
            unsigned int q1;

            q1 =
                domains->channel_ids[vars[1]][i1];

            for (i2 = 0;
                 i2 < domains->domain_size[vars[2]];
                 ++i2)
            {
                unsigned int q2;

                q2 =
                    domains->channel_ids[vars[2]][i2];

                for (i3 = 0;
                     i3 < domains->domain_size[vars[3]];
                     ++i3)
                {
                    unsigned int q3;
                    int value;

                    uint8_t assignment[4];

                    q3 =
                        domains->channel_ids[vars[3]][i3];

                    value =
                        n5_even_partition_tensor_value(
                            &tensor,
                            q0,
                            q1,
                            q2,
                            q3
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
                            n5_small_rational_make(
                                value,
                                1
                            )
                        ))
                    {
                        n5_sparse_factor_free(factor);

                        n5_free_even_partition_tensor(
                            &tensor
                        );

                        return 0;
                    }
                }
            }
        }
    }

    n5_free_even_partition_tensor(
        &tensor
    );

    return 1;
}

/*
 * 构造 fixed-flow inverse-Gram unary metric factor。
 *
 * metric(q)=1/norm(q)。
 */
static int n5_contraction_build_metric_sparse_factor(
    const N5SizeMatrix* size_matrix,
    unsigned int cell,
    unsigned int side,
    const N5FixedFlowDomains* domains,
    N5SparseFactor* factor
)
{
    N5EvenFourierChannel channels[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];

    unsigned int k;
    unsigned int channel_count;
    unsigned int var;
    unsigned int index;

    uint8_t vars[1];

    if (size_matrix == NULL ||
        domains == NULL ||
        factor == NULL ||
        cell >= N5_CONTRACTION_CELL_COUNT ||
        side > N5_CONTRACTION_SIDE_C)
    {
        return 0;
    }

    k =
        n5_contraction_cell_size(
            size_matrix,
            cell
        );

    channel_count =
        n5_generate_even_fourier_channels(
            k,
            channels,
            N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
        );

    if (channel_count == 0)
    {
        return 0;
    }

    var =
        n5_contraction_incidence_var(
            cell,
            side
        );

    vars[0] = (uint8_t)var;

    if (!n5_sparse_factor_init(
            factor,
            vars,
            1,
            domains->domain_size[var]
        ))
    {
        return 0;
    }

    for (index = 0;
         index < domains->domain_size[var];
         ++index)
    {
        unsigned int q;
        uint8_t assignment[1];

        q =
            domains->channel_ids[var][index];

        if (q >= channel_count ||
            channels[q].norm == 0)
        {
            n5_sparse_factor_free(factor);
            return 0;
        }

        assignment[0] = (uint8_t)q;

        if (!n5_sparse_factor_append(
                factor,
                assignment,
                n5_small_rational_make(
                    1,
                    channels[q].norm
                )
            ))
        {
            n5_sparse_factor_free(factor);
            return 0;
        }
    }

    return 1;
}

/*
 * 构造完整 fixed-flow raw sparse network。
 */
int n5_contraction_build_fixed_flow_sparse_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5SparseNetwork* network
)
{
    N5FixedFlowDomains domains;

    unsigned int cell;
    unsigned int line_id;

    if (network == NULL ||
        !n5_contraction_build_fixed_flow_domains(
            size_matrix,
            flow,
            &domains
        ))
    {
        return 0;
    }

    n5_sparse_network_init(network);

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        N5SparseFactor factor;

        if (!n5_contraction_build_cell_sparse_factor(
                size_matrix,
                cell,
                &domains,
                &factor
            ) ||
            !n5_sparse_network_add_move(
                network,
                &factor
            ))
        {
            n5_sparse_factor_free(&factor);
            n5_sparse_network_free(network);
            return 0;
        }
    }

    for (line_id = 0;
         line_id < N5_CONTRACTION_LINE_COUNT;
         ++line_id)
    {
        N5SparseFactor factor;

        if (!n5_contraction_build_line_sparse_factor(
                size_matrix,
                line_id,
                &domains,
                &factor
            ) ||
            !n5_sparse_network_add_move(
                network,
                &factor
            ))
        {
            n5_sparse_factor_free(&factor);
            n5_sparse_network_free(network);
            return 0;
        }
    }

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int side;

        for (side = 0;
             side <= N5_CONTRACTION_SIDE_C;
             ++side)
        {
            N5SparseFactor factor;

            if (!n5_contraction_build_metric_sparse_factor(
                    size_matrix,
                    cell,
                    side,
                    &domains,
                    &factor
                ) ||
                !n5_sparse_network_add_move(
                    network,
                    &factor
                ))
            {
                n5_sparse_factor_free(&factor);
                n5_sparse_network_free(network);
                return 0;
            }
        }
    }

    if (network->factor_count !=
        N5_CONTRACTION_TOTAL_FACTOR_COUNT)
    {
        n5_sparse_network_free(network);
        return 0;
    }

    return 1;
}

/*
 * internal k=2 second-Walsh decomposition。
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
} N5K2WalshDecomposition;

/*
 * 找 local decomposition 中 q 的位置。
 */
static int n5_contraction_k2_walsh_q_position(
    const N5K2WalshDecomposition* decomposition,
    unsigned int q
)
{
    unsigned int index;

    for (index = 0;
         index < decomposition->channel_count;
         ++index)
    {
        if (decomposition->channel_ids[index] == q)
        {
            return (int)index;
        }
    }

    return -1;
}

/*
 * definition-derived leg value：
 *
 *   gauge[q] * (-1)^(lambda·s_q)。
 */
static int n5_contraction_k2_walsh_leg_value(
    const N5K2WalshDecomposition* decomposition,
    const N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT],
    unsigned int q,
    N5V3 lambda
)
{
    int position;

    position =
        n5_contraction_k2_walsh_q_position(
            decomposition,
            q
        );

    if (position < 0)
    {
        return 0;
    }

    return
        (int)decomposition->gauge[position] *
        n5_walsh_character(
            lambda,
            channels[q].s
        );
}

/*
 * 从 transformed k=2 cell tensor 自动恢复固定 r 的
 * one-leg canonical gauge signs。
 */
static int n5_contraction_build_k2_walsh_decomposition(
    N5V3 r,
    N5K2WalshDecomposition* decomposition
)
{
    N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT];

    unsigned int channel_count;
    unsigned int q;

    unsigned int candidate_mask;
    unsigned int candidate_count;

    if (decomposition == NULL ||
        !n5_v3_is_valid(r))
    {
        return 0;
    }

    channel_count =
        n5_generate_even_fourier_channels(
            2,
            channels,
            N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT
        );

    if (channel_count != N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT)
    {
        return 0;
    }

    decomposition->r = r;
    decomposition->channel_count = 0;

    for (q = 0; q < channel_count; ++q)
    {
        if (channels[q].r == r)
        {
            unsigned int slot;

            slot =
                decomposition->channel_count;

            if (slot >= 7)
            {
                return 0;
            }

            decomposition->channel_ids[slot] =
                (uint8_t)q;

            ++decomposition->channel_count;
        }
    }

    if (decomposition->channel_count == 0)
    {
        return 0;
    }

    /*
     * Fix first gauge=+1，枚举其余 ±1。
     */
    candidate_count =
        1u <<
        (decomposition->channel_count - 1u);

    for (candidate_mask = 0;
         candidate_mask < candidate_count;
         ++candidate_mask)
    {
        unsigned int index;
        int ok;
        unsigned int nonzero_checked;

        decomposition->gauge[0] = 1;

        for (index = 1;
             index < decomposition->channel_count;
             ++index)
        {
            decomposition->gauge[index] =
                ((candidate_mask >> (index - 1u)) & 1u)
                    ? -1
                    : 1;
        }

        ok = 1;
        nonzero_checked = 0;

        {
            unsigned int ia;
            unsigned int ib;
            unsigned int ic;

            for (ia = 0;
                 ia < decomposition->channel_count && ok;
                 ++ia)
            {
                for (ib = 0;
                     ib < decomposition->channel_count && ok;
                     ++ib)
                {
                    for (ic = 0;
                         ic < decomposition->channel_count;
                         ++ic)
                    {
                        unsigned int qa;
                        unsigned int qb;
                        unsigned int qc;

                        int actual;
                        int actual_sign;
                        int predicted_sign;

                        qa = decomposition->channel_ids[ia];
                        qb = decomposition->channel_ids[ib];
                        qc = decomposition->channel_ids[ic];

                        actual =
                            n5_even_cell_tensor_value(
                                &channels[qa],
                                &channels[qb],
                                &channels[qc]
                            );

                        if (actual == 0)
                        {
                            continue;
                        }

                        ++nonzero_checked;

                        if (actual != 32 &&
                            actual != -32)
                        {
                            return 0;
                        }

                        actual_sign =
                            actual > 0 ? 1 : -1;

                        predicted_sign =
                            decomposition->gauge[ia] *
                            decomposition->gauge[ib] *
                            decomposition->gauge[ic];

                        if (actual_sign !=
                            predicted_sign)
                        {
                            ok = 0;
                            break;
                        }
                    }
                }
            }
        }

        if (ok && nonzero_checked > 0)
        {
            break;
        }
    }

    if (candidate_mask == candidate_count)
    {
        return 0;
    }

    /*
     * Canonical representatives of V3/<r>。
     */
    decomposition->lambda_count = 0;

    if (r == 0)
    {
        unsigned int lambda;

        for (lambda = 0;
             lambda < N5_V3_SIZE;
             ++lambda)
        {
            decomposition->lambda_values[
                decomposition->lambda_count
            ] =
                (N5V3)lambda;

            ++decomposition->lambda_count;
        }

        decomposition->scale = 4;
    }
    else
    {
        unsigned int seen[N5_V3_SIZE];
        unsigned int lambda;

        for (lambda = 0;
             lambda < N5_V3_SIZE;
             ++lambda)
        {
            seen[lambda] = 0;
        }

        for (lambda = 0;
             lambda < N5_V3_SIZE;
             ++lambda)
        {
            unsigned int partner;
            unsigned int representative;

            if (seen[lambda] != 0)
            {
                continue;
            }

            partner =
                lambda ^ r;

            representative =
                lambda < partner
                    ? lambda
                    : partner;

            decomposition->lambda_values[
                decomposition->lambda_count
            ] =
                (N5V3)representative;

            ++decomposition->lambda_count;

            seen[lambda] = 1;
            seen[partner] = 1;
        }

        if (decomposition->lambda_count != 4)
        {
            return 0;
        }

        decomposition->scale = 8;
    }

    return 1;
}

/*
 * 构造一个 k=2 Walsh incidence factor：
 *
 *   L(q,lambda)/norm(q)。
 */
static int n5_contraction_build_k2_walsh_incidence_factor(
    unsigned int cell,
    unsigned int side,
    const N5K2WalshDecomposition* decomposition,
    N5SparseFactor* factor
)
{
    N5EvenFourierChannel channels[N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT];
    unsigned int channel_count;

    uint8_t vars[2];

    unsigned int q_index;
    unsigned int lambda_index;

    if (factor == NULL ||
        decomposition == NULL ||
        cell >= N5_CONTRACTION_CELL_COUNT ||
        side > N5_CONTRACTION_SIDE_C)
    {
        return 0;
    }

    channel_count =
        n5_generate_even_fourier_channels(
            2,
            channels,
            N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT
        );

    if (channel_count != N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT)
    {
        return 0;
    }

    vars[0] =
        (uint8_t)n5_contraction_incidence_var(
            cell,
            side
        );

    vars[1] =
        (uint8_t)(48u + cell);

    if (!n5_sparse_factor_init(
            factor,
            vars,
            2,
            (unsigned int)decomposition->channel_count *
            (unsigned int)decomposition->lambda_count
        ))
    {
        return 0;
    }

    for (q_index = 0;
         q_index < decomposition->channel_count;
         ++q_index)
    {
        unsigned int q;

        q =
            decomposition->channel_ids[q_index];

        for (lambda_index = 0;
             lambda_index < decomposition->lambda_count;
             ++lambda_index)
        {
            N5V3 lambda;
            int leg_value;
            uint8_t assignment[2];

            lambda =
                decomposition->lambda_values[lambda_index];

            leg_value =
                n5_contraction_k2_walsh_leg_value(
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
                    n5_small_rational_make(
                        leg_value,
                        channels[q].norm
                    )
                ))
            {
                n5_sparse_factor_free(factor);
                return 0;
            }
        }
    }

    return 1;
}

/*
 * 对一条 partition factor 局部应用所有 incident k=2 Walsh legs。
 */
static int n5_contraction_transform_line_second_walsh(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    unsigned int line_id,
    const N5FixedFlowDomains* domains,
    N5SparseFactor* output
)
{
    unsigned int cells[4];
    unsigned int side;

    N5SparseFactor current;

    unsigned int index;

    if (!n5_contraction_line_cells(
            line_id,
            cells
        ))
    {
        return 0;
    }

    side =
        n5_contraction_line_side(
            line_id
        );

    if (!n5_contraction_build_line_sparse_factor(
            size_matrix,
            line_id,
            domains,
            &current
        ))
    {
        return 0;
    }

    for (index = 0; index < 4; ++index)
    {
        unsigned int cell;

        cell = cells[index];

        if (n5_contraction_cell_size(
                size_matrix,
                cell
            ) != 2)
        {
            continue;
        }

        {
            N5K2WalshDecomposition decomposition;
            N5SparseFactor leg;
            N5SparseFactor multiplied;
            N5SparseFactor summed;

            unsigned int q_var;

            if (!n5_contraction_build_k2_walsh_decomposition(
                    flow[cell],
                    &decomposition
                ))
            {
                n5_sparse_factor_free(&current);
                return 0;
            }

            if (!n5_contraction_build_k2_walsh_incidence_factor(
                    cell,
                    side,
                    &decomposition,
                    &leg
                ))
            {
                n5_sparse_factor_free(&current);
                return 0;
            }

            if (!n5_sparse_factor_multiply(
                    &current,
                    &leg,
                    &multiplied
                ))
            {
                n5_sparse_factor_free(&leg);
                n5_sparse_factor_free(&current);
                return 0;
            }

            q_var =
                n5_contraction_incidence_var(
                    cell,
                    side
                );

            if (!n5_sparse_factor_sum_out(
                    &multiplied,
                    q_var,
                    &summed
                ))
            {
                n5_sparse_factor_free(&multiplied);
                n5_sparse_factor_free(&leg);
                n5_sparse_factor_free(&current);
                return 0;
            }

            n5_sparse_factor_free(&multiplied);
            n5_sparse_factor_free(&leg);
            n5_sparse_factor_free(&current);

            current = summed;
        }
    }

    *output = current;

    return 1;
}

/*
 * 构造 exact second-Walsh rewritten fixed-flow sparse network。
 */
int n5_contraction_build_second_walsh_sparse_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5SparseNetwork* network
)
{
    N5FixedFlowDomains domains;

    int64_t total_scale;

    unsigned int cell;
    unsigned int line_id;

    if (network == NULL ||
        !n5_contraction_build_fixed_flow_domains(
            size_matrix,
            flow,
            &domains
        ))
    {
        return 0;
    }

    n5_sparse_network_init(network);

    total_scale = 1;

    /*
     * k=2 cell factors 被 decomposition 替换；
     * 其他 cell factors 原样保留。
     */
    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int k;

        k =
            n5_contraction_cell_size(
                size_matrix,
                cell
            );

        if (k == 2)
        {
            N5K2WalshDecomposition decomposition;

            if (!n5_contraction_build_k2_walsh_decomposition(
                    flow[cell],
                    &decomposition
                ))
            {
                n5_sparse_network_free(network);
                return 0;
            }

            total_scale *=
                decomposition.scale;

            continue;
        }

        {
            N5SparseFactor factor;

            if (!n5_contraction_build_cell_sparse_factor(
                    size_matrix,
                    cell,
                    &domains,
                    &factor
                ) ||
                !n5_sparse_network_add_move(
                    network,
                    &factor
                ))
            {
                n5_sparse_factor_free(&factor);
                n5_sparse_network_free(network);
                return 0;
            }
        }
    }

    /*
     * 12 条 partition lines：
     * 所有 incident k=2 q-index 被局部替换成 shared lambda_cell。
     */
    for (line_id = 0;
         line_id < N5_CONTRACTION_LINE_COUNT;
         ++line_id)
    {
        N5SparseFactor factor;

        if (!n5_contraction_transform_line_second_walsh(
                size_matrix,
                flow,
                line_id,
                &domains,
                &factor
            ) ||
            !n5_sparse_network_add_move(
                network,
                &factor
            ))
        {
            n5_sparse_factor_free(&factor);
            n5_sparse_network_free(network);
            return 0;
        }
    }

    /*
     * k=2 inverse-Gram 已吸收到 Walsh incidence legs。
     * 只保留 non-k2 metrics。
     */
    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        unsigned int k;
        unsigned int side;

        k =
            n5_contraction_cell_size(
                size_matrix,
                cell
            );

        if (k == 2)
        {
            continue;
        }

        for (side = 0;
             side <= N5_CONTRACTION_SIDE_C;
             ++side)
        {
            N5SparseFactor factor;

            if (!n5_contraction_build_metric_sparse_factor(
                    size_matrix,
                    cell,
                    side,
                    &domains,
                    &factor
                ) ||
                !n5_sparse_network_add_move(
                    network,
                    &factor
                ))
            {
                n5_sparse_factor_free(&factor);
                n5_sparse_network_free(network);
                return 0;
            }
        }
    }

    /*
     * 与 reference network 一致，global scale 作为 scalar factor。
     */
    if (total_scale != 1)
    {
        N5SparseFactor scale_factor;

        if (!n5_sparse_factor_make_scalar(
                n5_small_rational_make(
                    total_scale,
                    1
                ),
                &scale_factor
            ) ||
            !n5_sparse_network_add_move(
                network,
                &scale_factor
            ))
        {
            n5_sparse_factor_free(&scale_factor);
            n5_sparse_network_free(network);
            return 0;
        }
    }

    return 1;
}

/*
 * scalar-subspace backend 内部使用的 canonical subspace record。
 *
 * rows 保存相对于 scalar nullspace basis 的 RREF rows；
 * key 用于 deterministic sort / binary search。
 */
typedef struct
{
    uint8_t dimension;
    uint16_t rows[3];
    uint64_t key;
} N5FlowSubspaceRecord;

/*
 * 动态 subspace list。
 */
typedef struct
{
    N5FlowSubspaceRecord* data;
    unsigned int count;
    unsigned int capacity;
} N5FlowSubspaceList;

/*
 * 对 unsigned int 做 qsort 升序比较。
 */
static int n5_contraction_compare_unsigned(
    const void* left,
    const void* right
)
{
    unsigned int a;
    unsigned int b;

    a = *(const unsigned int*)left;
    b = *(const unsigned int*)right;

    if (a < b)
    {
        return -1;
    }

    if (a > b)
    {
        return 1;
    }

    return 0;
}

/*
 * 对 N5FlowSubspaceRecord 按 canonical key 排序。
 */
static int n5_contraction_compare_subspace_record(
    const void* left,
    const void* right
)
{
    const N5FlowSubspaceRecord* a;
    const N5FlowSubspaceRecord* b;

    a = (const N5FlowSubspaceRecord*)left;
    b = (const N5FlowSubspaceRecord*)right;

    if (a->key < b->key)
    {
        return -1;
    }

    if (a->key > b->key)
    {
        return 1;
    }

    return 0;
}

/*
 * 将 dim<=3 的 canonical RREF basis 编码成唯一 64-bit key。
 *
 * scalar nullity 在 N5 all-even 中最多为 8，
 * 因而每个 coordinate row 都安全落在 uint16_t 中。
 */
static uint64_t n5_contraction_subspace_key(
    unsigned int dimension,
    const uint16_t rows[3]
)
{
    return
        ((uint64_t)dimension << 48) |
        ((uint64_t)rows[0]) |
        ((uint64_t)rows[1] << 16) |
        ((uint64_t)rows[2] << 32);
}

/*
 * 对最多三个 GF(2)^d coordinate rows 做 reduced row echelon form。
 *
 * pivot columns 按 0..d-1 扫描；
 * 输出 rows 具有唯一 canonical 顺序。
 *
 * 返回：
 *   rank。
 */
static unsigned int n5_contraction_rref_coordinates(
    const uint16_t* input_rows,
    unsigned int input_count,
    unsigned int ambient_dimension,
    uint16_t output_rows[3]
)
{
    uint16_t rows[3];
    unsigned int row_count;
    unsigned int index;
    unsigned int rank;
    unsigned int column;

    row_count = 0;

    for (index = 0;
         index < input_count && index < 3;
         ++index)
    {
        if (input_rows[index] != 0)
        {
            rows[row_count] = input_rows[index];
            ++row_count;
        }
    }

    rank = 0;

    for (column = 0;
         column < ambient_dimension && rank < row_count;
         ++column)
    {
        unsigned int pivot_row;
        unsigned int search_row;

        pivot_row = row_count;

        for (search_row = rank;
             search_row < row_count;
             ++search_row)
        {
            if ((rows[search_row] &
                 (uint16_t)(1u << column)) != 0)
            {
                pivot_row = search_row;
                break;
            }
        }

        if (pivot_row == row_count)
        {
            continue;
        }

        {
            uint16_t temporary;

            temporary = rows[rank];
            rows[rank] = rows[pivot_row];
            rows[pivot_row] = temporary;
        }

        for (search_row = 0;
             search_row < row_count;
             ++search_row)
        {
            if (search_row != rank &&
                (rows[search_row] &
                 (uint16_t)(1u << column)) != 0)
            {
                rows[search_row] ^= rows[rank];
            }
        }

        ++rank;
    }

    for (index = 0; index < 3; ++index)
    {
        output_rows[index] =
            index < rank
                ? rows[index]
                : 0;
    }

    return rank;
}

/*
 * 向动态 subspace list 追加一个 canonical record。
 */
static int n5_contraction_subspace_list_append(
    N5FlowSubspaceList* list,
    unsigned int dimension,
    const uint16_t rows[3]
)
{
    N5FlowSubspaceRecord* enlarged;
    unsigned int new_capacity;
    N5FlowSubspaceRecord* record;

    if (list == NULL || rows == NULL)
    {
        return 0;
    }

    if (list->count == list->capacity)
    {
        new_capacity =
            list->capacity == 0
                ? 128u
                : 2u * list->capacity;

        enlarged =
            (N5FlowSubspaceRecord*)realloc(
                list->data,
                (size_t)new_capacity *
                    sizeof(N5FlowSubspaceRecord)
            );

        if (enlarged == NULL)
        {
            return 0;
        }

        list->data = enlarged;
        list->capacity = new_capacity;
    }

    record = &list->data[list->count];

    record->dimension = (uint8_t)dimension;
    record->rows[0] = rows[0];
    record->rows[1] = rows[1];
    record->rows[2] = rows[2];
    record->key =
        n5_contraction_subspace_key(
            dimension,
            rows
        );

    ++list->count;

    return 1;
}

/*
 * 对一个固定 pivot set 枚举全部 canonical RREF subspaces。
 *
 * RREF 中只允许在 non-pivot column 且 pivot<column 的位置
 * 放置自由 bit，这与 reference enumerator 的 convention 一致。
 */
static int n5_contraction_append_rref_subspaces_for_pivots(
    unsigned int ambient_dimension,
    unsigned int dimension,
    const unsigned int pivots[3],
    N5FlowSubspaceList* list
)
{
    unsigned int free_rows[24];
    unsigned int free_columns[24];
    unsigned int free_count;

    unsigned int column;
    unsigned int row;
    unsigned int bits;
    unsigned int bit_limit;

    free_count = 0;

    for (column = 0;
         column < ambient_dimension;
         ++column)
    {
        int is_pivot;

        is_pivot = 0;

        for (row = 0; row < dimension; ++row)
        {
            if (pivots[row] == column)
            {
                is_pivot = 1;
                break;
            }
        }

        if (is_pivot)
        {
            continue;
        }

        for (row = 0; row < dimension; ++row)
        {
            if (pivots[row] < column)
            {
                free_rows[free_count] = row;
                free_columns[free_count] = column;
                ++free_count;
            }
        }
    }

    if (free_count >= 31)
    {
        return 0;
    }

    bit_limit = 1u << free_count;

    for (bits = 0; bits < bit_limit; ++bits)
    {
        uint16_t rows[3];
        unsigned int free_index;

        rows[0] = 0;
        rows[1] = 0;
        rows[2] = 0;

        for (row = 0; row < dimension; ++row)
        {
            rows[row] =
                (uint16_t)(1u << pivots[row]);
        }

        for (free_index = 0;
             free_index < free_count;
             ++free_index)
        {
            if ((bits &
                 (1u << free_index)) != 0)
            {
                rows[free_rows[free_index]] |=
                    (uint16_t)(
                        1u << free_columns[free_index]
                    );
            }
        }

        if (!n5_contraction_subspace_list_append(
                list,
                dimension,
                rows
            ))
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 枚举 GF(2)^d 中全部 dimension<=3 subspaces，
 * 每个 subspace 恰好以唯一 RREF basis 出现一次。
 */
static int n5_contraction_enumerate_flow_subspaces(
    unsigned int ambient_dimension,
    N5FlowSubspaceList* list
)
{
    uint16_t zero_rows[3];
    unsigned int p0;
    unsigned int p1;
    unsigned int p2;

    if (list == NULL ||
        ambient_dimension > 8)
    {
        return 0;
    }

    list->data = NULL;
    list->count = 0;
    list->capacity = 0;

    zero_rows[0] = 0;
    zero_rows[1] = 0;
    zero_rows[2] = 0;

    if (!n5_contraction_subspace_list_append(
            list,
            0,
            zero_rows
        ))
    {
        return 0;
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        unsigned int pivots[3];

        pivots[0] = p0;
        pivots[1] = 0;
        pivots[2] = 0;

        if (!n5_contraction_append_rref_subspaces_for_pivots(
                ambient_dimension,
                1,
                pivots,
                list
            ))
        {
            free(list->data);
            list->data = NULL;
            return 0;
        }
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        for (p1 = p0 + 1;
             p1 < ambient_dimension;
             ++p1)
        {
            unsigned int pivots[3];

            pivots[0] = p0;
            pivots[1] = p1;
            pivots[2] = 0;

            if (!n5_contraction_append_rref_subspaces_for_pivots(
                    ambient_dimension,
                    2,
                    pivots,
                    list
                ))
            {
                free(list->data);
                list->data = NULL;
                return 0;
            }
        }
    }

    for (p0 = 0; p0 < ambient_dimension; ++p0)
    {
        for (p1 = p0 + 1;
             p1 < ambient_dimension;
             ++p1)
        {
            for (p2 = p1 + 1;
                 p2 < ambient_dimension;
                 ++p2)
            {
                unsigned int pivots[3];

                pivots[0] = p0;
                pivots[1] = p1;
                pivots[2] = p2;

                if (!n5_contraction_append_rref_subspaces_for_pivots(
                        ambient_dimension,
                        3,
                        pivots,
                        list
                    ))
                {
                    free(list->data);
                    list->data = NULL;
                    return 0;
                }
            }
        }
    }

    qsort(
        list->data,
        list->count,
        sizeof(N5FlowSubspaceRecord),
        n5_contraction_compare_subspace_record
    );

    return 1;
}

/*
 * 在已按 key 排序的 subspace list 中做 binary search。
 *
 * 返回：
 *   找到时返回 index；
 *   未找到时返回 N5_CONTRACTION_INVALID_ID。
 */
static unsigned int n5_contraction_find_subspace_key(
    const N5FlowSubspaceList* list,
    uint64_t key
)
{
    unsigned int low;
    unsigned int high;

    low = 0;
    high = list->count;

    while (low < high)
    {
        unsigned int middle;
        uint64_t middle_key;

        middle = low + ((high - low) >> 1u);
        middle_key = list->data[middle].key;

        if (middle_key < key)
        {
            low = middle + 1u;
        }
        else
        {
            high = middle;
        }
    }

    if (low < list->count &&
        list->data[low].key == key)
    {
        return low;
    }

    return N5_CONTRACTION_INVALID_ID;
}

/*
 * 收集固定 K 的全部顶层 stabilizer permutations。
 */
static unsigned int n5_contraction_collect_top_stabilizer(
    const N5SizeMatrix* size_matrix,
    N5CellPermutation output[N5_TOP_GROUP_SIZE]
)
{
    N5CellPermutation group[N5_TOP_GROUP_SIZE];
    unsigned int group_count;
    unsigned int count;
    unsigned int index;

    group_count =
        n5_generate_top_group(
            group
        );

    count = 0;

    for (index = 0; index < group_count; ++index)
    {
        N5SizeMatrix transformed;

        n5_apply_cell_permutation_to_size_matrix(
            size_matrix,
            &group[index],
            &transformed
        );

        if (n5_size_matrix_equal(
                size_matrix,
                &transformed
            ))
        {
            output[count] = group[index];
            ++count;
        }
    }

    return count;
}

/*
 * 对一个 16-bit scalar cell mask 应用顶层 cell permutation。
 */
static uint16_t n5_contraction_permute_cell_mask(
    uint16_t mask,
    const N5CellPermutation* permutation
)
{
    uint16_t output;
    unsigned int old_cell;

    output = 0;

    for (old_cell = 0;
         old_cell < N5_CONTRACTION_CELL_COUNT;
         ++old_cell)
    {
        if ((mask &
             (uint16_t)(1u << old_cell)) != 0)
        {
            unsigned int new_cell;

            new_cell =
                permutation->destination[old_cell];

            output |=
                (uint16_t)(1u << new_cell);
        }
    }

    return output;
}

/*
 * 根据 scalar nullspace basis 建立：
 *
 *   coordinate -> 16-bit cell mask
 *
 * 与反向：
 *
 *   cell mask -> coordinate。
 *
 * reverse map 长度固定为 2^16；
 * 不属于 scalar flow space 的 mask 保存 -1。
 */
static int n5_contraction_build_scalar_coordinate_maps(
    const N5ScalarFlowBasis* basis,
    uint16_t coord_to_mask[256],
    int* mask_to_coord,
    unsigned int* coord_count_out
)
{
    unsigned int coord_count;
    unsigned int coord;
    unsigned int mask_value;

    if (basis == NULL ||
        mask_to_coord == NULL ||
        coord_count_out == NULL ||
        basis->nullity > 8)
    {
        return 0;
    }

    coord_count =
        1u << basis->nullity;

    for (mask_value = 0;
         mask_value < (1u << 16);
         ++mask_value)
    {
        mask_to_coord[mask_value] = -1;
    }

    for (coord = 0;
         coord < coord_count;
         ++coord)
    {
        uint16_t mask;
        unsigned int basis_index;

        mask = 0;

        for (basis_index = 0;
             basis_index < basis->nullity;
             ++basis_index)
        {
            if ((coord &
                 (1u << basis_index)) != 0)
            {
                mask ^=
                    basis->basis_masks[basis_index];
            }
        }

        coord_to_mask[coord] = mask;

        if (mask_to_coord[mask] >= 0)
        {
            return 0;
        }

        mask_to_coord[mask] =
            (int)coord;
    }

    *coord_count_out = coord_count;

    return 1;
}

/*
 * 构造 Stab(K) 在 scalar flow coordinates GF(2)^d 上的
 * 不同 induced linear actions。
 *
 * 每个 action 以完整 lookup table 保存：
 *
 *   action[coord] = transformed_coord。
 */
static int n5_contraction_build_scalar_stabilizer_actions(
    const N5SizeMatrix* size_matrix,
    const N5ScalarFlowBasis* basis,
    N5CellPermutation stabilizer[N5_TOP_GROUP_SIZE],
    unsigned int* stabilizer_count_out,
    uint16_t** actions_out,
    unsigned int* action_count_out,
    unsigned int* coord_count_out,
    uint16_t coord_to_mask[256]
)
{
    int* mask_to_coord;
    unsigned int coord_count;
    unsigned int stabilizer_count;

    uint16_t* actions;
    unsigned int action_count;
    unsigned int perm_index;

    mask_to_coord =
        (int*)malloc(
            (size_t)(1u << 16) * sizeof(int)
        );

    if (mask_to_coord == NULL)
    {
        return 0;
    }

    if (!n5_contraction_build_scalar_coordinate_maps(
            basis,
            coord_to_mask,
            mask_to_coord,
            &coord_count
        ))
    {
        free(mask_to_coord);
        return 0;
    }

    stabilizer_count =
        n5_contraction_collect_top_stabilizer(
            size_matrix,
            stabilizer
        );

    if (stabilizer_count == 0)
    {
        free(mask_to_coord);
        return 0;
    }

    actions =
        (uint16_t*)malloc(
            (size_t)stabilizer_count *
            (size_t)coord_count *
            sizeof(uint16_t)
        );

    if (actions == NULL)
    {
        free(mask_to_coord);
        return 0;
    }

    action_count = 0;

    for (perm_index = 0;
         perm_index < stabilizer_count;
         ++perm_index)
    {
        uint16_t temporary[256];
        unsigned int coord;
        int valid;
        int duplicate;
        unsigned int existing;

        valid = 1;

        for (coord = 0;
             coord < coord_count;
             ++coord)
        {
            uint16_t transformed_mask;
            int transformed_coord;

            transformed_mask =
                n5_contraction_permute_cell_mask(
                    coord_to_mask[coord],
                    &stabilizer[perm_index]
                );

            transformed_coord =
                mask_to_coord[transformed_mask];

            if (transformed_coord < 0)
            {
                valid = 0;
                break;
            }

            temporary[coord] =
                (uint16_t)transformed_coord;
        }

        if (!valid)
        {
            free(actions);
            free(mask_to_coord);
            return 0;
        }

        duplicate = 0;

        for (existing = 0;
             existing < action_count;
             ++existing)
        {
            if (memcmp(
                    &actions[existing * coord_count],
                    temporary,
                    (size_t)coord_count * sizeof(uint16_t)
                ) == 0)
            {
                duplicate = 1;
                break;
            }
        }

        if (!duplicate)
        {
            memcpy(
                &actions[action_count * coord_count],
                temporary,
                (size_t)coord_count * sizeof(uint16_t)
            );

            ++action_count;
        }
    }

    free(mask_to_coord);

    *stabilizer_count_out =
        stabilizer_count;

    *actions_out =
        actions;

    *action_count_out =
        action_count;

    *coord_count_out =
        coord_count;

    return 1;
}

/*
 * spanning-triple multiplicity：
 *
 *   dim 0 -> 1
 *   dim 1 -> 7
 *   dim 2 -> 42
 *   dim 3 -> 168。
 */
static uint64_t n5_contraction_spanning_triple_count(
    unsigned int dimension
)
{
    switch (dimension)
    {
    case 0:
        return UINT64_C(1);
    case 1:
        return UINT64_C(7);
    case 2:
        return UINT64_C(42);
    case 3:
        return UINT64_C(168);
    default:
        return UINT64_C(0);
    }
}

/*
 * 从一个 canonical scalar subspace basis
 * 构造 concrete V3-valued representative flow。
 */
static void n5_contraction_representative_flow_from_subspace(
    const N5FlowSubspaceRecord* subspace,
    const uint16_t coord_to_mask[256],
    N5V3 flow[N5_CONTRACTION_CELL_COUNT]
)
{
    unsigned int cell;
    unsigned int component;

    for (cell = 0;
         cell < N5_CONTRACTION_CELL_COUNT;
         ++cell)
    {
        flow[cell] = 0;
    }

    for (component = 0;
         component < subspace->dimension;
         ++component)
    {
        uint16_t cell_mask;
        N5V3 component_value;

        cell_mask =
            coord_to_mask[
                subspace->rows[component]
            ];

        component_value =
            (N5V3)(1u << component);

        for (cell = 0;
             cell < N5_CONTRACTION_CELL_COUNT;
             ++cell)
        {
            if ((cell_mask &
                 (uint16_t)(1u << cell)) != 0)
            {
                flow[cell] ^=
                    component_value;
            }
        }
    }
}

void n5_contraction_reduced_flow_catalog_free(N5ReducedFlowCatalog* catalog)
{
    if (catalog != NULL)
    {
        free(catalog->classes);
        memset(catalog, 0, sizeof(*catalog));
    }
}

/*
 * 使用 scalar-subspace quotient 构造 reduced V3-flow catalog。
 */
int n5_contraction_reduce_v3_flows_subspace(
    const N5SizeMatrix* size_matrix,
    N5ReducedFlowCatalog* catalog
)
{
    N5ScalarFlowBasis basis;
    N5CellPermutation stabilizer[N5_TOP_GROUP_SIZE];

    uint16_t* actions;
    unsigned int stabilizer_count;
    unsigned int action_count;
    unsigned int coord_count;
    uint16_t coord_to_mask[256];

    N5FlowSubspaceList subspaces;
    unsigned char* visited;
    unsigned int* image_indices;

    unsigned int seed_index;

    if (catalog == NULL)
    {
        return 0;
    }
    memset(catalog, 0, sizeof(*catalog));

    if (!n5_contraction_build_scalar_flow_basis(
            size_matrix,
            &basis
        ) ||
        basis.nullity > 8)
    {
        return 0;
    }

    actions = NULL;

    if (!n5_contraction_build_scalar_stabilizer_actions(
            size_matrix,
            &basis,
            stabilizer,
            &stabilizer_count,
            &actions,
            &action_count,
            &coord_count,
            coord_to_mask
        ))
    {
        return 0;
    }

    if (!n5_contraction_enumerate_flow_subspaces(
            basis.nullity,
            &subspaces
        ))
    {
        free(actions);
        return 0;
    }

    visited =
        (unsigned char*)calloc(
            subspaces.count,
            sizeof(unsigned char)
        );

    image_indices =
        (unsigned int*)malloc(
            (size_t)action_count *
            sizeof(unsigned int)
        );

    if (visited == NULL ||
        image_indices == NULL)
    {
        free(image_indices);
        free(visited);
        free(subspaces.data);
        free(actions);
        return 0;
    }

    catalog->classes = (N5ReducedFlowClass*)calloc(
        (size_t)subspaces.count, sizeof(*catalog->classes));
    if (catalog->classes == NULL)
    {
        free(image_indices);
        free(visited);
        free(subspaces.data);
        free(actions);
        return 0;
    }

    catalog->scalar_nullity =
        basis.nullity;

    catalog->stabilizer_size =
        stabilizer_count;

    catalog->induced_action_count =
        action_count;

    catalog->subspace_count =
        subspaces.count;

    catalog->raw_flow_count =
        n5_contraction_count_v3_flows(
            &basis
        );

    catalog->weighted_flow_count = 0;
    catalog->class_count = 0;

    for (seed_index = 0;
         seed_index < subspaces.count;
         ++seed_index)
    {
        const N5FlowSubspaceRecord* seed;

        unsigned int action_index;
        unsigned int image_count;
        unsigned int unique_count;
        unsigned int image_index;
        unsigned int canonical_index;

        N5ReducedFlowClass* flow_class;
        uint64_t spanning_count;

        if (visited[seed_index] != 0)
        {
            continue;
        }

        seed =
            &subspaces.data[seed_index];

        image_count = 0;

        for (action_index = 0;
             action_index < action_count;
             ++action_index)
        {
            uint16_t transformed_input[3];
            uint16_t transformed_rows[3];
            unsigned int row;
            unsigned int transformed_rank;
            uint64_t transformed_key;
            unsigned int found_index;

            transformed_input[0] = 0;
            transformed_input[1] = 0;
            transformed_input[2] = 0;

            for (row = 0;
                 row < seed->dimension;
                 ++row)
            {
                transformed_input[row] =
                    actions[
                        action_index * coord_count +
                        seed->rows[row]
                    ];
            }

            transformed_rank =
                n5_contraction_rref_coordinates(
                    transformed_input,
                    seed->dimension,
                    basis.nullity,
                    transformed_rows
                );

            if (transformed_rank !=
                seed->dimension)
            {
                free(image_indices);
                free(visited);
                free(subspaces.data);
                free(actions);
                free(catalog->classes);
                memset(catalog, 0, sizeof(*catalog));
                return 0;
            }

            transformed_key =
                n5_contraction_subspace_key(
                    transformed_rank,
                    transformed_rows
                );

            found_index =
                n5_contraction_find_subspace_key(
                    &subspaces,
                    transformed_key
                );

            if (found_index ==
                N5_CONTRACTION_INVALID_ID)
            {
                free(image_indices);
                free(visited);
                free(subspaces.data);
                free(actions);
                free(catalog->classes);
                memset(catalog, 0, sizeof(*catalog));
                return 0;
            }

            image_indices[image_count] =
                found_index;

            ++image_count;
        }

        qsort(
            image_indices,
            image_count,
            sizeof(unsigned int),
            n5_contraction_compare_unsigned
        );

        unique_count = 0;
        canonical_index =
            N5_CONTRACTION_INVALID_ID;

        for (image_index = 0;
             image_index < image_count;
             ++image_index)
        {
            unsigned int current;

            current =
                image_indices[image_index];

            if (image_index > 0 &&
                current ==
                    image_indices[image_index - 1])
            {
                continue;
            }

            if (canonical_index ==
                N5_CONTRACTION_INVALID_ID)
            {
                canonical_index = current;
            }

            visited[current] = 1;
            ++unique_count;
        }

        if (canonical_index == N5_CONTRACTION_INVALID_ID)
        {
            free(image_indices);
            free(visited);
            free(subspaces.data);
            free(actions);
            free(catalog->classes);
            memset(catalog, 0, sizeof(*catalog));
            return 0;
        }

        flow_class =
            &catalog->classes[
                catalog->class_count
            ];

        flow_class->dimension =
            subspaces.data[canonical_index].dimension;

        flow_class->basis_coordinates[0] =
            subspaces.data[canonical_index].rows[0];

        flow_class->basis_coordinates[1] =
            subspaces.data[canonical_index].rows[1];

        flow_class->basis_coordinates[2] =
            subspaces.data[canonical_index].rows[2];

        flow_class->subspace_orbit_size =
            unique_count;

        spanning_count =
            n5_contraction_spanning_triple_count(
                flow_class->dimension
            );

        flow_class->flow_orbit_size =
            (uint64_t)unique_count *
            spanning_count;

        n5_contraction_representative_flow_from_subspace(
            &subspaces.data[canonical_index],
            coord_to_mask,
            flow_class->representative_flow
        );

        catalog->weighted_flow_count +=
            flow_class->flow_orbit_size;

        ++catalog->class_count;
    }

    free(image_indices);
    free(visited);
    free(subspaces.data);
    free(actions);

    return 1;
}


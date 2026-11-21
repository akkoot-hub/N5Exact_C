#ifndef N5_ODD_ADAPTIVE_H
#define N5_ODD_ADAPTIVE_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_odd_kernel.h"
#include "n5_sparse.h"
#include "n5_slicing.h"
#include "n5_crt.h"


/*
 * fixed-charge adaptive contraction 采用的实际表示/路径。
 *
 * DIRECT:
 *   ordinary OddRepresentationNetwork 直接 contraction。
 *
 * DIRECT_REWRITTEN:
 *   k=2 second-Walsh rewrite 后直接 contraction。
 *
 * SLICED:
 *   ordinary network structural slicing 后求和。
 *
 * SLICED_REWRITTEN:
 *   second-Walsh rewrite 后再 structural slicing。
 */
typedef enum
{
    N5_ODD_ADAPTIVE_MODE_ZERO = 0,
    N5_ODD_ADAPTIVE_MODE_DIRECT = 1,
    N5_ODD_ADAPTIVE_MODE_DIRECT_REWRITTEN = 2,
    N5_ODD_ADAPTIVE_MODE_SLICED = 3,
    N5_ODD_ADAPTIVE_MODE_SLICED_REWRITTEN = 4
} N5OddAdaptiveMode;

/*
 * 一个 fixed reduced-charge class 的 six-prime adaptive contraction 结果。
 *
 * original_* 描述 Stage31 ordinary normalized network。
 * selected_* 描述实际用于 contraction 的 ordinary/rewrite network。
 * slice_* 描述 outer structural slicing；direct path 的 raw assignments=1。
 */
typedef struct
{
    int success;
    N5OddAdaptiveMode mode;
    int rewrite_applied;

    uint32_t residues[N5_CRT_PRIME_COUNT];

    unsigned int original_factor_count;
    unsigned int original_primitive_nnz;
    uint64_t original_peak;

    unsigned int selected_factor_count;
    unsigned int selected_primitive_nnz;
    uint64_t selected_peak;

    uint8_t slice_variable_count;
    uint8_t slice_variables[N5_SLICING_MAX_VARIABLES];
    uint64_t raw_slice_assignments;
    uint64_t zero_slice_assignments;
    uint64_t sliced_peak;

    uint64_t max_dense_table_entries;
} N5OddAdaptiveClassResult;

/*
 * 把 Stage31 normalized fixed-charge odd network 中全部 k=2 cells
 * 做 exact second-Walsh structural rewrite。
 *
 * 关键等式：
 *
 *   T2(qA,qB,qC)/(normA normB normC)
 *     = scale * sum_lambda
 *         product_side [gauge(q_side) chi(lambda,s_q)/norm(q_side)]。
 *
 * 因而：
 *   - 删除 k=2 normalized cell factor；
 *   - 将三个 incidence q variables 分别局部吸收到对应 line factor；
 *   - 三条 legs 共用 lambda_cell=48+cell；
 *   - 保留 exact definition-derived gauge 和 scale。
 *
 * output 拥有 factor memory，调用者必须 n5_sparse_network_free(output)。
 */
int n5_odd_adaptive_build_second_walsh_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[16],
    const N5OddKernelLibrary* library,
    const N5SparseNetwork* ordinary_network,
    N5SparseNetwork* output
);

/*
 * generic fixed-charge production path：
 *
 *   ordinary network + plan
 *      -> 若 peak>target 且存在 k=2，尝试 exact second-Walsh rewrite
 *      -> 仅当 rewrite peak 严格更小才采用
 *      -> peak<=target 时 six-prime shared dense direct contraction
 *      -> 否则 greedy structural slicing
 *      -> 每个 conditioned slice 一次 shared six-prime dense contraction
 *      -> 精确求和 residues。
 *
 * 不检查 family/global orbit ID。
 */
int n5_odd_adaptive_contract_fixed_charge(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[16],
    const N5OddKernelLibrary* library,
    uint64_t target_peak,
    unsigned int max_slices,
    N5OddAdaptiveClassResult* result
);

#endif

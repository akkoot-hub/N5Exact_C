#ifndef N5_CONTRACTION_H
#define N5_CONTRACTION_H

#include <stdint.h>

#include "n5_core.h"

/*
 * Forward declaration of the sparse network type.
 *
 * n5_contraction.h only needs pointers to this type.
 * The complete definition lives in n5_sparse.h.
 * Keeping the headers decoupled avoids circular/stale IntelliSense parsing.
 */
#ifndef N5_SPARSE_H
typedef struct N5SparseNetwork N5SparseNetwork;
#endif

/*
 * fixed-K transformed factor network 中的三种 cell side。
 *
 * 数学意义：
 *
 *   side A 对应固定 alpha 的 row partition；
 *   side B 对应固定 beta 的 column partition；
 *   side C 对应固定 gamma=alpha XOR beta 的 xor-diagonal partition。
 */
#define N5_CONTRACTION_SIDE_A 0
#define N5_CONTRACTION_SIDE_B 1
#define N5_CONTRACTION_SIDE_C 2

/*
 * 顶层固定 K 网络中的基础数量。
 */
#define N5_CONTRACTION_CELL_COUNT (N5_V2_SIZE * N5_V2_SIZE)
#define N5_CONTRACTION_LINE_COUNT (3u * N5_V2_SIZE)
#define N5_CONTRACTION_INCIDENCE_COUNT (N5_CONTRACTION_CELL_COUNT * 3u)

/*
 * 完整 transformed fixed-K network 的基础 factor 数量。
 *
 *   16 个 cell factors；
 *   12 个 partition-line factors；
 *   48 个 inverse-Gram metric factors。
 */
#define N5_CONTRACTION_TOTAL_FACTOR_COUNT (N5_CONTRACTION_CELL_COUNT + N5_CONTRACTION_LINE_COUNT + N5_CONTRACTION_INCIDENCE_COUNT)

/*
 * fixed-flow 后单个 incidence variable 的最大 channel domain。
 *
 * 在 even canonical basis 中最大的是：
 *
 *   k=4, r=0 -> 14 channels。
 */
#define N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN ((N5_V3_SIZE * (N5_V3_SIZE - 1u) * (N5_V3_SIZE - 2u) * (N5_V3_SIZE - 3u)) / (4u * 3u * 2u * 1u * 5u))

/*
 * 无效 ID。
 */
#define N5_CONTRACTION_INVALID_ID 0xFFFFFFFFu

/*
 * 一个固定 all-even K 的 scalar GF(2) r-flow nullspace 基。
 *
 * 数学来源：
 *
 *   每个能够携带非零 Fourier charge 的 cell
 *   对应一个 scalar GF(2) variable。
 *
 *   k=0 和 k=8 只有 r=0 channel，
 *   因此不属于 active cells。
 *
 *   每条 top-level line 都施加：
 *
 *     XOR_{c in line} x_c = 0。
 *
 * 字段：
 *
 *   active_count：
 *     active cells 数量。
 *
 *   active_cells：
 *     active cell 的 global cell IDs。
 *
 *   rank：
 *     12 条 line equations 在 active variables 上的 GF(2) rank。
 *
 *   nullity：
 *     scalar flow space 维数。
 *
 *   basis_masks：
 *     scalar nullspace basis。
 *
 *     每个 basis vector 使用 16-bit global cell mask 表示；
 *     第 c 位为 1 表示该 basis vector 在 cell c 上取 1。
 */
typedef struct
{
    uint8_t active_count;
    uint8_t active_cells[N5_CONTRACTION_CELL_COUNT];

    uint8_t rank;
    uint8_t nullity;

    uint16_t basis_masks[N5_CONTRACTION_CELL_COUNT];
} N5ScalarFlowBasis;

/*
 * 一个 GL(3,2) x Stab(K) reduced V3-valued r-flow class。
 *
 * dimension：
 *   对应 scalar flow space 中 subspace U 的维数，范围 0..3。
 *
 * basis_coordinates：
 *   U 的 canonical GF(2) RREF basis，坐标相对于
 *   N5ScalarFlowBasis::basis_masks。
 *
 * subspace_orbit_size：
 *   U 在 Stab(K) 下的 scalar-subspace orbit 大小。
 *
 * flow_orbit_size：
 *   完整 GL(3,2) x Stab(K) 下代表的 raw V3-flow 数量。
 *
 * representative_flow：
 *   用 basis vector 0/1/2 分别映射到 V3 bit 1/2/4
 *   得到的一个 concrete 16-cell representative flow。
 */
typedef struct
{
    uint8_t dimension;

    uint16_t basis_coordinates[3];

    unsigned int subspace_orbit_size;
    uint64_t flow_orbit_size;

    N5V3 representative_flow[N5_CONTRACTION_CELL_COUNT];
} N5ReducedFlowClass;

/*
 * subspace backend 生成的完整 flow quotient catalog。
 *
 * scalar_nullity：
 *   scalar GF(2) flow space dimension d。
 *
 * stabilizer_size：
 *   顶层 Stab(K) 大小。
 *
 * induced_action_count：
 *   Stab(K) 在 scalar flow coordinate space 上诱导出的
 *   不同线性 action 数量。
 *
 * subspace_count：
 *   GL(3,2) 先约化后，需要考虑的全部 dim<=3 scalar subspaces 数量。
 *
 * raw_flow_count：
 *   8^d。
 *
 * weighted_flow_count：
 *   全部 reduced class flow_orbit_size 的和；
 *   正确时必须等于 raw_flow_count。
 */
typedef struct
{
    uint8_t scalar_nullity;

    unsigned int stabilizer_size;
    unsigned int induced_action_count;
    unsigned int subspace_count;

    uint64_t raw_flow_count;
    uint64_t weighted_flow_count;

    unsigned int class_count;
    N5ReducedFlowClass* classes;
} N5ReducedFlowCatalog;

/*
 * 一个 fixed V3-valued r-flow 下的 48 个 incidence domains。
 *
 * 数学意义：
 *
 *   每个 cell 的三条 incidence 共享同一个 Fourier charge：
 *
 *     r_{c,A}=r_{c,B}=r_{c,C}=r_c。
 *
 *   对 incidence variable：
 *
 *     var = 3*cell + side
 *
 *   只保留 canonical channel library 中满足：
 *
 *     channel.r == r_c
 *
 *   的 channel IDs。
 *
 * 字段：
 *
 *   cell_r：
 *     16 个 cell 的固定 V3 charge。
 *
 *   domain_size：
 *     48 个 incidence variables 的 restricted domain size。
 *
 *   channel_ids：
 *     对每个 incidence 保存 surviving channel IDs。
 *
 * 注意：
 *   channel ID 保持原 cardinality-k channel library 的全局编号，
 *   不在 fixed-flow 后重新编号。
 */
typedef struct
{
    N5V3 cell_r[N5_CONTRACTION_CELL_COUNT];

    uint8_t domain_size[N5_CONTRACTION_INCIDENCE_COUNT];

    uint8_t channel_ids
        [N5_CONTRACTION_INCIDENCE_COUNT]
        [N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN];
} N5FixedFlowDomains;

/*
 * 将一个 cell-side incidence 编码成全局变量 ID。
 *
 * 固定约定：
 *
 *   incidence_var = 3 * cell + side。
 *
 * 因此 48 个 incidence variables 的 ID 为：
 *
 *   0..47。
 *
 * 参数：
 *   cell - global cell ID，范围 0..15；
 *   side - A/B/C，分别为 0/1/2。
 *
 * 返回：
 *   合法时返回 0..47；
 *   输入非法时返回 N5_CONTRACTION_INVALID_ID。
 */
unsigned int n5_contraction_incidence_var(
    unsigned int cell,
    unsigned int side
);

/*
 * 返回指定 top-level line 对应的 side。
 *
 * line IDs：
 *
 *   0..3   -> A rows；
 *   4..7   -> B columns；
 *   8..11  -> C xor-diagonals。
 *
 * 返回：
 *   0/1/2；
 *   输入非法时返回 N5_CONTRACTION_INVALID_ID。
 */
unsigned int n5_contraction_line_side(
    unsigned int line_id
);

/*
 * 生成指定 top-level line 上的四个 cell IDs。
 *
 * 参数：
 *   line_id - 0..11；
 *   cells   - 长度至少为 4 的输出数组。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法。
 */
int n5_contraction_line_cells(
    unsigned int line_id,
    unsigned int cells[4]
);

/*
 * 判断一个固定 K 是否属于当前 all-even contraction domain。
 *
 * 条件：
 *
 *   1. K 满足完整 size-matrix 合法性；
 *   2. 16 个 cell size 全部为偶数。
 *
 * 返回：
 *   1 表示可进入 all-even contraction engine；
 *   0 表示不可。
 */
int n5_contraction_k_is_all_even_legal(
    const N5SizeMatrix* size_matrix
);

/*
 * 构造固定 all-even K 的 scalar GF(2) r-flow nullspace basis。
 *
 * cell-level Fourier charge 满足：
 *
 *   r_{c,A}=r_{c,B}=r_{c,C}=r_c。
 *
 * 每条 line 满足：
 *
 *   XOR_{c in line} r_c = 0。
 *
 * 三个 V3 bit component 共享同一个 scalar GF(2) system。
 *
 * 参数：
 *   size_matrix - 固定 all-even K；
 *   basis       - 保存 scalar nullspace 信息。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法。
 */
int n5_contraction_build_scalar_flow_basis(
    const N5SizeMatrix* size_matrix,
    N5ScalarFlowBasis* basis
);

/*
 * 根据 scalar nullity 计算 V3-valued r-flow 总数量。
 *
 * 若 scalar nullity = d，
 * 三个 bit component 独立，因此：
 *
 *   number of V3 flows = 8^d。
 *
 * 返回：
 *   8^nullity。
 */
uint64_t n5_contraction_count_v3_flows(
    const N5ScalarFlowBasis* basis
);

/*
 * 检查一个完整 V3-valued cell flow 是否满足 fixed-K 约束。
 *
 * 条件：
 *
 *   1. 每个 r_c ∈ V3；
 *   2. k=0/8 cells 强制 r_c=0；
 *   3. 每条 top-level line 满足：
 *
 *        XOR_{c in line} r_c = 0。
 *
 * 返回：
 *   1 表示合法；
 *   0 表示非法。
 */
int n5_contraction_v3_flow_is_valid(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT]
);

/*
 * 对一个 fixed V3 r-flow 构造 48 个 restricted channel domains。
 *
 * 对每个 incidence：
 *
 *   保留 channel.r == r_cell。
 *
 * channel IDs 不重新编号。
 *
 * 参数：
 *   size_matrix - fixed all-even K；
 *   flow        - 合法 V3-valued r-flow；
 *   domains     - 保存 48 个 restricted domains。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示 flow 非法或某个 incidence 没有 surviving channel。
 */
int n5_contraction_build_fixed_flow_domains(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5FixedFlowDomains* domains
);

/*
 * 为 fixed all-even K 与合法 r-flow 构造完整 raw sparse network。
 *
 * 网络包含：
 *
 *   16 cell factors；
 *   12 ordered-partition factors；
 *   48 inverse-Gram metrics。
 *
 * 所有 factor 在构造时已经限制到 fixed-flow channel domains。
 *
 * exact values 使用 N5SmallRational，仅用于 local/network
 * correctness 与 rewrite diagnostics；它不是最终 production
 * 大整数 contraction backend。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示构造失败。
 */
int n5_contraction_build_fixed_flow_sparse_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5SparseNetwork* network
);

/*
 * 对 fixed-flow sparse network 应用 exact k=2 second-Walsh rewrite。
 *
 * 对每个 k=2 cell：
 *
 *   - 删除原三腿 cell factor；
 *   - 删除三条 k=2 inverse-Gram unary metric；
 *   - 引入一个 shared lambda_cell variable，ID=48+cell；
 *   - 将每条 incident partition factor 的 q-index
 *     通过 leg(q,lambda)/norm(q) 局部求和；
 *   - 保存 definition-derived local scale 与 gauge sign。
 *
 * 非 k=2 primitives 保持 canonical channel basis。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示构造失败。
 */
int n5_contraction_build_second_walsh_sparse_network(
    const N5SizeMatrix* size_matrix,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    N5SparseNetwork* network
);

/*
 * 使用 GL(3,2)-first scalar-subspace method
 * 构造 fixed all-even K 的 reduced r-flow catalog。
 *
 * 数学步骤：
 *
 *   1. 计算 scalar nullspace F2^d；
 *   2. 枚举其中全部 dim<=3 subspaces U；
 *   3. 用 Stab(K) quotient 这些 subspaces；
 *   4. 对 dim(U)=m，乘 spanning-triple multiplicity：
 *
 *        m=0 -> 1
 *        m=1 -> 7
 *        m=2 -> 42
 *        m=3 -> 168；
 *
 *   5. 生成一个 concrete V3-valued representative flow。
 *
 * 这与直接枚举全部 8^d flows 数学等价，
 * 但避免大 flow space 的 raw enumeration。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法、容量不足或内部一致性失败。
 */
int n5_contraction_reduce_v3_flows_subspace(
    const N5SizeMatrix* size_matrix,
    N5ReducedFlowCatalog* catalog
);

void n5_contraction_reduced_flow_catalog_free(
    N5ReducedFlowCatalog* catalog
);



#endif

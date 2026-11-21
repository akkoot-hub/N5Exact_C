#ifndef N5_SPARSE_H
#define N5_SPARSE_H

#include <stdint.h>

#define N5_SPARSE_MAX_ARITY 24
#define N5_SPARSE_MAX_VARIABLES 64
#define N5_SPARSE_MAX_VALUE 255
#define N5_SPARSE_NETWORK_MAX_FACTORS 96

/*
 * 小型 exact rational，仅用于 local primitive / sparse-factor
 * 结构变换阶段。
 *
 * 约定：
 *   denominator 始终为正；
 *   fraction 始终约分；
 *   0 统一表示为 0/1。
 *
 * 注意：
 *   这不是最终 production 大整数后端。
 *   它只承担本阶段 local tensor / second-Walsh rewrite
 *   中规模很小的 exact rational 运算。
 */
typedef struct
{
    int64_t numerator;
    int64_t denominator;
} N5SmallRational;

/*
 * SparseFactor 的一个非零 entry。
 *
 * assignment 的前 factor.arity 个位置有效。
 */
typedef struct
{
    uint8_t assignment[N5_SPARSE_MAX_ARITY];
    N5SmallRational value;
} N5SparseEntry;

/*
 * 通用小型 sparse factor。
 *
 * vars：
 *   factor scope，变量 ID 均小于 64。
 *
 * entries：
 *   只保存非零 assignment。
 *
 * entries 内部按 assignment 的 packed key 排序；
 * 对外函数会维持这一 canonical form。
 */
typedef struct
{
    uint8_t arity;
    uint8_t vars[N5_SPARSE_MAX_ARITY];

    unsigned int nnz;
    unsigned int capacity;

    N5SparseEntry* entries;
} N5SparseFactor;

/*
 * 一组拥有 factor ownership 的 sparse network。
 */
typedef struct N5SparseNetwork
{
    unsigned int factor_count;
    N5SparseFactor factors[N5_SPARSE_NETWORK_MAX_FACTORS];
} N5SparseNetwork;

/*
 * 从 concrete sparse factor data 推导出的 structural plan。
 */
typedef struct
{
    uint8_t variable_count;
    uint8_t elimination_order[N5_SPARSE_MAX_VARIABLES];
    uint64_t estimated_peak;
} N5SparseContractionPlan;

/* exact rational helpers */
N5SmallRational n5_small_rational_make(int64_t numerator, int64_t denominator);
N5SmallRational n5_small_rational_add(N5SmallRational a, N5SmallRational b);
N5SmallRational n5_small_rational_multiply(N5SmallRational a, N5SmallRational b);
int n5_small_rational_is_zero(N5SmallRational value);
int n5_small_rational_equal(N5SmallRational a, N5SmallRational b);

/* factor lifecycle */
int n5_sparse_factor_init(
    N5SparseFactor* factor,
    const uint8_t* vars,
    unsigned int arity,
    unsigned int initial_capacity
);
void n5_sparse_factor_free(N5SparseFactor* factor);

/*
 * 追加一个 entry；允许暂时重复 assignment。
 * 调用 n5_sparse_factor_compress() 后会合并同 key 并删除零项。
 */
int n5_sparse_factor_append(
    N5SparseFactor* factor,
    const uint8_t* assignment,
    N5SmallRational value
);
int n5_sparse_factor_compress(N5SparseFactor* factor);

/* exact sparse algebra */
int n5_sparse_factor_multiply(
    const N5SparseFactor* a,
    const N5SparseFactor* b,
    N5SparseFactor* output
);
int n5_sparse_factor_sum_out(
    const N5SparseFactor* input,
    unsigned int variable,
    N5SparseFactor* output
);

/* scalar / clone helpers */
int n5_sparse_factor_make_scalar(
    N5SmallRational value,
    N5SparseFactor* factor
);
int n5_sparse_factor_clone_with_vars(
    const N5SparseFactor* source,
    const uint8_t* vars,
    N5SparseFactor* output
);

/* network lifecycle */
void n5_sparse_network_init(N5SparseNetwork* network);
void n5_sparse_network_free(N5SparseNetwork* network);
int n5_sparse_network_add_move(
    N5SparseNetwork* network,
    N5SparseFactor* factor
);
unsigned int n5_sparse_network_primitive_nnz(
    const N5SparseNetwork* network
);

/* structural planner inferred from actual sparse data */
int n5_sparse_choose_plan(
    const N5SparseNetwork* network,
    N5SparseContractionPlan* plan
);
int n5_sparse_plan_is_valid(
    const N5SparseContractionPlan* plan
);

#endif

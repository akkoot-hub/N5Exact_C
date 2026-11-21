#ifndef N5_PRODUCTION_H
#define N5_PRODUCTION_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_sparse.h"
#include "n5_contraction.h"
#include "n5_dense.h"

#define N5_PRODUCTION_MAX_SLICE_VARIABLES 6

/*
 * 只允许 second-Walsh lambda variables (48..63) 的 structural slice plan。
 *
 * 这样 slice variable 的 concrete values 都属于 V3 / <r_cell>，
 * 可以安全利用 fixed-flow 的 lower-GL(3,2) stabilizer 做 assignment quotient。
 */
typedef struct
{
    uint8_t variable_count;
    uint8_t variables[N5_PRODUCTION_MAX_SLICE_VARIABLES];

    uint64_t unsliced_peak;
    uint64_t estimated_peak;
    uint64_t raw_assignment_count;
    uint64_t target_peak;

    int target_reached;
} N5LambdaSlicePlan;

/*
 * 一个 symmetry-reduced lambda-slice assignment orbit。
 *
 * values：
 *   与 plan.variables 一一对应的 canonical concrete lambda values。
 *
 * orbit_size：
 *   在 fixed-flow lower-GL(3,2) stabilizer 下，该 assignment orbit 的大小。
 */
typedef struct
{
    uint8_t values[N5_PRODUCTION_MAX_SLICE_VARIABLES];
    unsigned int orbit_size;
} N5LambdaSliceOrbit;

/*
 * 一个 fixed flow 的完整 lambda-slice quotient catalog。
 *
 * orbits 使用动态内存；调用者最终必须调用：
 *
 *   n5_production_free_lambda_slice_catalog()。
 */
typedef struct
{
    N5LambdaSlicePlan plan;

    unsigned int symmetry_count;
    unsigned int orbit_count;

    uint64_t weighted_assignment_count;

    N5LambdaSliceOrbit* orbits;
} N5LambdaSliceCatalog;

/*
 * 多模数 symmetry-sliced dense contraction 结果。
 *
 * values[lane] 对应调用时 moduli[lane]。
 */
typedef struct
{
    uint8_t modulus_count;
    uint32_t values[N5_DENSE_MAX_MODULI];

    unsigned int representative_contractions;
    uint64_t represented_raw_assignments;
    uint64_t zero_representatives;

    uint64_t max_table_entries;
    double condition_seconds;
    double plan_seconds;
    double dense_seconds;
} N5SymmetrySlicedDenseMultiResult;

/*
 * 对 second-Walsh network 只从 lambda variables 48..63 中
 * greedy 选择 slice variables。
 *
 * 每一步选择使 structural peak 最小的 lambda variable；
 * 平局时选择 variable ID 更小者。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法。
 */
int n5_production_choose_lambda_slices(
    const N5SparseNetwork* network,
    uint64_t target_peak,
    unsigned int max_slices,
    N5LambdaSlicePlan* plan
);

/*
 * 构造 lambda slice assignments 在 fixed-flow lower-GL(3,2)
 * stabilizer 下的 exact quotient。
 *
 * lower map 必须逐 cell 固定当前 r-flow：
 *
 *   g(r_c) = r_c。
 *
 * 对 r_c != 0 的 lambda，值代表 quotient V3/<r_c>，
 * 所以变换后重新 canonicalize：
 *
 *   lambda -> min(g(lambda), g(lambda) XOR r_c)。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示 plan/flow/network 不一致或内存失败。
 */
int n5_production_build_lambda_slice_catalog(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const N5LambdaSlicePlan* plan,
    N5LambdaSliceCatalog* catalog
);

/*
 * 释放 catalog 的动态 orbit storage。
 */
void n5_production_free_lambda_slice_catalog(
    N5LambdaSliceCatalog* catalog
);

/*
 * 对 pure-lambda second-Walsh network 构造更强的 exact affine quotient。
 *
 * 候选 action 为：
 *
 *   lambda_c -> canonical( L(lambda_c) XOR t mod <r_c> )
 *
 * 其中：
 *
 *   L in GL(3,2),
 *   t in V3,
 *   L(r_c)=r_c for every cell。
 *
 * 与普通 lambda catalog 不同，本函数不会仅凭群论假设候选 action
 * 一定保持当前 factor network；它会逐 factor、逐 nonzero entry
 * 做 exact invariance test，只保留真正的 network automorphisms。
 *
 * 该路径特别适合 uniform-k=2 网络：
 * 所有 dynamic variables 都是 lambda variables，
 * 因此 translation symmetry 能显著压缩 slice representatives。
 */
int n5_production_build_affine_lambda_slice_catalog(
    const N5SparseNetwork* network,
    const N5V3 flow[N5_CONTRACTION_CELL_COUNT],
    const N5LambdaSlicePlan* plan,
    N5LambdaSliceCatalog* catalog
);

/*
 * 对一个 symmetry-reduced slice catalog，一次完成多个 CRT primes。
 *
 * 每个 slice representative：
 *
 *   condition 一次；
 *   choose elimination plan 一次；
 *   dense structural traversal 一次；
 *   arithmetic lanes 独立更新多个 residues。
 */
int n5_production_contract_lambda_slice_catalog_dense_multi(
    const N5SparseNetwork* network,
    const N5LambdaSliceCatalog* catalog,
    const uint32_t* moduli,
    unsigned int modulus_count,
    N5SymmetrySlicedDenseMultiResult* result
);

#endif

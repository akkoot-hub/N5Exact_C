#ifndef N5_ODD_PRODUCTION_H
#define N5_ODD_PRODUCTION_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_crt.h"
#include "n5_odd_adaptive.h"
#include "n5_odd_kernel.h"

#define N5_ODD_PRODUCTION_DEFAULT_TARGET_PEAK UINT64_C(262144)
#define N5_ODD_PRODUCTION_DEFAULT_MAX_SLICES 6u

/*
 * 一个 fixed odd-sector K 的 adaptive production 汇总结果。
 *
 * class_count：ChargeReducer 输出的 reduced charge class 数量。
 * raw_charge_count：全部 class weights 之和，等于 8^d。
 * worker_threads：本次 class-level parallel solve 的总参与线程数（含主线程）。
 * residues / crt：fixed-K contribution C(K) 的 six-prime residues 与 CRT 值。
 *
 * mode_count：按 reduced class 统计 zero/direct/direct-rewritten/
 *             sliced/sliced-rewritten 五种实际 contraction path。
 *
 * max_* / total_*：用于最终 memory / slicing audit，不参与数学求值。
 */
typedef struct
{
    unsigned int class_count;
    uint64_t raw_charge_count;
    unsigned int worker_threads;

    uint32_t residues[N5_CRT_PRIME_COUNT];
    N5CRTResult crt;

    unsigned int mode_count[5];

    unsigned int max_original_primitive_nnz;
    unsigned int max_selected_primitive_nnz;
    uint64_t max_original_peak;
    uint64_t max_selected_peak;
    uint64_t max_sliced_peak;

    uint64_t total_raw_slice_assignments;
    uint64_t total_zero_slice_assignments;
    uint64_t max_dense_table_entries;

    double wall_seconds;
} N5OddProductionResult;

/*
 * 对一个 fixed odd-sector K 做完整 adaptive production solve。
 *
 * pipeline：
 *   ChargeReducer
 *     -> reduced charge classes
 *     -> 全局 CPU worker pool parallel-for
 *     -> 每个 class 调用 Stage33 adaptive contractor
 *          ordinary
 *          -> optional k=2 second-Walsh
 *          -> optional structural slicing
 *          -> shared six-prime dense contraction
 *     -> 按 class index deterministic weighted residue aggregation
 *     -> six-prime CRT reconstruction of C(K)。
 *
 * 该函数不读取 global orbit ID、family ID 或历史 C(K)。
 * 唯一影响 fast path 的输入是当前 K / charge 所产生的 factor graph structure。
 *
 * requested_threads：
 *   0 表示自动使用硬件 logical processors；
 *   >0 表示请求总参与线程数（含调用线程）。
 */
int n5_odd_production_solve_adaptive_modular(
    const N5SizeMatrix* size_matrix,
    const N5OddKernelLibrary* library,
    unsigned int requested_threads,
    uint64_t target_peak,
    unsigned int max_slices,
    N5OddProductionResult* result
);

#endif

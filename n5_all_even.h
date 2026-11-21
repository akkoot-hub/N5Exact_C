#ifndef N5_ALL_EVEN_H
#define N5_ALL_EVEN_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_crt.h"

/*
 * generic all-even canonical-K exact solver 的 production result。
 *
 * residues：
 *   C(K) 在六个 certified odd primes 下的 exact residues。
 *
 * crt：
 *   六模 CRT reconstruction。
 *
 * decimal_value：
 *   reconstructed C(K) 的十进制字符串。
 *
 * diagnostics：
 *   flow quotient、lambda slicing、dense contraction 与多线程 wall time。
 */
typedef struct
{
    uint32_t residues[N5_CRT_PRIME_COUNT];
    N5CRTResult crt;
    char decimal_value[N5_CRT_DECIMAL_BUFFER];

    unsigned int reduced_flow_count;
    uint64_t raw_flow_count;
    uint64_t weighted_flow_count;

    uint64_t total_raw_lambda_slices;
    uint64_t total_slice_representatives;
    uint64_t total_dense_contractions;
    uint64_t total_shared_dense_eliminations;
    uint64_t total_zero_representatives;
    uint64_t max_dense_table_entries;

    unsigned int direct_flow_classes;
    unsigned int rewrite_selected_classes;
    unsigned int lambda_symmetry_sliced_classes;
    unsigned int generic_sliced_classes;
    unsigned int affine_lambda_sliced_classes;
    unsigned int max_slice_symmetry_count;

    unsigned int worker_threads; /* total participating threads, including the caller */
    double wall_seconds;
    double summed_task_seconds;
} N5AllEvenSolveResult;

/*
 * 从一个 canonical all-even size matrix K 出发，执行完整 exact solve：
 *
 *   scalar-flow nullspace
 *     -> GL(3,2) x Stab(K) reduced flows
 *     -> second-Walsh representation
 *     -> adaptive lambda slicing
 *     -> fixed-flow GL(3,2) slice quotient
 *     -> six-prime fixed-memory dense modular contraction
 *     -> deterministic flow-orbit weighting
 *     -> CRT reconstruction。
 *
 * requested_threads：
 *   0 表示自动检测当前可用逻辑处理器；
 *   >0 表示请求总参与线程数（含调用线程）。
 *
 * 不读取 orbit ID 或任何 frozen C(K) result。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示 K 非 all-even 或内部计算失败。
 */
int n5_all_even_solve(
    const N5SizeMatrix* size_matrix,
    unsigned int requested_threads,
    N5AllEvenSolveResult* result
);


/*
 * n5_all_even_solve() 的 verbose production variant。
 *
 * report_progress != 0 时，每个 reduced-flow worker 完成后打印一行
 * flow-level progress。完成顺序可能因线程调度不同而变化；
 * 最终 exact aggregation 仍严格按 class index 固定顺序进行。
 */
int n5_all_even_solve_ex(
    const N5SizeMatrix* size_matrix,
    unsigned int requested_threads,
    int report_progress,
    N5AllEvenSolveResult* result
);

#endif

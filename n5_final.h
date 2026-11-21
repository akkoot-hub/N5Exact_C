#ifndef N5_FINAL_H
#define N5_FINAL_H

#include <stdint.h>

#include "n5_crt.h"

/*
 * 最终 one-shot N5 production 的运行选项。
 *
 * requested_threads：
 *   0 表示自动使用当前可用逻辑处理器；
 *   >0 表示显式请求总参与线程数。
 *
 * print_catalog：
 *   非零时打印从输入定义重新生成的顶层轨道目录。
 *
 * print_orbit_results：
 *   非零时打印每个生成轨道的求值摘要。
 */
typedef struct
{
    unsigned int requested_threads;
    int print_catalog;
    int print_orbit_results;
} N5FinalRunOptions;

/*
 * 最终 one-shot N5 production 的汇总结果。
 *
 * 三组 residues 分别对应 all-even、odd 与全部生成轨道的 weighted sum。
 * 三个 CRT result 分别恢复两个 sector subtotal 与总计数。
 *
 * total_wall_seconds 从进入 n5_final_run() 后立即开始计时，直到最终聚合完成。
 */
typedef struct
{
    unsigned int processed_orbits;
    unsigned int processed_all_even_orbits;
    unsigned int processed_odd_orbits;

    unsigned int total_all_even_reduced_flows;
    unsigned int total_odd_reduced_classes;

    uint64_t total_odd_raw_charge_count;
    unsigned int odd_mode_count[5];
    uint64_t odd_total_raw_slice_assignments;
    uint64_t odd_total_zero_slice_assignments;

    uint64_t odd_max_original_peak;
    uint64_t odd_max_selected_peak;
    uint64_t odd_max_sliced_peak;
    uint64_t global_max_dense_table_entries;

    uint32_t all_even_residues[N5_CRT_PRIME_COUNT];
    uint32_t odd_residues[N5_CRT_PRIME_COUNT];
    uint32_t total_residues[N5_CRT_PRIME_COUNT];

    N5CRTResult all_even_crt;
    N5CRTResult odd_crt;
    N5CRTResult total_crt;

    char all_even_decimal[N5_CRT_DECIMAL_BUFFER];
    char odd_decimal[N5_CRT_DECIMAL_BUFFER];
    char total_decimal[N5_CRT_DECIMAL_BUFFER];
    char modulus_decimal[N5_CRT_DECIMAL_BUFFER];

    double catalog_seconds;
    double all_even_orbit_seconds;
    double odd_orbit_seconds;
    double total_wall_seconds;
} N5FinalRunResult;

/*
 * 从 n=5 顶层定义开始执行最终完整 production：
 *
 *   enumerate the six-parameter candidates
 *     -> retain legal size matrices
 *     -> generate the top symmetry action
 *     -> construct canonical top-level orbits in canonical order
 *     -> dispatch each generated orbit:
 *          all-even -> n5_all_even_solve()
 *          odd      -> n5_odd_production_solve_adaptive_modular()
 *     -> 乘 top orbit size
 *     -> deterministic six-prime residue aggregation
 *     -> CRT reconstruction of N5_all_even / N5_odd / N5。
 *
 * 该函数不读取 frozen orbit table、历史 C(K)、checkpoint 或 CSV。
 * 返回：
 *   1 表示完整计算成功；
 *   0 表示输入、内存分配、solver 或 CRT 计算无法完成。
 */
int n5_final_run(
    const N5FinalRunOptions* options,
    N5FinalRunResult* result
);

#endif

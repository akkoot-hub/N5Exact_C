#ifndef N5_ODD_NETWORK_H
#define N5_ODD_NETWORK_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_odd_kernel.h"
#include "n5_sparse.h"

#define N5_ODD_NETWORK_CELL_COUNT (N5_V2_SIZE * N5_V2_SIZE)
#define N5_ODD_NETWORK_INCIDENCE_COUNT (N5_ODD_NETWORK_CELL_COUNT * 3u)
#define N5_ODD_NETWORK_MAX_DOMAIN (N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY / 5u)

/*
 * 一个 fixed translation-charge assignment 下 48 个 subset-incidence
 * variables 的 unified channel domains。
 *
 * 对 cell c / side s：var = 3*c+s。
 * domain 只保留 KernelLibrary 中 channel.r == r_c 的 channel IDs。
 */
typedef struct
{
    uint8_t domain_size[N5_ODD_NETWORK_INCIDENCE_COUNT];
    uint8_t channel_ids[N5_ODD_NETWORK_INCIDENCE_COUNT]
                       [N5_ODD_NETWORK_MAX_DOMAIN];
} N5OddFixedChargeDomains;

/*
 * 一个 generic odd fixed-charge exact representation network。
 *
 * sparse network 包含：
 *   16 个 cell factors；
 *   12 个 ordered partition-line factors。
 *
 * 每个 subset-incidence 的 inverse-Gram 1/norm 已精确吸收到
 * 所属 cell factor 的对应 entry 中，因此无需额外保留 48 个 unary
 * metric factors。这样在 singleton-domain 网络中也不会先乘出巨大的
 * 未约分整数，再等待 metric factors 抵消。
 *
 * 数学上这与“16 cell + 12 line + 48 metric”完全等价；只是把每个
 * metric factor 与唯一 incident cell factor 预先相乘。
 *
 * 这是 definition-level normalized representation；不使用 orbit/family ID。
 */
typedef struct
{
    N5SparseNetwork sparse;
    N5OddFixedChargeDomains domains;

    unsigned int odd_cell_count;
    unsigned int even_cell_count;

    unsigned int line_odd_degree_0;
    unsigned int line_odd_degree_2;
    unsigned int line_odd_degree_4;

    unsigned int primitive_nnz;
} N5OddRepresentationNetwork;

/* 初始化/释放 representation network。 */
void n5_odd_network_init(N5OddRepresentationNetwork* network);
void n5_odd_network_free(N5OddRepresentationNetwork* network);

/*
 * 从 K + fixed charge + KernelLibrary 构造完整 normalized exact network。
 * 返回 1 成功，0 表示输入、charge、内存或 factor 构造失败。
 */
int n5_odd_network_build(
    const N5SizeMatrix* size_matrix,
    const N5V3 r_by_cell[N5_ODD_NETWORK_CELL_COUNT],
    const N5OddKernelLibrary* library,
    N5OddRepresentationNetwork* network
);

#endif

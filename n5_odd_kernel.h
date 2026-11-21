#ifndef N5_ODD_KERNEL_H
#define N5_ODD_KERNEL_H

#include <stdint.h>

#include "n5_core.h"
#include "n5_kernel.h"

#define N5_ODD_KERNEL_CARDINALITY_COUNT (N5_V3_SIZE + 1u)
#define N5_ODD_KERNEL_MAX_CHANNELS_PER_K N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
#define N5_ODD_KERNEL_FANO_COUNT (N5_V3_SIZE - 1u)

/*
 * odd/even 统一 subset-Fourier basis channel。
 *
 * 该结构把现有 all-even canonical channel 与 odd translation-Fourier
 * channel 投到同一个 definition-level 表示：一个 channel 只需要给出
 * 它在 concrete subset basis 上的非零 support 与整数系数。
 *
 * k：subset cardinality。
 * r：translation Fourier charge。
 * s：even channel 的 XOR charge；odd channel 固定为 0。
 * fano_label：k=3/5 odd channel 的 Fano normal n；其他情况为 0。
 * exceptional：k=4,s=0,r=0 的 exceptional plane-sum channel 标志。
 * norm：unnormalized orthogonal basis 的 norm-square。
 */
typedef struct
{
    uint8_t k;
    N5V3 r;
    N5V3 s;
    uint8_t fano_label;
    uint8_t exceptional;
    uint8_t support_count;
    uint8_t norm;
    N5SubsetMask masks[8];
    int8_t coefficients[8];
} N5OddBasisChannel;

/*
 * odd KernelLibrary 的 definition-level 核心数据。
 *
 * channels[k][q]：cardinality k 的完整 orthogonal basis。
 * channel_count[k]：对应 basis dimension，理论上应为 C(8,k)。
 * m3：flattened 7x7x7 Fano multiplicity tensor M3。
 * odd_scalar_k1 / odd_scalar_k7：k=1,7 zero-sum representative 的 T3 scalar。
 */
typedef struct
{
    unsigned int channel_count[N5_ODD_KERNEL_CARDINALITY_COUNT];
    N5OddBasisChannel channels[N5_ODD_KERNEL_CARDINALITY_COUNT]
                              [N5_ODD_KERNEL_MAX_CHANNELS_PER_K];
    uint8_t m3[N5_ODD_KERNEL_FANO_COUNT *
               N5_ODD_KERNEL_FANO_COUNT *
               N5_ODD_KERNEL_FANO_COUNT];
    unsigned int odd_scalar_k1;
    unsigned int odd_scalar_k7;
} N5OddKernelLibrary;

/*
 * 从当前 all-even basis provider 与 odd zero-sum representatives 出发，
 * 构造统一 KernelLibrary。
 *
 * 该函数不读取 frozen kernel table。
 */
int n5_odd_kernel_build_library(
    N5OddKernelLibrary* library
);

/*
 * 返回 M3(a,b,c)，其中 a,b,c in {1,...,7}。
 */
int n5_odd_kernel_m3(
    const N5OddKernelLibrary* library,
    unsigned int a,
    unsigned int b,
    unsigned int c
);

/*
 * 返回 M5(a,b,c) = 4 - 2 M3(a,b,c)。
 */
int n5_odd_kernel_m5(
    const N5OddKernelLibrary* library,
    unsigned int a,
    unsigned int b,
    unsigned int c
);

/*
 * 计算 ordered four-leg partition line 的 exact unnormalized transformed
 * kernel：
 *
 *   Lhat(q0,q1,q2,q3)
 *   = sum_{S0 sqcup S1 sqcup S2 sqcup S3 = V3}
 *       product_i b_qi(Si).
 *
 * 返回值是精确整数。
 */
int n5_odd_kernel_line4(
    const N5OddBasisChannel* q0,
    const N5OddBasisChannel* q1,
    const N5OddBasisChannel* q2,
    const N5OddBasisChannel* q3
);

#endif

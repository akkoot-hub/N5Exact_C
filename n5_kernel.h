#ifndef N5_KERNEL_H
#define N5_KERNEL_H

#include <stdint.h>

#include "n5_core.h"

#define N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY ((N5_V3_SIZE * (N5_V3_SIZE - 1u) * (N5_V3_SIZE - 2u) * (N5_V3_SIZE - 3u)) / (4u * 3u * 2u * 1u))
#define N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT ((N5_V3_SIZE * (N5_V3_SIZE - 1u)) / 2u)

/*
 * 一个 canonical unnormalized translation-Fourier channel。
 *
 * 数学对象：
 *
 *   对固定 cardinality k 的 subset translation orbit，
 *   选择最小 subset mask 作为 representative。
 *
 *   对满足 stabilizer 正交条件的 r ∈ V3，
 *   定义：
 *
 *     e_q =
 *       sum_{S = representative + t}
 *       (-1)^(r·t) |S>。
 *
 * 字段：
 *
 *   k：
 *     subset cardinality。
 *
 *   representative：
 *     translation orbit 中最小的 subset mask。
 *
 *   r：
 *     translation Fourier label。
 *
 *   s：
 *     representative 的 XOR charge sigma(S)。
 *
 *   orbit_size：
 *     translation orbit 大小。
 *
 *   stabilizer_size：
 *     translation stabilizer 大小。
 *
 *   orbit_masks：
 *     orbit 中实际出现的 subset masks。
 *
 *   coefficients：
 *     与 orbit_masks 一一对应的 Fourier coefficient，
 *     每个值均为 +1 或 -1。
 *
 *   norm：
 *     unnormalized basis vector 的 Gram norm。
 *     在该 canonical basis 中等于 orbit_size。
 */
typedef struct
{
    uint8_t k;

    N5SubsetMask representative;

    N5V3 r;
    N5V3 s;

    uint8_t orbit_size;
    uint8_t stabilizer_size;

    N5SubsetMask orbit_masks[8];
    int8_t coefficients[8];

    N5V3 stabilizer[8];

    uint8_t norm;
} N5EvenFourierChannel;

/*
 * 一个由 definition-level generator 构造出的
 * all-even ordered-partition tensor。
 *
 * sizes[leg]：
 *   四条腿的 subset cardinality。
 *
 * channel_counts[leg]：
 *   对应 cardinality 的 canonical Fourier channel 数量。
 *
 * ordered_partitions：
 *   满足指定四个 cardinalities 的 concrete ordered partitions 数量。
 *
 * entry_count：
 *   transformed tensor 的总 entry 数量，
 *
 *     product_leg channel_counts[leg]。
 *
 * values：
 *   row-major 扁平化保存的 exact integer tensor。
 *
 * 当前 V3/all-even 五种 motif 中，
 * raw partition 数最多为 2520，
 * 因此 int16_t 足够保存每个 transformed entry。
 */
typedef struct
{
    unsigned int sizes[4];
    unsigned int channel_counts[4];

    unsigned int ordered_partitions;
    unsigned int entry_count;

    int16_t* values;
} N5EvenPartitionTensor;

/*
 * 计算 V3 上的二进制点积。
 *
 * 数学定义：
 *
 *   r · x
 *   =
 *   parity(popcount(r & x))
 *
 * 结果属于 F2。
 *
 * 返回：
 *   0 或 1。
 */
unsigned int n5_v3_dot(
    N5V3 r,
    N5V3 x
);

/*
 * 计算 Walsh character：
 *
 *   chi_r(x) = (-1)^(r·x)。
 *
 * 返回：
 *   +1 或 -1。
 */
int n5_walsh_character(
    N5V3 r,
    N5V3 x
);

/*
 * 计算局部计数 T3(A,B,C)。
 *
 * 数学定义：
 *
 *   对 A,B,C ⊆ V3，
 *
 *   T3(A,B,C)
 *
 * 等于满足以下条件的双射 f : A -> B 的数量：
 *
 *   { a XOR f(a) : a ∈ A } = C。
 *
 * 返回：
 *   精确的 T3(A,B,C)。
 */
unsigned int n5_t3(
    N5SubsetMask a_mask,
    N5SubsetMask b_mask,
    N5SubsetMask c_mask
);

/*
 * 从定义生成固定 even cardinality k 的
 * canonical translation-Fourier channels。
 *
 * 支持：
 *
 *   k = 0,2,4,6,8。
 *
 * 生成过程：
 *
 *   1. 枚举 cardinality 为 k 的全部 subset；
 *   2. 按 V3 translation 分成 orbit；
 *   3. 每个 orbit 取最小 mask 为 representative；
 *   4. 计算 translation stabilizer H；
 *   5. 枚举满足
 *
 *        r·h = 0  对全部 h ∈ H
 *
 *      的 Fourier label r；
 *   6. 生成 unnormalized Fourier channel。
 *
 * 参数：
 *   k        - even cardinality；
 *   output   - 保存 channel 的数组。
 *              如果为 NULL，则只统计数量；
 *   capacity - output 的容量。
 *
 * 返回：
 *   生成的 channel 数量。
 *
 *   输入非法或容量不足时返回 0。
 */
unsigned int n5_generate_even_fourier_channels(
    unsigned int k,
    N5EvenFourierChannel* output,
    unsigned int capacity
);

/*
 * 从 raw T3 定义直接计算一个 transformed even cell
 * tensor coefficient。
 *
 * 数学定义：
 *
 *   T_hat(qA,qB,qC)
 *
 *   =
 *
 *   sum_{A,B,C}
 *     c_qA(A)
 *     c_qB(B)
 *     c_qC(C)
 *     T3(A,B,C)。
 *
 * 这里使用的是 canonical unnormalized Fourier basis，
 * 因此本函数不插入 inverse Gram factor。
 *
 * 参数：
 *   channel_a - A 腿 Fourier channel；
 *   channel_b - B 腿 Fourier channel；
 *   channel_c - C 腿 Fourier channel。
 *
 * 返回：
 *   精确整数 transformed coefficient。
 */
int n5_even_cell_tensor_value(
    const N5EvenFourierChannel* channel_a,
    const N5EvenFourierChannel* channel_b,
    const N5EvenFourierChannel* channel_c
);

/*
 * 对一个 k=4 plane Fourier channel 求对应的
 * nonzero Fano normal n ∈ {1,...,7}。
 *
 * 数学定义：
 *
 *   plane translation stabilizer 是一个 2D subspace H；
 *   normal n 是唯一满足
 *
 *     n · h = 0
 *
 *   对全部 h ∈ H 成立的非零向量。
 *
 * 参数：
 *   channel - k=4、s=0 的 plane channel。
 *
 * 返回：
 *   唯一 nonzero normal n；
 *   输入不是合法 plane channel 时返回 0。
 */
N5V3 n5_k4_plane_normal(
    const N5EvenFourierChannel* channel
);

/*
 * 从 ordered-partition 定义直接构造任意 all-even
 * 四腿 partition tensor。
 *
 * 支持的 cardinality：
 *
 *   0,2,4,6,8
 *
 * 且必须满足：
 *
 *   k0+k1+k2+k3 = 8。
 *
 * generator 的 source of truth：
 *
 *   逐个枚举 V3 的 4^8 个 point-colorings，
 *   保留 cardinality 正确的 ordered partitions，
 *   再使用 canonical Fourier coefficients
 *   累加 transformed tensor。
 *
 * 参数：
 *   sizes  - 四条腿的 ordered cardinalities；
 *   tensor - 保存构造结果。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法或内存分配失败。
 */
int n5_build_even_partition_tensor(
    const unsigned int sizes[4],
    N5EvenPartitionTensor* tensor
);

/*
 * 释放 n5_build_even_partition_tensor()
 * 为 tensor->values 分配的内存。
 *
 * 参数：
 *   tensor - 待释放的 tensor。
 *
 * 返回：
 *   无。
 */
void n5_free_even_partition_tensor(
    N5EvenPartitionTensor* tensor
);

/*
 * 读取 generic partition tensor 的一个 entry。
 *
 * 参数：
 *   tensor      - 已构造好的 tensor；
 *   q0,q1,q2,q3- 四条腿的 local channel index。
 *
 * 返回：
 *   对应 transformed tensor value。
 *
 * 若输入越界则返回 0。
 */
int n5_even_partition_tensor_value(
    const N5EvenPartitionTensor* tensor,
    unsigned int q0,
    unsigned int q1,
    unsigned int q2,
    unsigned int q3
);

#endif

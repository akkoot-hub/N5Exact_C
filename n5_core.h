#ifndef N5_CORE_H
#define N5_CORE_H

#include <stdint.h>

/*
 * V2 = F_2^2 中一个向量的整数编码。
 *
 * 合法取值：
 *   0, 1, 2, 3
 *
 * 二进制的两个 bit 直接对应 F_2^2 的两个坐标。
 */
typedef uint8_t N5V2;

/*
 * V3 = F_2^3 中一个向量的整数编码。
 *
 * 合法取值：
 *   0, 1, 2, 3, 4, 5, 6, 7
 *
 * 二进制的三个 bit 直接对应 F_2^3 的三个坐标。
 */
typedef uint8_t N5V3;

/*
 * V3 的一个子集。
 *
 * 使用一个 8-bit mask 表示。
 * 第 x 位为 1，当且仅当元素 x 属于该子集。
 */
typedef uint8_t N5SubsetMask;

/*
 * 4×4 顶层 size matrix K。
 *
 * 数学意义：
 *
 *   K = (k_{a,b}),
 *   a,b ∈ V2。
 *
 * matrix[a][b] 表示 high cell (a,b)
 * 中包含的低位 V3 元素数量。
 */
typedef struct
{
    uint8_t matrix[4][4];
} N5SizeMatrix;

/*
 * V2 中元素的总数。
 */
#define N5_V2_SIZE 4

 /*
  * V3 中元素的总数。
  */
#define N5_V3_SIZE 8

  /*
   * 一个 high cell 中最多可以包含的低位元素数量。
   */
#define N5_CELL_MAX_SIZE 8

   /*
    * V3 全集对应的 subset mask。
    */

/*
 * 判断一个整数编码是否属于 V3。
 *
 * 参数：
 *   value - 待检查的整数编码。
 *
 * 返回：
 *   1 表示 value 属于 V3；
 *   0 表示 value 不属于 V3。
 */
int n5_v3_is_valid(unsigned int value);

/*
 * 计算 V3 中两个向量的和。
 *
 * 数学意义：
 *   V3 = F_2^3，因此向量加法就是逐 bit XOR。
 *
 * 参数：
 *   a - 第一个 V3 向量；
 *   b - 第二个 V3 向量。
 *
 * 返回：
 *   a XOR b。
 */
N5V3 n5_v3_add(
    N5V3 a,
    N5V3 b
);

/*
 * 判断元素 x 是否属于给定的 V3 子集。
 *
 * 参数：
 *   subset - 8-bit subset mask；
 *   x      - V3 中的元素。
 *
 * 返回：
 *   1 表示 x 属于 subset；
 *   0 表示 x 不属于 subset。
 */
int n5_subset_contains(
    N5SubsetMask subset,
    N5V3 x
);

/*
 * 计算 V3 子集的元素个数。
 *
 * 数学意义：
 *   即 subset 的 cardinality |S|。
 *
 * 参数：
 *   subset - 8-bit subset mask。
 *
 * 返回：
 *   子集中的元素数量，范围为 0~8。
 */
unsigned int n5_subset_size(
    N5SubsetMask subset
);

/*
 * 计算 V3 子集的 XOR charge。
 *
 * 数学定义：
 *
 *   sigma(S) = XOR_{x in S} x
 *
 * 参数：
 *   subset - 8-bit subset mask。
 *
 * 返回：
 *   sigma(S)，属于 V3。
 */
N5V3 n5_subset_xor_charge(
    N5SubsetMask subset
);

/*
 * 对 V3 子集执行 translation。
 *
 * 数学定义：
 *
 *   S + t = {x XOR t : x in S}
 *
 * 参数：
 *   subset - 原始 V3 子集；
 *   t      - translation 向量。
 *
 * 返回：
 *   translation 后的 subset mask。
 */
N5SubsetMask n5_subset_translate(
    N5SubsetMask subset,
    N5V3 t
);

/*
 * 判断一个 4×4 size matrix 是否满足全部合法性条件。
 *
 * 数学条件：
 *
 *   1. 每个 cell size 位于 0~8；
 *   2. 每一行和为 8；
 *   3. 每一列和为 8；
 *   4. 每条 XOR diagonal 的和为 8。
 *
 * 参数：
 *   size_matrix - 待检查的 size matrix。
 *
 * 返回：
 *   1 表示 K 合法；
 *   0 表示 K 不合法。
 */
int n5_size_matrix_is_legal(
    const N5SizeMatrix* size_matrix
);

#endif

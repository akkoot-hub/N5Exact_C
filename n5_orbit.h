#ifndef N5_ORBIT_H
#define N5_ORBIT_H

#include <stdint.h>

#include "n5_core.h"

#define N5_TOP_CELL_COUNT (N5_V2_SIZE * N5_V2_SIZE)
#define N5_SIDE_PERMUTATION_COUNT (3u * 2u * 1u)
#define N5_TOP_GROUP_SIZE ((N5_V2_SIZE - 1u) * (N5_V2_SIZE - 2u) * \
                           N5_V2_SIZE * N5_V2_SIZE * N5_SIDE_PERMUTATION_COUNT)

/*
 * 顶层 size-matrix orbit 的分类。
 */
typedef enum
{
    N5_FAMILY_UNKNOWN = 0,
    N5_FAMILY_ALL_EVEN,
    N5_FAMILY_P6,
    N5_FAMILY_P8A,
    N5_FAMILY_P8B,
    N5_FAMILY_P10,
    N5_FAMILY_P16
} N5OrbitFamily;

/*
 * 一个合法 K 的 12 条顶层 line 的 parity profile。
 */
typedef struct
{
    unsigned int zero_odd_lines;
    unsigned int two_odd_lines;
    unsigned int four_odd_lines;
} N5ParityLineProfile;

/*
 * 一个作用在 16 个 high cells 上的置换。
 * destination[c] 表示旧 cell c 被送到的新 cell。
 */
typedef struct
{
    uint8_t destination[N5_TOP_CELL_COUNT];
} N5CellPermutation;

/*
 * 一个顶层 size-matrix orbit 的完整基础信息。
 */
typedef struct
{
    N5SizeMatrix canonical;
    unsigned int member_count;
    unsigned int stabilizer_size;
    N5OrbitFamily family;
} N5TopOrbit;

/*
 * 在 V2 上应用一个 GL(2,2) 线性变换。
 */
N5V2 n5_v2_apply_linear_map(
    N5V2 x,
    N5V2 image_e0,
    N5V2 image_e1
);

/*
 * 生成由 V2 线性变换、平移和三侧置换组成的顶层 cell permutations。
 */
unsigned int n5_generate_top_group(
    N5CellPermutation permutations[N5_TOP_GROUP_SIZE]
);

/*
 * 判断两个 16-cell permutation 是否完全相同。
 */
int n5_cell_permutation_equal(
    const N5CellPermutation* a,
    const N5CellPermutation* b
);

/*
 * 统计 permutation 数组中不同 permutation 的数量。
 */
unsigned int n5_count_unique_cell_permutations(
    const N5CellPermutation* permutations,
    unsigned int count
);

/*
 * 对一个 size matrix K 应用指定的 16-cell permutation。
 */
void n5_apply_cell_permutation_to_size_matrix(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutation,
    N5SizeMatrix* output
);

/*
 * 判断两个 size matrix 是否完全相同。
 */
int n5_size_matrix_equal(
    const N5SizeMatrix* a,
    const N5SizeMatrix* b
);

/*
 * 按 row-major 顺序比较两个 size matrix。
 */
int n5_size_matrix_lex_compare(
    const N5SizeMatrix* a,
    const N5SizeMatrix* b
);

/*
 * 计算一个 size matrix 在顶层群作用下的 canonical representative。
 */
void n5_size_matrix_canonicalize(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutations,
    unsigned int permutation_count,
    N5SizeMatrix* canonical
);

/*
 * 计算一个 size matrix 的 stabilizer 大小。
 */
unsigned int n5_size_matrix_stabilizer_size(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutations,
    unsigned int permutation_count
);

/*
 * 判断一个 size matrix 是否属于 all-even sector。
 */
int n5_size_matrix_is_all_even(
    const N5SizeMatrix* size_matrix
);

/*
 * 计算一个合法 K 的 parity-line profile。
 */
int n5_size_matrix_parity_line_profile(
    const N5SizeMatrix* size_matrix,
    N5ParityLineProfile* profile
);

/*
 * 根据合法 K 的 parity topology 自动分类。
 */
N5OrbitFamily n5_classify_size_matrix_family(
    const N5SizeMatrix* size_matrix
);

/*
 * 返回 orbit family 的可读名称。
 */
const char* n5_orbit_family_name(
    N5OrbitFamily family
);

/*
 * 从全部合法 K 构造顶层 orbit catalog。
 */
unsigned int n5_build_top_orbits(
    const N5SizeMatrix* matrices,
    unsigned int matrix_count,
    const N5CellPermutation* permutations,
    unsigned int permutation_count,
    N5TopOrbit* orbits,
    unsigned int orbit_capacity
);

#endif

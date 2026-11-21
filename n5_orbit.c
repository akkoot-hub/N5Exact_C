#include <stddef.h>
#include "n5_orbit.h"

/*
 * 把 high cell (alpha,beta) 编码成 0~15。
 */
static unsigned int n5_cell_encode(
    N5V2 alpha,
    N5V2 beta
)
{
    return (unsigned int)alpha * N5_V2_SIZE + (unsigned int)beta;
}

/*
 * 对三元组的三个坐标执行一个 S3 置换。
 */
static void n5_apply_side_permutation(
    unsigned int permutation_index,
    N5V2 x0,
    N5V2 x1,
    N5V2 x2,
    N5V2* y0,
    N5V2* y1,
    N5V2* y2
)
{
    switch (permutation_index)
    {
    case 0:
        *y0 = x0; *y1 = x1; *y2 = x2;
        break;
    case 1:
        *y0 = x0; *y1 = x2; *y2 = x1;
        break;
    case 2:
        *y0 = x1; *y1 = x0; *y2 = x2;
        break;
    case 3:
        *y0 = x1; *y1 = x2; *y2 = x0;
        break;
    case 4:
        *y0 = x2; *y1 = x0; *y2 = x1;
        break;
    default:
        *y0 = x2; *y1 = x1; *y2 = x0;
        break;
    }
}

/*
 * 根据一条 line 上的 odd cell 数量更新 parity-line profile。
 */
static int n5_parity_profile_add_line(
    N5ParityLineProfile* profile,
    unsigned int odd_count
)
{
    if (profile == NULL)
    {
        return 0;
    }

    switch (odd_count)
    {
    case 0:
        ++profile->zero_odd_lines;
        return 1;
    case 2:
        ++profile->two_odd_lines;
        return 1;
    case 4:
        ++profile->four_odd_lines;
        return 1;
    default:
        return 0;
    }
}

/*
 * 将 orbit catalog 按 canonical K 的 row-major 字典序排序。
 */
static void n5_sort_top_orbits(
    N5TopOrbit* orbits,
    unsigned int orbit_count
)
{
    unsigned int i;

    for (i = 1; i < orbit_count; ++i)
    {
        N5TopOrbit current;
        unsigned int j;

        current = orbits[i];
        j = i;

        while (j > 0 &&
               n5_size_matrix_lex_compare(
                   &current.canonical,
                   &orbits[j - 1].canonical
               ) < 0)
        {
            orbits[j] = orbits[j - 1];
            --j;
        }

        orbits[j] = current;
    }
}

/*
 * 在 V2 上应用一个 GL(2,2) 线性变换。
 */
N5V2 n5_v2_apply_linear_map(
    N5V2 x,
    N5V2 image_e0,
    N5V2 image_e1
)
{
    N5V2 result;

    result = 0;

    if ((x & 1u) != 0)
    {
        result ^= image_e0;
    }

    if ((x & 2u) != 0)
    {
        result ^= image_e1;
    }

    return result;
}

/*
 * 根据 V2 线性变换、平移和三侧置换生成顶层 cell permutations。
 */
unsigned int n5_generate_top_group(
    N5CellPermutation permutations[N5_TOP_GROUP_SIZE]
)
{
    unsigned int count;
    N5V2 image_e0;
    N5V2 image_e1;
    unsigned int u;
    unsigned int v;
    unsigned int side;
    unsigned int alpha;
    unsigned int beta;

    count = 0;

    for (image_e0 = 1; image_e0 < N5_V2_SIZE; ++image_e0)
    {
        for (image_e1 = 1; image_e1 < N5_V2_SIZE; ++image_e1)
        {
            if (image_e0 == image_e1)
            {
                continue;
            }

            for (u = 0; u < N5_V2_SIZE; ++u)
            {
                for (v = 0; v < N5_V2_SIZE; ++v)
                {
                    for (side = 0; side < N5_SIDE_PERMUTATION_COUNT; ++side)
                    {
                        N5CellPermutation* current;

                        current = &permutations[count];

                        for (alpha = 0; alpha < N5_V2_SIZE; ++alpha)
                        {
                            for (beta = 0; beta < N5_V2_SIZE; ++beta)
                            {
                                N5V2 gamma;
                                N5V2 transformed_alpha;
                                N5V2 transformed_beta;
                                N5V2 transformed_gamma;
                                N5V2 side_alpha;
                                N5V2 side_beta;
                                N5V2 side_gamma;
                                unsigned int old_cell;
                                unsigned int new_cell;

                                gamma = (N5V2)(alpha ^ beta);

                                transformed_alpha =
                                    (N5V2)(
                                        n5_v2_apply_linear_map(
                                            (N5V2)alpha,
                                            image_e0,
                                            image_e1
                                        )
                                        ^ (N5V2)u
                                        ^ (N5V2)v
                                    );

                                transformed_beta =
                                    (N5V2)(
                                        n5_v2_apply_linear_map(
                                            (N5V2)beta,
                                            image_e0,
                                            image_e1
                                        )
                                        ^ (N5V2)u
                                    );

                                transformed_gamma =
                                    (N5V2)(
                                        n5_v2_apply_linear_map(
                                            gamma,
                                            image_e0,
                                            image_e1
                                        )
                                        ^ (N5V2)v
                                    );

                                n5_apply_side_permutation(
                                    side,
                                    transformed_alpha,
                                    transformed_beta,
                                    transformed_gamma,
                                    &side_alpha,
                                    &side_beta,
                                    &side_gamma
                                );

                                (void)side_gamma;

                                old_cell =
                                    n5_cell_encode(
                                        (N5V2)alpha,
                                        (N5V2)beta
                                    );

                                new_cell =
                                    n5_cell_encode(
                                        side_alpha,
                                        side_beta
                                    );

                                current->destination[old_cell] =
                                    (uint8_t)new_cell;
                            }
                        }

                        ++count;
                    }
                }
            }
        }
    }

    return count;
}

/*
 * 判断两个 16-cell permutation 是否完全相同。
 */
int n5_cell_permutation_equal(
    const N5CellPermutation* a,
    const N5CellPermutation* b
)
{
    unsigned int cell;

    if (a == NULL || b == NULL)
    {
        return 0;
    }

    for (cell = 0; cell < N5_TOP_CELL_COUNT; ++cell)
    {
        if (a->destination[cell] != b->destination[cell])
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 统计 permutation 数组中不同 permutation 的数量。
 */
unsigned int n5_count_unique_cell_permutations(
    const N5CellPermutation* permutations,
    unsigned int count
)
{
    unsigned int unique_count;
    unsigned int i;
    unsigned int j;
    int already_seen;

    if (permutations == NULL)
    {
        return 0;
    }

    unique_count = 0;

    for (i = 0; i < count; ++i)
    {
        already_seen = 0;

        for (j = 0; j < i; ++j)
        {
            if (n5_cell_permutation_equal(
                    &permutations[i],
                    &permutations[j]
                ))
            {
                already_seen = 1;
                break;
            }
        }

        if (!already_seen)
        {
            ++unique_count;
        }
    }

    return unique_count;
}

/*
 * 对 size matrix K 应用指定的 cell permutation。
 */
void n5_apply_cell_permutation_to_size_matrix(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutation,
    N5SizeMatrix* output
)
{
    unsigned int old_cell;

    if (input == NULL || permutation == NULL || output == NULL)
    {
        return;
    }

    for (old_cell = 0; old_cell < N5_TOP_CELL_COUNT; ++old_cell)
    {
        unsigned int new_cell;
        unsigned int old_row;
        unsigned int old_column;
        unsigned int new_row;
        unsigned int new_column;

        new_cell = permutation->destination[old_cell];

        old_row = old_cell >> 2u;
        old_column = old_cell & 3u;
        new_row = new_cell >> 2u;
        new_column = new_cell & 3u;

        output->matrix[new_row][new_column] =
            input->matrix[old_row][old_column];
    }
}

/*
 * 判断两个 size matrix 是否完全相同。
 */
int n5_size_matrix_equal(
    const N5SizeMatrix* a,
    const N5SizeMatrix* b
)
{
    return n5_size_matrix_lex_compare(a, b) == 0;
}

/*
 * 按 row-major 顺序比较两个 size matrix。
 */
int n5_size_matrix_lex_compare(
    const N5SizeMatrix* a,
    const N5SizeMatrix* b
)
{
    unsigned int row;
    unsigned int column;

    if (a == NULL || b == NULL)
    {
        return 0;
    }

    for (row = 0; row < N5_V2_SIZE; ++row)
    {
        for (column = 0; column < N5_V2_SIZE; ++column)
        {
            if (a->matrix[row][column] < b->matrix[row][column])
            {
                return -1;
            }

            if (a->matrix[row][column] > b->matrix[row][column])
            {
                return 1;
            }
        }
    }

    return 0;
}

/*
 * 计算 K 在全部顶层群作用下的 canonical representative。
 */
void n5_size_matrix_canonicalize(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutations,
    unsigned int permutation_count,
    N5SizeMatrix* canonical
)
{
    unsigned int index;
    N5SizeMatrix transformed;

    if (input == NULL || permutations == NULL ||
        canonical == NULL || permutation_count == 0)
    {
        return;
    }

    *canonical = *input;

    for (index = 0; index < permutation_count; ++index)
    {
        n5_apply_cell_permutation_to_size_matrix(
            input,
            &permutations[index],
            &transformed
        );

        if (n5_size_matrix_lex_compare(
                &transformed,
                canonical
            ) < 0)
        {
            *canonical = transformed;
        }
    }
}

/*
 * 计算 K 的 stabilizer 大小。
 */
unsigned int n5_size_matrix_stabilizer_size(
    const N5SizeMatrix* input,
    const N5CellPermutation* permutations,
    unsigned int permutation_count
)
{
    unsigned int stabilizer_size;
    unsigned int index;
    N5SizeMatrix transformed;

    if (input == NULL || permutations == NULL)
    {
        return 0;
    }

    stabilizer_size = 0;

    for (index = 0; index < permutation_count; ++index)
    {
        n5_apply_cell_permutation_to_size_matrix(
            input,
            &permutations[index],
            &transformed
        );

        if (n5_size_matrix_equal(input, &transformed))
        {
            ++stabilizer_size;
        }
    }

    return stabilizer_size;
}

/*
 * 判断一个 size matrix 是否属于 all-even sector。
 */
int n5_size_matrix_is_all_even(
    const N5SizeMatrix* size_matrix
)
{
    unsigned int row;
    unsigned int column;

    if (size_matrix == NULL)
    {
        return 0;
    }

    for (row = 0; row < N5_V2_SIZE; ++row)
    {
        for (column = 0; column < N5_V2_SIZE; ++column)
        {
            if ((size_matrix->matrix[row][column] & 1u) != 0)
            {
                return 0;
            }
        }
    }

    return 1;
}

/*
 * 计算一个合法 K 的 12-line parity profile。
 */
int n5_size_matrix_parity_line_profile(
    const N5SizeMatrix* size_matrix,
    N5ParityLineProfile* profile
)
{
    unsigned int fixed;
    unsigned int variable;

    if (size_matrix == NULL || profile == NULL)
    {
        return 0;
    }

    if (!n5_size_matrix_is_legal(size_matrix))
    {
        return 0;
    }

    profile->zero_odd_lines = 0;
    profile->two_odd_lines = 0;
    profile->four_odd_lines = 0;

    for (fixed = 0; fixed < N5_V2_SIZE; ++fixed)
    {
        unsigned int odd_count;
        odd_count = 0;

        for (variable = 0; variable < N5_V2_SIZE; ++variable)
        {
            odd_count += size_matrix->matrix[fixed][variable] & 1u;
        }

        if (!n5_parity_profile_add_line(profile, odd_count))
        {
            return 0;
        }
    }

    for (fixed = 0; fixed < N5_V2_SIZE; ++fixed)
    {
        unsigned int odd_count;
        odd_count = 0;

        for (variable = 0; variable < N5_V2_SIZE; ++variable)
        {
            odd_count += size_matrix->matrix[variable][fixed] & 1u;
        }

        if (!n5_parity_profile_add_line(profile, odd_count))
        {
            return 0;
        }
    }

    for (fixed = 0; fixed < N5_V2_SIZE; ++fixed)
    {
        unsigned int odd_count;
        odd_count = 0;

        for (variable = 0; variable < N5_V2_SIZE; ++variable)
        {
            unsigned int beta;
            beta = variable ^ fixed;
            odd_count += size_matrix->matrix[variable][beta] & 1u;
        }

        if (!n5_parity_profile_add_line(profile, odd_count))
        {
            return 0;
        }
    }

    return profile->zero_odd_lines +
           profile->two_odd_lines +
           profile->four_odd_lines == 12;
}

/*
 * 根据合法 K 的 parity-line topology 自动判断其 family。
 */
N5OrbitFamily n5_classify_size_matrix_family(
    const N5SizeMatrix* size_matrix
)
{
    N5ParityLineProfile profile;

    if (size_matrix == NULL)
    {
        return N5_FAMILY_UNKNOWN;
    }

    if (!n5_size_matrix_is_legal(size_matrix))
    {
        return N5_FAMILY_UNKNOWN;
    }

    if (n5_size_matrix_is_all_even(size_matrix))
    {
        return N5_FAMILY_ALL_EVEN;
    }

    if (!n5_size_matrix_parity_line_profile(size_matrix, &profile))
    {
        return N5_FAMILY_UNKNOWN;
    }

    if (profile.zero_odd_lines == 3 &&
        profile.two_odd_lines == 9 &&
        profile.four_odd_lines == 0)
    {
        return N5_FAMILY_P6;
    }

    if (profile.zero_odd_lines == 0 &&
        profile.two_odd_lines == 12 &&
        profile.four_odd_lines == 0)
    {
        return N5_FAMILY_P8A;
    }

    if (profile.zero_odd_lines == 2 &&
        profile.two_odd_lines == 8 &&
        profile.four_odd_lines == 2)
    {
        return N5_FAMILY_P8B;
    }

    if (profile.zero_odd_lines == 0 &&
        profile.two_odd_lines == 9 &&
        profile.four_odd_lines == 3)
    {
        return N5_FAMILY_P10;
    }

    if (profile.zero_odd_lines == 0 &&
        profile.two_odd_lines == 0 &&
        profile.four_odd_lines == 12)
    {
        return N5_FAMILY_P16;
    }

    return N5_FAMILY_UNKNOWN;
}

/*
 * 返回 orbit family 的可读名称。
 */
const char* n5_orbit_family_name(
    N5OrbitFamily family
)
{
    switch (family)
    {
    case N5_FAMILY_ALL_EVEN:
        return "all-even";
    case N5_FAMILY_P6:
        return "P6";
    case N5_FAMILY_P8A:
        return "P8a";
    case N5_FAMILY_P8B:
        return "P8b";
    case N5_FAMILY_P10:
        return "P10";
    case N5_FAMILY_P16:
        return "P16";
    default:
        return "unknown";
    }
}

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
)
{
    unsigned int orbit_count;
    unsigned int matrix_index;

    if (matrices == NULL || permutations == NULL ||
        orbits == NULL || orbit_capacity == 0)
    {
        return 0;
    }

    orbit_count = 0;

    for (matrix_index = 0; matrix_index < matrix_count; ++matrix_index)
    {
        N5SizeMatrix canonical;
        unsigned int orbit_index;
        int found;

        n5_size_matrix_canonicalize(
            &matrices[matrix_index],
            permutations,
            permutation_count,
            &canonical
        );

        found = 0;

        for (orbit_index = 0; orbit_index < orbit_count; ++orbit_index)
        {
            if (n5_size_matrix_equal(
                    &canonical,
                    &orbits[orbit_index].canonical
                ))
            {
                ++orbits[orbit_index].member_count;
                found = 1;
                break;
            }
        }

        if (!found)
        {
            if (orbit_count >= orbit_capacity)
            {
                return 0;
            }

            orbits[orbit_count].canonical = canonical;
            orbits[orbit_count].member_count = 1;
            orbits[orbit_count].stabilizer_size = 0;
            orbits[orbit_count].family = N5_FAMILY_UNKNOWN;
            ++orbit_count;
        }
    }

    for (matrix_index = 0; matrix_index < orbit_count; ++matrix_index)
    {
        orbits[matrix_index].stabilizer_size =
            n5_size_matrix_stabilizer_size(
                &orbits[matrix_index].canonical,
                permutations,
                permutation_count
            );

        orbits[matrix_index].family =
            n5_classify_size_matrix_family(
                &orbits[matrix_index].canonical
            );
    }

    n5_sort_top_orbits(orbits, orbit_count);

    return orbit_count;
}

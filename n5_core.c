#include <stdio.h>

#include "n5_core.h"

/*
 * 计算 size matrix 指定一行的元素总和。
 */
static unsigned int n5_size_matrix_row_sum(
    const N5SizeMatrix* size_matrix,
    unsigned int row
)
{
    unsigned int sum;
    unsigned int column;

    sum = 0;

    for (column = 0; column < N5_V2_SIZE; ++column)
    {
        sum += size_matrix->matrix[row][column];
    }

    return sum;
}

/*
 * 计算 size matrix 指定一列的元素总和。
 */
static unsigned int n5_size_matrix_column_sum(
    const N5SizeMatrix* size_matrix,
    unsigned int column
)
{
    unsigned int sum;
    unsigned int row;

    sum = 0;

    for (row = 0; row < N5_V2_SIZE; ++row)
    {
        sum += size_matrix->matrix[row][column];
    }

    return sum;
}

/*
 * 计算指定 XOR diagonal 上的 cell size 总和。
 *
 * 数学定义：
 *
 *   sum_{a XOR b = w} k[a][b]。
 */
static unsigned int n5_size_matrix_xor_sum(
    const N5SizeMatrix* size_matrix,
    unsigned int w
)
{
    unsigned int sum;
    unsigned int a;
    unsigned int b;

    sum = 0;

    for (a = 0; a < N5_V2_SIZE; ++a)
    {
        for (b = 0; b < N5_V2_SIZE; ++b)
        {
            if ((a ^ b) == w)
            {
                sum += size_matrix->matrix[a][b];
            }
        }
    }

    return sum;
}

/*
 * 判断一个整数编码是否属于 V3。
 */
int n5_v3_is_valid(unsigned int value)
{
    return value < N5_V3_SIZE;
}

/*
 * 计算 V3 中两个向量的和。
 *
 * V3 = F_2^3，因此向量加法就是 XOR。
 */
N5V3 n5_v3_add(
    N5V3 a,
    N5V3 b
)
{
    return (N5V3)(a ^ b);
}

/*
 * 判断元素 x 是否属于给定的 V3 子集。
 */
int n5_subset_contains(
    N5SubsetMask subset,
    N5V3 x
)
{
    N5SubsetMask bit;

    bit = (N5SubsetMask)(1u << x);

    return (subset & bit) != 0;
}

/*
 * 计算 V3 子集的元素个数。
 */
unsigned int n5_subset_size(
    N5SubsetMask subset
)
{
    unsigned int count;
    unsigned int x;

    count = 0;

    for (x = 0; x < N5_V3_SIZE; ++x)
    {
        if (n5_subset_contains(
            subset,
            (N5V3)x
        ))
        {
            ++count;
        }
    }

    return count;
}

/*
 * 计算 V3 子集的 XOR charge。
 *
 *   sigma(S) = XOR_{x in S} x
 */
N5V3 n5_subset_xor_charge(
    N5SubsetMask subset
)
{
    N5V3 charge;
    unsigned int x;

    charge = 0;

    for (x = 0; x < N5_V3_SIZE; ++x)
    {
        if (n5_subset_contains(
            subset,
            (N5V3)x
        ))
        {
            charge = n5_v3_add(
                charge,
                (N5V3)x
            );
        }
    }

    return charge;
}

/*
 * 对 V3 子集执行 translation。
 *
 *   S + t = {x XOR t : x in S}
 */
N5SubsetMask n5_subset_translate(
    N5SubsetMask subset,
    N5V3 t
)
{
    N5SubsetMask translated;
    unsigned int x;

    translated = 0;

    for (x = 0; x < N5_V3_SIZE; ++x)
    {
        if (n5_subset_contains(
            subset,
            (N5V3)x
        ))
        {
            N5V3 translated_x;

            translated_x =
                n5_v3_add(
                    (N5V3)x,
                    t
                );

            translated |=
                (N5SubsetMask)(
                    1u << translated_x
                    );
        }
    }

    return translated;
}

/*
 * 判断一个 4×4 size matrix 是否满足全部合法性条件。
 */
int n5_size_matrix_is_legal(
    const N5SizeMatrix* size_matrix
)
{
    unsigned int row;
    unsigned int column;
    unsigned int w;

    if (size_matrix == NULL)
    {
        return 0;
    }

    for (row = 0; row < N5_V2_SIZE; ++row)
    {
        for (column = 0;
            column < N5_V2_SIZE;
            ++column)
        {
            if (size_matrix->matrix[row][column] >
                N5_CELL_MAX_SIZE)
            {
                return 0;
            }
        }
    }

    for (row = 0; row < N5_V2_SIZE; ++row)
    {
        if (n5_size_matrix_row_sum(
            size_matrix,
            row
        ) != N5_V3_SIZE)
        {
            return 0;
        }
    }

    for (column = 0;
        column < N5_V2_SIZE;
        ++column)
    {
        if (n5_size_matrix_column_sum(
            size_matrix,
            column
        ) != N5_V3_SIZE)
        {
            return 0;
        }
    }

    for (w = 0; w < N5_V2_SIZE; ++w)
    {
        if (n5_size_matrix_xor_sum(
            size_matrix,
            w
        ) != N5_V3_SIZE)
        {
            return 0;
        }
    }

    return 1;
}


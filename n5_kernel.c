#include <stdlib.h>

#include "n5_kernel.h"

/*
 * 判断指定 r 是否允许出现在给定 stabilizer 上。
 *
 * 数学条件：
 *
 *   r · h = 0
 *
 * 对 stabilizer 中全部 h 都成立。
 */
static int n5_fourier_label_allowed(
    N5V3 r,
    const N5V3* stabilizer,
    unsigned int stabilizer_size
)
{
    unsigned int index;

    for (index = 0;
         index < stabilizer_size;
         ++index)
    {
        if (n5_v3_dot(
                r,
                stabilizer[index]
            ) != 0)
        {
            return 0;
        }
    }

    return 1;
}

/*
 * 对 orbit mask 和对应 translation 做同步排序。
 *
 * 排序规则：
 *   subset mask 从小到大。
 */
static void n5_sort_orbit_support(
    N5SubsetMask* masks,
    N5V3* translations,
    unsigned int count
)
{
    unsigned int i;

    for (i = 1; i < count; ++i)
    {
        N5SubsetMask current_mask;
        N5V3 current_translation;

        unsigned int j;

        current_mask =
            masks[i];

        current_translation =
            translations[i];

        j = i;

        while (j > 0 &&
               current_mask < masks[j - 1])
        {
            masks[j] =
                masks[j - 1];

            translations[j] =
                translations[j - 1];

            --j;
        }

        masks[j] =
            current_mask;

        translations[j] =
            current_translation;
    }
}

/*
 * 计算 V3 上的二进制点积。
 *
 *   r · x = parity(popcount(r & x))。
 */
unsigned int n5_v3_dot(
    N5V3 r,
    N5V3 x
)
{
    unsigned int value;
    unsigned int parity;

    value =
        (unsigned int)(r & x);

    parity = 0;

    while (value != 0)
    {
        parity ^= value & 1u;
        value >>= 1;
    }

    return parity;
}

/*
 * 计算 Walsh character：
 *
 *   chi_r(x) = (-1)^(r·x)。
 */
int n5_walsh_character(
    N5V3 r,
    N5V3 x
)
{
    return
        n5_v3_dot(r, x) == 0
            ? 1
            : -1;
}

/*
 * 计算局部计数 T3(A,B,C)。
 *
 * 固定 A 中最小元素 a，
 * 枚举 b ∈ B，
 * 令 c = a XOR b，
 * 然后递归删除 a、b、c。
 */
unsigned int n5_t3(
    N5SubsetMask a_mask,
    N5SubsetMask b_mask,
    N5SubsetMask c_mask
)
{
    unsigned int size_a;
    unsigned int size_b;
    unsigned int size_c;
    unsigned int x;

    N5V3 a;

    unsigned int total;
    unsigned int b;

    x = (unsigned int)a_mask;
    x = x - ((x >> 1u) & 0x55u);
    x = (x & 0x33u) + ((x >> 2u) & 0x33u);
    size_a = (x + (x >> 4u)) & 0x0Fu;

    x = (unsigned int)b_mask;
    x = x - ((x >> 1u) & 0x55u);
    x = (x & 0x33u) + ((x >> 2u) & 0x33u);
    size_b = (x + (x >> 4u)) & 0x0Fu;

    x = (unsigned int)c_mask;
    x = x - ((x >> 1u) & 0x55u);
    x = (x & 0x33u) + ((x >> 2u) & 0x33u);
    size_c = (x + (x >> 4u)) & 0x0Fu;

    if (size_a != size_b ||
        size_a != size_c)
    {
        return 0;
    }

    if (a_mask == 0)
    {
        return 1;
    }

    /* Inline the former one-use n5_subset_first_element helper. */
    a = 0;
    while ((a_mask & (N5SubsetMask)(1u << a)) == 0u)
    {
        ++a;
    }

    total = 0;

    for (b = 0; b < N5_V3_SIZE; ++b)
    {
        N5SubsetMask b_bit;

        N5V3 c;
        N5SubsetMask c_bit;

        b_bit =
            (N5SubsetMask)(1u << b);

        if ((b_mask & b_bit) == 0)
        {
            continue;
        }

        c =
            (N5V3)(a ^ (N5V3)b);

        c_bit =
            (N5SubsetMask)(1u << c);

        if ((c_mask & c_bit) == 0)
        {
            continue;
        }

        total +=
            n5_t3(
                (N5SubsetMask)(
                    a_mask ^
                    (N5SubsetMask)(1u << a)
                ),
                (N5SubsetMask)(
                    b_mask ^ b_bit
                ),
                (N5SubsetMask)(
                    c_mask ^ c_bit
                )
            );
    }

    return total;
}

/*
 * 从定义生成固定 even cardinality k 的
 * canonical translation-Fourier channels。
 */
unsigned int n5_generate_even_fourier_channels(
    unsigned int k,
    N5EvenFourierChannel* output,
    unsigned int capacity
)
{
    unsigned int seen[256];

    unsigned int channel_count;
    unsigned int mask_value;

    if (k > N5_V3_SIZE ||
        (k & 1u) != 0)
    {
        return 0;
    }

    for (mask_value = 0;
         mask_value < 256;
         ++mask_value)
    {
        seen[mask_value] = 0;
    }

    channel_count = 0;

    for (mask_value = 0;
         mask_value < 256;
         ++mask_value)
    {
        N5SubsetMask representative;

        N5SubsetMask orbit_masks[8];
        N5V3 orbit_translations[8];

        N5V3 stabilizer[8];

        unsigned int orbit_size;
        unsigned int stabilizer_size;

        unsigned int t;
        unsigned int r;

        if (n5_subset_size(
                (N5SubsetMask)mask_value
            ) != k)
        {
            continue;
        }

        if (seen[mask_value] != 0)
        {
            continue;
        }

        representative =
            (N5SubsetMask)mask_value;

        orbit_size = 0;

        for (t = 0; t < N5_V3_SIZE; ++t)
        {
            N5SubsetMask translated;
            unsigned int index;
            int already_present;

            translated =
                n5_subset_translate(
                    representative,
                    (N5V3)t
                );

            already_present = 0;

            for (index = 0;
                 index < orbit_size;
                 ++index)
            {
                if (orbit_masks[index] ==
                    translated)
                {
                    already_present = 1;
                    break;
                }
            }

            if (!already_present)
            {
                orbit_masks[orbit_size] =
                    translated;

                orbit_translations[orbit_size] =
                    (N5V3)t;

                ++orbit_size;
            }

            seen[translated] = 1;
        }

        n5_sort_orbit_support(
            orbit_masks,
            orbit_translations,
            orbit_size
        );

        stabilizer_size = 0;

        for (t = 0; t < N5_V3_SIZE; ++t)
        {
            if (n5_subset_translate(
                    representative,
                    (N5V3)t
                ) == representative)
            {
                stabilizer[stabilizer_size] =
                    (N5V3)t;

                ++stabilizer_size;
            }
        }

        if (orbit_size * stabilizer_size !=
            N5_V3_SIZE)
        {
            return 0;
        }

        for (r = 0; r < N5_V3_SIZE; ++r)
        {
            unsigned int index;

            if (!n5_fourier_label_allowed(
                    (N5V3)r,
                    stabilizer,
                    stabilizer_size
                ))
            {
                continue;
            }

            if (output != NULL)
            {
                N5EvenFourierChannel* channel;

                if (channel_count >= capacity)
                {
                    return 0;
                }

                channel =
                    &output[channel_count];

                channel->k =
                    (uint8_t)k;

                channel->representative =
                    representative;

                channel->r =
                    (N5V3)r;

                channel->s =
                    n5_subset_xor_charge(
                        representative
                    );

                channel->orbit_size =
                    (uint8_t)orbit_size;

                channel->stabilizer_size =
                    (uint8_t)stabilizer_size;

                channel->norm =
                    (uint8_t)orbit_size;

                for (index = 0;
                     index < orbit_size;
                     ++index)
                {
                    channel->orbit_masks[index] =
                        orbit_masks[index];

                    channel->coefficients[index] =
                        (int8_t)n5_walsh_character(
                            (N5V3)r,
                            orbit_translations[index]
                        );
                }

                for (index = orbit_size;
                     index < 8;
                     ++index)
                {
                    channel->orbit_masks[index] = 0;
                    channel->coefficients[index] = 0;
                }

                for (index = 0;
                     index < stabilizer_size;
                     ++index)
                {
                    channel->stabilizer[index] =
                        stabilizer[index];
                }

                for (index = stabilizer_size;
                     index < 8;
                     ++index)
                {
                    channel->stabilizer[index] = 0;
                }
            }

            ++channel_count;
        }
    }

    return channel_count;
}

/*
 * 从 raw T3 定义直接计算一个 transformed even cell
 * tensor coefficient。
 *
 * 这里不使用任何 hand-written support rule。
 * coefficient 完全由：
 *
 *   Fourier coefficients
 *
 * 与：
 *
 *   raw T3(A,B,C)
 *
 * 的有限求和得到。
 */
int n5_even_cell_tensor_value(
    const N5EvenFourierChannel* channel_a,
    const N5EvenFourierChannel* channel_b,
    const N5EvenFourierChannel* channel_c
)
{
    int total;

    unsigned int index_a;
    unsigned int index_b;
    unsigned int index_c;

    if (channel_a == NULL ||
        channel_b == NULL ||
        channel_c == NULL)
    {
        return 0;
    }

    if (channel_a->k != channel_b->k ||
        channel_a->k != channel_c->k)
    {
        return 0;
    }

    total = 0;

    for (index_a = 0;
         index_a < channel_a->orbit_size;
         ++index_a)
    {
        for (index_b = 0;
             index_b < channel_b->orbit_size;
             ++index_b)
        {
            for (index_c = 0;
                 index_c < channel_c->orbit_size;
                 ++index_c)
            {
                int coefficient_product;
                unsigned int raw_t3;

                coefficient_product =
                    (int)channel_a->coefficients[index_a] *
                    (int)channel_b->coefficients[index_b] *
                    (int)channel_c->coefficients[index_c];

                raw_t3 =
                    n5_t3(
                        channel_a->orbit_masks[index_a],
                        channel_b->orbit_masks[index_b],
                        channel_c->orbit_masks[index_c]
                    );

                total +=
                    coefficient_product *
                    (int)raw_t3;
            }
        }
    }

    return total;
}

/*
 * 对一个 k=4 plane Fourier channel 求其 Fano normal。
 *
 * plane stabilizer 为 2D subspace，
 * 因此其正交补恰好包含：
 *
 *   {0,n}
 *
 * 其中 n 是唯一 nonzero normal。
 */
N5V3 n5_k4_plane_normal(
    const N5EvenFourierChannel* channel
)
{
    unsigned int n;

    if (channel == NULL ||
        channel->k != 4 ||
        channel->s != 0 ||
        channel->stabilizer_size != 4 ||
        channel->orbit_size != 2)
    {
        return 0;
    }

    for (n = 1; n < N5_V3_SIZE; ++n)
    {
        unsigned int index;
        int orthogonal;

        orthogonal = 1;

        for (index = 0;
             index < channel->stabilizer_size;
             ++index)
        {
            if (n5_v3_dot(
                    (N5V3)n,
                    channel->stabilizer[index]
                ) != 0)
            {
                orthogonal = 0;
                break;
            }
        }

        if (orthogonal)
        {
            return (N5V3)n;
        }
    }

    return 0;
}

/*
 * 根据四个 local channel indices 计算
 * generic partition tensor 的 flat index。
 *
 * 参数：
 *   tensor - partition tensor；
 *   q0..q3 - 四条腿 channel index。
 *
 * 返回：
 *   flat index。
 */
static unsigned int n5_even_partition_flat_index(
    const N5EvenPartitionTensor* tensor,
    unsigned int q0,
    unsigned int q1,
    unsigned int q2,
    unsigned int q3
)
{
    return
        (((q0 * tensor->channel_counts[1] + q1)
           * tensor->channel_counts[2] + q2)
          * tensor->channel_counts[3] + q3);
}

/*
 * 判断一个 cardinality 是否属于当前 all-even
 * channel library 支持集合。
 *
 * 返回：
 *   1 表示 k∈{0,2,4,6,8}；
 *   0 表示不支持。
 */
static int n5_even_cardinality_is_supported(
    unsigned int k
)
{
    return
        k <= N5_V3_SIZE &&
        (k & 1u) == 0;
}

/*
 * 从 ordered-partition 定义直接构造任意 all-even
 * 四腿 partition tensor。
 *
 * 实现严格按照：
 *
 *   Pi(q0,q1,q2,q3)
 *   =
 *   sum_{S0 sqcup S1 sqcup S2 sqcup S3=V3}
 *     product_i c_qi(Si)
 *
 * 执行。
 *
 * 不使用任何 complement/Fano/charge closed form
 * 来生成 tensor。
 */
int n5_build_even_partition_tensor(
    const unsigned int sizes[4],
    N5EvenPartitionTensor* tensor
)
{
    N5EvenFourierChannel channels[4][N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];

    uint8_t membership_count[4][256];
    uint8_t membership_channel[4][256][8];
    int8_t membership_coefficient[4][256][8];

    unsigned int leg;
    unsigned int mask_value;

    unsigned int size_sum;
    unsigned int coloring;

    if (sizes == NULL ||
        tensor == NULL)
    {
        return 0;
    }

    tensor->values = NULL;
    tensor->ordered_partitions = 0;
    tensor->entry_count = 0;

    size_sum = 0;

    for (leg = 0; leg < 4; ++leg)
    {
        if (!n5_even_cardinality_is_supported(
                sizes[leg]
            ))
        {
            return 0;
        }

        tensor->sizes[leg] =
            sizes[leg];

        size_sum +=
            sizes[leg];

        tensor->channel_counts[leg] =
            n5_generate_even_fourier_channels(
                sizes[leg],
                channels[leg],
                N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
            );

        if (tensor->channel_counts[leg] == 0)
        {
            return 0;
        }
    }

    if (size_sum != N5_V3_SIZE)
    {
        return 0;
    }

    tensor->entry_count = 1;

    for (leg = 0; leg < 4; ++leg)
    {
        tensor->entry_count *=
            tensor->channel_counts[leg];
    }

    tensor->values =
        (int16_t*)calloc(
            tensor->entry_count,
            sizeof(int16_t)
        );

    if (tensor->values == NULL)
    {
        tensor->entry_count = 0;
        return 0;
    }

    /*
     * 初始化 concrete subset -> channel membership。
     */
    for (leg = 0; leg < 4; ++leg)
    {
        for (mask_value = 0;
             mask_value < 256;
             ++mask_value)
        {
            membership_count[leg][mask_value] = 0;
        }
    }

    /*
     * 对每条腿分别建立 Fourier membership。
     */
    for (leg = 0; leg < 4; ++leg)
    {
        unsigned int channel_index;

        for (channel_index = 0;
             channel_index < tensor->channel_counts[leg];
             ++channel_index)
        {
            unsigned int support_index;

            for (support_index = 0;
                 support_index <
                    channels[leg][channel_index].orbit_size;
                 ++support_index)
            {
                N5SubsetMask subset;
                unsigned int slot;

                subset =
                    channels[leg][channel_index]
                        .orbit_masks[support_index];

                slot =
                    membership_count[leg][subset];

                if (slot >= 8)
                {
                    n5_free_even_partition_tensor(
                        tensor
                    );

                    return 0;
                }

                membership_channel[leg][subset][slot] =
                    (uint8_t)channel_index;

                membership_coefficient[leg][subset][slot] =
                    channels[leg][channel_index]
                        .coefficients[support_index];

                ++membership_count[leg][subset];
            }
        }
    }

    /*
     * 完全枚举 4^8 个 point-colorings。
     *
     * 每个 V3 point 被赋给 leg 0..3 中恰好一个。
     */
    for (coloring = 0;
         coloring < (1u << 16);
         ++coloring)
    {
        unsigned int code;

        N5SubsetMask subsets[4];
        unsigned int counts[4];

        unsigned int x;

        int cardinalities_match;

        code = coloring;

        for (leg = 0; leg < 4; ++leg)
        {
            subsets[leg] = 0;
            counts[leg] = 0;
        }

        for (x = 0; x < N5_V3_SIZE; ++x)
        {
            unsigned int color;

            color =
                code & 3u;

            code >>= 2;

            subsets[color] |=
                (N5SubsetMask)(1u << x);

            ++counts[color];
        }

        cardinalities_match = 1;

        for (leg = 0; leg < 4; ++leg)
        {
            if (counts[leg] != sizes[leg])
            {
                cardinalities_match = 0;
                break;
            }
        }

        if (!cardinalities_match)
        {
            continue;
        }

        ++tensor->ordered_partitions;

        /*
         * 对当前 concrete ordered partition，
         * 枚举四条腿 Fourier memberships 的
         * Cartesian product 并累加 exact coefficient。
         */
        {
            unsigned int m0;
            unsigned int m1;
            unsigned int m2;
            unsigned int m3;

            for (m0 = 0;
                 m0 < membership_count[0][subsets[0]];
                 ++m0)
            {
                unsigned int q0;
                int c0;

                q0 =
                    membership_channel[0][subsets[0]][m0];

                c0 =
                    membership_coefficient[0][subsets[0]][m0];

                for (m1 = 0;
                     m1 < membership_count[1][subsets[1]];
                     ++m1)
                {
                    unsigned int q1;
                    int c1;

                    q1 =
                        membership_channel[1][subsets[1]][m1];

                    c1 =
                        membership_coefficient[1][subsets[1]][m1];

                    for (m2 = 0;
                         m2 < membership_count[2][subsets[2]];
                         ++m2)
                    {
                        unsigned int q2;
                        int c2;

                        q2 =
                            membership_channel[2][subsets[2]][m2];

                        c2 =
                            membership_coefficient[2][subsets[2]][m2];

                        for (m3 = 0;
                             m3 < membership_count[3][subsets[3]];
                             ++m3)
                        {
                            unsigned int q3;
                            int c3;

                            unsigned int flat_index;

                            q3 =
                                membership_channel[3][subsets[3]][m3];

                            c3 =
                                membership_coefficient[3][subsets[3]][m3];

                            flat_index =
                                n5_even_partition_flat_index(
                                    tensor,
                                    q0,
                                    q1,
                                    q2,
                                    q3
                                );

                            tensor->values[flat_index] +=
                                (int16_t)(
                                    c0 * c1 * c2 * c3
                                );
                        }
                    }
                }
            }
        }
    }

    return 1;
}

/*
 * 释放 generic partition tensor 内部动态内存。
 */
void n5_free_even_partition_tensor(
    N5EvenPartitionTensor* tensor
)
{
    unsigned int leg;

    if (tensor == NULL)
    {
        return;
    }

    if (tensor->values != NULL)
    {
        free(tensor->values);
    }

    tensor->values = NULL;
    tensor->ordered_partitions = 0;
    tensor->entry_count = 0;

    for (leg = 0; leg < 4; ++leg)
    {
        tensor->sizes[leg] = 0;
        tensor->channel_counts[leg] = 0;
    }
}

/*
 * 读取 generic partition tensor 的一个 entry。
 */
int n5_even_partition_tensor_value(
    const N5EvenPartitionTensor* tensor,
    unsigned int q0,
    unsigned int q1,
    unsigned int q2,
    unsigned int q3
)
{
    unsigned int flat_index;

    if (tensor == NULL ||
        tensor->values == NULL)
    {
        return 0;
    }

    if (q0 >= tensor->channel_counts[0] ||
        q1 >= tensor->channel_counts[1] ||
        q2 >= tensor->channel_counts[2] ||
        q3 >= tensor->channel_counts[3])
    {
        return 0;
    }

    flat_index =
        n5_even_partition_flat_index(
            tensor,
            q0,
            q1,
            q2,
            q3
        );

    return
        (int)tensor->values[flat_index];
}


#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "n5_odd_kernel.h"

#define N5_ODD_FULL_MASK ((N5SubsetMask)0xFFu)

/*
 * 计算 subset mask 的 translation S -> S+t。
 */
static N5SubsetMask n5_odd_translate_mask(
    N5SubsetMask mask,
    N5V3 t
)
{
    N5SubsetMask result;
    unsigned int x;

    result = 0;

    for (x = 0; x < 8; ++x)
    {
        if ((mask & (N5SubsetMask)(1u << x)) != 0)
        {
            result = (N5SubsetMask)(
                result |
                (N5SubsetMask)(1u << (x ^ (unsigned int)t))
            );
        }
    }

    return result;
}

/*
 * 生成 Fano 3-subset S3(n)=n^perp\{0}。
 */
static N5SubsetMask n5_odd_fano3(
    unsigned int n
)
{
    N5SubsetMask mask;
    unsigned int x;

    mask = 0;

    if (n < 1 || n > 7)
    {
        return 0;
    }

    for (x = 1; x < 8; ++x)
    {
        if (n5_v3_dot((N5V3)n, (N5V3)x) == 0)
        {
            mask = (N5SubsetMask)(
                mask | (N5SubsetMask)(1u << x)
            );
        }
    }

    return mask;
}

/*
 * 返回 odd cardinality 的 zero-sum representative。
 */
static N5SubsetMask n5_odd_zero_sum_representative(
    unsigned int k,
    unsigned int fano_label
)
{
    if (k == 1)
    {
        return (N5SubsetMask)1u;
    }

    if (k == 7)
    {
        return (N5SubsetMask)(N5_ODD_FULL_MASK ^ 1u);
    }

    if (k == 3)
    {
        return n5_odd_fano3(fano_label);
    }

    if (k == 5)
    {
        return (N5SubsetMask)(
            N5_ODD_FULL_MASK ^ n5_odd_fano3(fano_label)
        );
    }

    return 0;
}

/*
 * 在统一 basis channel 中寻找指定 concrete subset 的系数。
 */
static int n5_odd_channel_coefficient(
    const N5OddBasisChannel* channel,
    N5SubsetMask mask
)
{
    unsigned int index;

    for (index = 0;
         index < channel->support_count;
         ++index)
    {
        if (channel->masks[index] == mask)
        {
            return channel->coefficients[index];
        }
    }

    return 0;
}

/*
 * 填充一个 odd translation-Fourier channel。
 */
static int n5_odd_fill_odd_channel(
    N5OddBasisChannel* channel,
    unsigned int k,
    N5V3 r,
    unsigned int fano_label
)
{
    N5SubsetMask base;
    unsigned int x;

    if (channel == NULL)
    {
        return 0;
    }

    base = n5_odd_zero_sum_representative(
        k,
        fano_label
    );

    if (base == 0 && k != 0)
    {
        return 0;
    }

    memset(channel, 0, sizeof(*channel));

    channel->k = (uint8_t)k;
    channel->r = r;
    channel->s = 0;
    channel->fano_label = (uint8_t)fano_label;
    channel->exceptional = 0;
    channel->support_count = 8;
    channel->norm = 8;

    for (x = 0; x < 8; ++x)
    {
        channel->masks[x] =
            n5_odd_translate_mask(
                base,
                (N5V3)x
            );

        channel->coefficients[x] =
            (int8_t)n5_walsh_character(
                r,
                (N5V3)x
            );
    }

    return 1;
}

/*
 * 将现有 all-even canonical channel 复制到统一 basis record。
 */
static int n5_odd_copy_even_channel(
    N5OddBasisChannel* target,
    const N5EvenFourierChannel* source
)
{
    unsigned int index;

    if (target == NULL || source == NULL)
    {
        return 0;
    }

    memset(target, 0, sizeof(*target));

    target->k = source->k;
    target->r = source->r;
    target->s = source->s;
    target->support_count = source->orbit_size;
    target->norm = source->norm;

    if (source->k == 4 &&
        source->s == 0 &&
        source->r == 0)
    {
        N5V3 normal;

        normal = n5_k4_plane_normal(
            source
        );

        if (normal != 0)
        {
            target->exceptional = 1;
            target->fano_label = normal;
        }
    }

    for (index = 0;
         index < source->orbit_size;
         ++index)
    {
        target->masks[index] =
            source->orbit_masks[index];

        target->coefficients[index] =
            source->coefficients[index];
    }

    return 1;
}

/*
 * 生成固定 cardinality 的 unified channel list。
 */
static unsigned int n5_odd_build_channels_for_k(
    unsigned int k,
    N5OddBasisChannel output[N5_ODD_KERNEL_MAX_CHANNELS_PER_K]
)
{
    unsigned int count;

    count = 0;

    if ((k & 1u) == 0)
    {
        N5EvenFourierChannel even[N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY];
        unsigned int even_count;
        unsigned int index;

        even_count = n5_generate_even_fourier_channels(
            k,
            even,
            N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY
        );

        if (even_count == 0 || even_count > N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY)
        {
            return 0;
        }

        for (index = 0; index < even_count; ++index)
        {
            if (!n5_odd_copy_even_channel(
                    &output[index],
                    &even[index]
                ))
            {
                return 0;
            }
        }

        return even_count;
    }

    if (k == 1 || k == 7)
    {
        unsigned int r;

        for (r = 0; r < 8; ++r)
        {
            if (!n5_odd_fill_odd_channel(
                    &output[count],
                    k,
                    (N5V3)r,
                    0
                ))
            {
                return 0;
            }

            ++count;
        }

        return count;
    }

    if (k == 3 || k == 5)
    {
        unsigned int r;
        unsigned int n;

        for (r = 0; r < 8; ++r)
        {
            for (n = 1; n <= 7; ++n)
            {
                if (!n5_odd_fill_odd_channel(
                        &output[count],
                        k,
                        (N5V3)r,
                        n
                    ))
                {
                    return 0;
                }

                ++count;
            }
        }

        return count;
    }

    return 0;
}

/*
 * 返回 flattened M3 table index。
 */
static unsigned int n5_odd_m3_index(
    unsigned int a,
    unsigned int b,
    unsigned int c
)
{
    return
        (a - 1u) * 49u +
        (b - 1u) * 7u +
        (c - 1u);
}

/*
 * 从当前 all-even basis provider 与 odd zero-sum representatives
 * 构造统一 KernelLibrary。
 */
int n5_odd_kernel_build_library(
    N5OddKernelLibrary* library
)
{
    unsigned int k;
    unsigned int a;
    unsigned int b;
    unsigned int c;

    if (library == NULL)
    {
        return 0;
    }

    memset(library, 0, sizeof(*library));

    for (k = 0; k <= 8; ++k)
    {
        library->channel_count[k] =
            n5_odd_build_channels_for_k(
                k,
                library->channels[k]
            );

        if (library->channel_count[k] == 0)
        {
            return 0;
        }
    }

    for (a = 1; a <= 7; ++a)
    {
        N5SubsetMask sa;

        sa = n5_odd_zero_sum_representative(3, a);

        for (b = 1; b <= 7; ++b)
        {
            N5SubsetMask sb;

            sb = n5_odd_zero_sum_representative(3, b);

            for (c = 1; c <= 7; ++c)
            {
                N5SubsetMask sc;
                unsigned int value;

                sc = n5_odd_zero_sum_representative(3, c);
                value = n5_t3(sa, sb, sc);

                if (value > 255)
                {
                    return 0;
                }

                library->m3[n5_odd_m3_index(a, b, c)] =
                    (uint8_t)value;
            }
        }
    }

    library->odd_scalar_k1 = n5_t3(
        n5_odd_zero_sum_representative(1, 0),
        n5_odd_zero_sum_representative(1, 0),
        n5_odd_zero_sum_representative(1, 0)
    );

    library->odd_scalar_k7 = n5_t3(
        n5_odd_zero_sum_representative(7, 0),
        n5_odd_zero_sum_representative(7, 0),
        n5_odd_zero_sum_representative(7, 0)
    );

    return 1;
}

/*
 * 返回 M3(a,b,c)。
 */
int n5_odd_kernel_m3(
    const N5OddKernelLibrary* library,
    unsigned int a,
    unsigned int b,
    unsigned int c
)
{
    if (library == NULL ||
        a < 1 || a > 7 ||
        b < 1 || b > 7 ||
        c < 1 || c > 7)
    {
        return 0;
    }

    return (int)library->m3[
        n5_odd_m3_index(a, b, c)
    ];
}

/*
 * 返回 M5=4-2*M3。
 */
int n5_odd_kernel_m5(
    const N5OddKernelLibrary* library,
    unsigned int a,
    unsigned int b,
    unsigned int c
)
{
    return 4 - 2 * n5_odd_kernel_m3(
        library,
        a,
        b,
        c
    );
}

/*
 * 计算 ordered four-leg partition line 的 exact unnormalized kernel。
 */
int n5_odd_kernel_line4(
    const N5OddBasisChannel* q0,
    const N5OddBasisChannel* q1,
    const N5OddBasisChannel* q2,
    const N5OddBasisChannel* q3
)
{
    unsigned int i0;
    unsigned int i1;
    unsigned int i2;
    int total;

    if (q0 == NULL || q1 == NULL || q2 == NULL || q3 == NULL)
    {
        return 0;
    }

    if ((unsigned int)q0->k +
        (unsigned int)q1->k +
        (unsigned int)q2->k +
        (unsigned int)q3->k != 8u)
    {
        return 0;
    }

    total = 0;

    for (i0 = 0; i0 < q0->support_count; ++i0)
    {
        N5SubsetMask m0;
        int c0;

        m0 = q0->masks[i0];
        c0 = q0->coefficients[i0];

        for (i1 = 0; i1 < q1->support_count; ++i1)
        {
            N5SubsetMask m1;
            N5SubsetMask u01;
            int c1;

            m1 = q1->masks[i1];

            if ((m0 & m1) != 0)
            {
                continue;
            }

            c1 = q1->coefficients[i1];
            u01 = (N5SubsetMask)(m0 | m1);

            for (i2 = 0; i2 < q2->support_count; ++i2)
            {
                N5SubsetMask m2;
                N5SubsetMask m3;
                int c2;
                int c3;

                m2 = q2->masks[i2];

                if ((u01 & m2) != 0)
                {
                    continue;
                }

                c2 = q2->coefficients[i2];
                m3 = (N5SubsetMask)(
                    N5_ODD_FULL_MASK ^
                    (N5SubsetMask)(u01 | m2)
                );

                c3 = n5_odd_channel_coefficient(q3, m3);

                if (c3 != 0)
                {
                    total += c0 * c1 * c2 * c3;
                }
            }
        }
    }

    return total;
}


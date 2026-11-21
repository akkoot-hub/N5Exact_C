#include <stdio.h>
#include <string.h>

#include "n5_crt.h"

#define N5_BIG_BASE UINT64_C(1000000000)

/*
 * GCC 14 reports a false-positive maybe-uninitialized diagnostic when this
 * bounded limb update is inlined into the CRT reconstruction caller.
 * Keep this small helper out of GCC's inlining path; other compilers use
 * their normal optimization behavior.
 */
#if defined(__GNUC__) && !defined(__clang__)
#define N5_NOINLINE __attribute__((noinline))
#else
#define N5_NOINLINE
#endif

/*
 * Single source of truth for the pairwise-coprime CRT prime basis.
 */
static const uint32_t N5_CRT_PRIMES[N5_CRT_PRIME_COUNT] =
{
    1000003u,
    1000033u,
    1000037u,
    1000039u,
    1000081u,
    1000099u
};

const uint32_t* n5_crt_prime_basis(void)
{
    return N5_CRT_PRIMES;
}

int n5_crt_get_primes(uint32_t* output, unsigned int count)
{
    if (output == NULL || count < N5_CRT_PRIME_COUNT)
    {
        return 0;
    }

    memcpy(output, n5_crt_prime_basis(), sizeof(N5_CRT_PRIMES));
    return 1;
}

/*
 * 规范化 used，保证 zero 使用 used=1, limb[0]=0。
 */
static void n5_big_normalize(
    N5BigUInt* value
)
{
    if (value == NULL)
    {
        return;
    }

    while (value->used > 1 &&
           value->limb[value->used - 1u] == 0)
    {
        --value->used;
    }

    if (value->used == 0)
    {
        value->used = 1;
        value->limb[0] = 0;
    }
}

/*
 * 从 uint64 初始化大整数。
 */
static void n5_big_from_u64(
    N5BigUInt* value,
    uint64_t input
)
{
    unsigned int index;

    memset(value, 0, sizeof(*value));
    index = 0;

    do
    {
        value->limb[index] = (uint32_t)(input % N5_BIG_BASE);
        input /= N5_BIG_BASE;
        ++index;
    }
    while (input != 0 && index < N5_CRT_MAX_LIMBS);

    value->used = index;
    n5_big_normalize(value);
}

/*
 * 计算 big integer modulo 一个 32-bit modulus。
 */
static uint32_t n5_big_mod_small(
    const N5BigUInt* value,
    uint32_t modulus
)
{
    uint64_t remainder;
    int index;

    remainder = 0;

    for (index = (int)value->used - 1;
         index >= 0;
         --index)
    {
        remainder =
            (remainder * N5_BIG_BASE + value->limb[index]) % modulus;
    }

    return (uint32_t)remainder;
}

/*
 * value *= multiplier。
 *
 * multiplier 是当前 CRT prime，约 10^6；
 * 64-bit intermediate 足够安全。
 */
static int n5_big_multiply_small(
    N5BigUInt* value,
    uint32_t multiplier
)
{
    uint64_t carry;
    unsigned int index;

    carry = 0;

    for (index = 0; index < value->used; ++index)
    {
        uint64_t current;

        current =
            (uint64_t)value->limb[index] * multiplier + carry;

        value->limb[index] = (uint32_t)(current % N5_BIG_BASE);
        carry = current / N5_BIG_BASE;
    }

    while (carry != 0)
    {
        if (value->used >= N5_CRT_MAX_LIMBS)
        {
            return 0;
        }

        value->limb[value->used] = (uint32_t)(carry % N5_BIG_BASE);
        carry /= N5_BIG_BASE;
        ++value->used;
    }

    n5_big_normalize(value);
    return 1;
}

/*
 * target += source * multiplier。
 *
 * source 与 target 可不同长度；multiplier < current prime。
 */
static N5_NOINLINE int n5_big_add_scaled(
    N5BigUInt* target,
    const N5BigUInt* source,
    uint32_t multiplier
)
{
    uint64_t carry;
    unsigned int index;
    unsigned int limit;
    unsigned int target_used;

    if (multiplier == 0)
    {
        return 1;
    }

    carry = 0;
    limit = source->used;
    target_used = target->used;

    if (target->used < limit)
    {
        target->used = limit;
    }

    for (index = 0; index < limit; ++index)
    {
        uint64_t current;

        current =
            (uint64_t)source->limb[index] * multiplier +
            (uint64_t)(index < target_used ? target->limb[index] : 0u) +
            carry;

        target->limb[index] = (uint32_t)(current % N5_BIG_BASE);
        carry = current / N5_BIG_BASE;
    }

    index = limit;

    while (carry != 0)
    {
        uint64_t current;

        if (index >= N5_CRT_MAX_LIMBS)
        {
            return 0;
        }

        current =
            (uint64_t)(index < target_used ? target->limb[index] : 0u) +
            carry;

        target->limb[index] = (uint32_t)(current % N5_BIG_BASE);
        carry = current / N5_BIG_BASE;

        if (index >= target->used)
        {
            target->used = index + 1u;
        }

        ++index;
    }

    n5_big_normalize(target);
    return 1;
}

/*
 * Extended-Euclid modular inverse。
 *
 * 返回：
 *   1 表示 inverse 存在；
 *   0 表示 gcd(a,m) != 1。
 */
static int n5_mod_inverse_u32(
    uint32_t a,
    uint32_t modulus,
    uint32_t* inverse
)
{
    int64_t old_r;
    int64_t r;
    int64_t old_t;
    int64_t t;

    if (inverse == NULL || modulus <= 1)
    {
        return 0;
    }

    old_r = (int64_t)modulus;
    r = (int64_t)(a % modulus);
    old_t = 0;
    t = 1;

    while (r != 0)
    {
        int64_t quotient;
        int64_t next_r;
        int64_t next_t;

        quotient = old_r / r;

        next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;

        next_t = old_t - quotient * t;
        old_t = t;
        t = next_t;
    }

    if (old_r != 1)
    {
        return 0;
    }

    old_t %= (int64_t)modulus;

    if (old_t < 0)
    {
        old_t += modulus;
    }

    *inverse = (uint32_t)old_t;
    return 1;
}

/*
 * 逐模扩张 CRT：
 *
 *   x' = x + M * t,
 *
 * 其中
 *
 *   t = (r - x) * M^{-1} mod p。
 */
int n5_crt_reconstruct(
    const uint32_t* residues,
    const uint32_t* moduli,
    unsigned int count,
    N5CRTResult* result
)
{
    N5BigUInt x = {{0u}, 1u};
    N5BigUInt product;
    unsigned int index;

    if (residues == NULL || moduli == NULL || result == NULL ||
        count == 0)
    {
        return 0;
    }

    n5_big_from_u64(&product, 1);

    for (index = 0; index < count; ++index)
    {
        uint32_t modulus;
        uint32_t x_mod;
        uint32_t product_mod;
        uint32_t inverse;
        uint32_t delta;
        uint32_t coefficient;

        modulus = moduli[index];

        if (modulus <= 1 || residues[index] >= modulus)
        {
            return 0;
        }

        x_mod = n5_big_mod_small(&x, modulus);
        product_mod = n5_big_mod_small(&product, modulus);

        if (!n5_mod_inverse_u32(
                product_mod,
                modulus,
                &inverse
            ))
        {
            return 0;
        }

        delta = (uint32_t)(
            ((uint64_t)residues[index] + modulus - x_mod) % modulus
        );

        coefficient = (uint32_t)(
            ((uint64_t)delta * inverse) % modulus
        );

        if (!n5_big_add_scaled(
                &x,
                &product,
                coefficient
            ) ||
            !n5_big_multiply_small(
                &product,
                modulus
            ))
        {
            return 0;
        }
    }

    result->value = x;
    result->modulus_product = product;

    return 1;
}

/*
 * base-10^9 limbs -> canonical decimal string。
 */
int n5_biguint_to_decimal(
    const N5BigUInt* value,
    char* buffer,
    unsigned int buffer_size
)
{
    int written;
    unsigned int offset;
    int index;

    if (value == NULL || buffer == NULL || buffer_size == 0 ||
        value->used == 0)
    {
        return 0;
    }

    written = snprintf(
        buffer,
        buffer_size,
        "%u",
        value->limb[value->used - 1u]
    );

    if (written < 0 || (unsigned int)written >= buffer_size)
    {
        return 0;
    }

    offset = (unsigned int)written;

    for (index = (int)value->used - 2;
         index >= 0;
         --index)
    {
        written = snprintf(
            buffer + offset,
            buffer_size - offset,
            "%09u",
            value->limb[index]
        );

        if (written < 0 ||
            (unsigned int)written >= buffer_size - offset)
        {
            return 0;
        }

        offset += (unsigned int)written;
    }

    return 1;
}

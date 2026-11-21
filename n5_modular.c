#include <stdlib.h>
#include <string.h>

#include "n5_modular.h"

/*
 * 把 signed int64 值 canonicalize 到 0..modulus-1。
 */
static uint32_t n5_modular_normalize_i64(
    int64_t value,
    uint32_t modulus
)
{
    int64_t remainder;

    remainder =
        value % (int64_t)modulus;

    if (remainder < 0)
    {
        remainder += modulus;
    }

    return (uint32_t)remainder;
}

/*
 * modular multiplication。
 *
 * 使用 uint64_t 保存中间乘积。
 * Stage 21 modulus 小于 2^32，故乘积不会超过 uint64_t。
 */
static uint32_t n5_modular_multiply(
    uint32_t a,
    uint32_t b,
    uint32_t modulus
)
{
    return
        (uint32_t)(
            ((uint64_t)a * (uint64_t)b) % modulus
        );
}

/*
 * 用 extended Euclid 求 modular inverse。
 *
 * 返回：
 *   1 表示 inverse 存在；
 *   0 表示 gcd(value,modulus) != 1。
 */
static int n5_modular_inverse(
    uint32_t value,
    uint32_t modulus,
    uint32_t* inverse
)
{
    int64_t old_r;
    int64_t r;
    int64_t old_t;
    int64_t t;

    if (inverse == NULL ||
        modulus <= 1)
    {
        return 0;
    }

    value %= modulus;

    if (value == 0)
    {
        return 0;
    }

    old_r = modulus;
    r = value;

    old_t = 0;
    t = 1;

    while (r != 0)
    {
        int64_t quotient;
        int64_t temporary;

        quotient = old_r / r;

        temporary = old_r - quotient * r;
        old_r = r;
        r = temporary;

        temporary = old_t - quotient * t;
        old_t = t;
        t = temporary;
    }

    if (old_r != 1)
    {
        return 0;
    }

    old_t %= modulus;

    if (old_t < 0)
    {
        old_t += modulus;
    }

    *inverse = (uint32_t)old_t;

    return 1;
}

/*
 * exact rational -> modular residue。
 */
int n5_modular_rational_to_residue(
    N5SmallRational value,
    uint32_t modulus,
    uint32_t* residue
)
{
    uint32_t numerator;
    uint32_t denominator;
    uint32_t inverse;

    if (residue == NULL ||
        modulus <= 1 ||
        value.denominator == 0)
    {
        return 0;
    }

    numerator =
        n5_modular_normalize_i64(
            value.numerator,
            modulus
        );

    denominator =
        n5_modular_normalize_i64(
            value.denominator,
            modulus
        );

    if (!n5_modular_inverse(
            denominator,
            modulus,
            &inverse
        ))
    {
        return 0;
    }

    *residue =
        n5_modular_multiply(
            numerator,
            inverse,
            modulus
        );

    return 1;
}


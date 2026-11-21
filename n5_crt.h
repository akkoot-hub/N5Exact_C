#ifndef N5_CRT_H
#define N5_CRT_H

#include <stdint.h>

#define N5_CRT_PRIME_COUNT 6u
#define N5_CRT_MAX_LIMBS 8
#define N5_CRT_DECIMAL_BUFFER 96

/*
 * 一个足够覆盖 N5 CRT modulus product 的无符号大整数。
 *
 * limb 使用 little-endian base 10^9：
 *
 *   value = sum limb[i] * (10^9)^i。
 *
 * Stage 24 只需要：
 *   - mod small prime；
 *   - multiply by small prime；
 *   - add M * small coefficient；
 *   - decimal serialization。
 *
 * 因而不引入完整 arbitrary-precision library。
 */
/*
 * Return the shared six-prime CRT basis used by every production sector.
 * count must be at least N5_CRT_PRIME_COUNT.
 */
const uint32_t* n5_crt_prime_basis(void);
int n5_crt_get_primes(uint32_t* output, unsigned int count);

typedef struct
{
    uint32_t limb[N5_CRT_MAX_LIMBS];
    unsigned int used;
} N5BigUInt;

/*
 * 多模 CRT reconstruction 的结果。
 *
 * value：唯一的 0 <= x < modulus_product 解。
 * modulus_product：全部 pairwise-coprime moduli 的乘积。
 */
typedef struct
{
    N5BigUInt value;
    N5BigUInt modulus_product;
} N5CRTResult;

/*
 * 用逐步 CRT / Garner-style update 重建非负 canonical residue。
 *
 * 输入：
 *   residues[i] mod moduli[i]。
 *
 * 要求：
 *   moduli 两两互素且 > 1。
 *
 * 返回：
 *   1 表示成功；
 *   0 表示输入非法、模逆不存在或固定 limb 容量不足。
 */
int n5_crt_reconstruct(
    const uint32_t* residues,
    const uint32_t* moduli,
    unsigned int count,
    N5CRTResult* result
);

/*
 * 把 N5BigUInt 写成十进制字符串。
 *
 * buffer 必须至少有 N5_CRT_DECIMAL_BUFFER bytes，
 * 或者调用者传入更大的 buffer。
 */
int n5_biguint_to_decimal(
    const N5BigUInt* value,
    char* buffer,
    unsigned int buffer_size
);

#endif

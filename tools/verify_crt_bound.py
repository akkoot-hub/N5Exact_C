#!/usr/bin/env python3
"""Verify that the production CRT modulus product exceeds the trivial 5-bit count bound.

This is an independent documentation/audit helper. It does not run the N5
enumerator and does not contain or compare an expected N5 answer.
"""
from __future__ import annotations

import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "n5_crt.h"
SOURCE = ROOT / "n5_crt.c"


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def main() -> int:
    header = HEADER.read_text(encoding="utf-8")
    source = SOURCE.read_text(encoding="utf-8")

    count_match = re.search(
        r"^\s*#define\s+N5_CRT_PRIME_COUNT\s+(\d+)u?\s*$",
        header,
        re.MULTILINE,
    )
    if not count_match:
        raise SystemExit("Could not find N5_CRT_PRIME_COUNT in n5_crt.h")
    expected_count = int(count_match.group(1))

    array_match = re.search(
        r"N5_CRT_PRIMES\s*\[\s*N5_CRT_PRIME_COUNT\s*\]\s*=\s*\{([^}]*)\}",
        source,
        re.DOTALL,
    )
    if not array_match:
        raise SystemExit("Could not find N5_CRT_PRIMES in n5_crt.c")

    primes = [int(x) for x in re.findall(r"(\d+)u", array_match.group(1))]
    if len(primes) != expected_count:
        raise SystemExit(
            f"CRT basis length mismatch: macro={expected_count}, values={len(primes)}"
        )
    if not all(is_prime(p) for p in primes):
        raise SystemExit("At least one CRT modulus is not prime")
    if len(set(primes)) != len(primes):
        raise SystemExit("CRT moduli are not distinct")
    for i, left in enumerate(primes):
        for right in primes[:i]:
            if math.gcd(left, right) != 1:
                raise SystemExit(f"CRT moduli are not pairwise coprime: {left}, {right}")

    modulus_product = math.prod(primes)
    permutation_upper_bound = math.factorial(1 << 5)

    print(f"CRT prime count: {len(primes)}")
    print(f"CRT moduli: {', '.join(map(str, primes))}")
    print(f"CRT modulus product M: {modulus_product}")
    print(f"5-bit domain size: {1 << 5}")
    print(f"Trivial permutation upper bound 32!: {permutation_upper_bound}")
    print(f"M > 32!: {modulus_product > permutation_upper_bound}")

    if modulus_product <= permutation_upper_bound:
        raise SystemExit("CRT modulus product does not exceed the trivial upper bound")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

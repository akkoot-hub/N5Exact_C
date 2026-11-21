#ifndef N5_MODULAR_H
#define N5_MODULAR_H

#include <stdint.h>

#include "n5_sparse.h"

/* Map an exact small rational to a residue when its denominator is invertible. */
int n5_modular_rational_to_residue(
    N5SmallRational value,
    uint32_t modulus,
    uint32_t* residue
);

#endif

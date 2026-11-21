# Mathematical parameters, algorithm settings, and implementation limits

This file distinguishes mathematical dimensions from algorithm tuning and storage limits. None of the items below is a pre-entered orbit count, orbit contribution, or final N5 answer.

## Mathematical dimensions

| Definition | Meaning |
|---|---|
| `N5_V2_SIZE = 4` | Cardinality of the V2 index set |
| `N5_V3_SIZE = 8` | Cardinality of the V3 index set |
| `N5_CELL_MAX_SIZE = 8` | Maximum allowed size of a cell in a size matrix |
| `N5_TOP_CELL_COUNT = N5_V2_SIZE * N5_V2_SIZE` | Number of cells in a top-level size matrix |
| `N5_CONTRACTION_CELL_COUNT` | Number of cell factors in the fixed-K network |
| `N5_CONTRACTION_LINE_COUNT = 3 * N5_V2_SIZE` | Three families of V2 partition lines |
| `N5_CONTRACTION_INCIDENCE_COUNT = 3 * N5_CONTRACTION_CELL_COUNT` | Cell/line incidences |
| `N5_CONTRACTION_TOTAL_FACTOR_COUNT` | Sum of cell, line, and incidence factor counts |
| `N5_TOP_GROUP_SIZE` | Formula-derived size of the generated top-level symmetry group |
| `N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY` | Formula-derived maximum number of even Fourier channels at one cardinality |
| `N5_KERNEL_TWO_SUBSET_CHANNEL_COUNT` | Formula-derived two-subset channel count |
| `N5_CONTRACTION_MAX_FIXED_FLOW_DOMAIN` | Formula-derived upper bound for a fixed-flow incidence domain |
| `N5_ODD_CHARGE_MAX_NULLITY = N5_V3_SIZE` | Structural maximum nullity dimension used by the odd charge reducer |
| `N5_ODD_NETWORK_MAX_DOMAIN` | Derived from the Fourier-channel capacity |

The source computes group/channel/domain capacities from the V2/V3 formulas where appropriate. For example, the even Fourier-channel capacity is expressed as a binomial coefficient rather than a pre-entered result count.

## CRT configuration

`N5_CRT_PRIME_COUNT` is the single shared count definition. The six pairwise-coprime prime moduli are stored once in `n5_crt.c` and shared by the all-even, odd, and final aggregation paths.

The production executable assembles the modulus product at runtime. It does not embed a precomputed N5 answer.

## Solver tuning

| Setting | Meaning |
|---|---|
| `N5_ODD_PRODUCTION_DEFAULT_TARGET_PEAK` | Default target peak used by the odd adaptive solver |
| `N5_ODD_PRODUCTION_DEFAULT_MAX_SLICES` | Default maximum number of structural slices |

These are solver configuration values. They influence the choice among exact contraction paths; they are not counts of solutions or orbits.

## Fixed storage and capacity limits

| Setting | Meaning |
|---|---|
| `N5_SPARSE_MAX_ARITY` | Maximum number of variables in one sparse factor |
| `N5_SPARSE_MAX_VARIABLES` | Maximum variable count represented by the sparse network |
| `N5_SPARSE_MAX_VALUE` | Maximum scalar value stored by the sparse-factor representation |
| `N5_SPARSE_NETWORK_MAX_FACTORS` | Maximum number of factors in the sparse network |
| `N5_DENSE_MAX_MODULI` | Capacity for modular lanes in dense contraction; the configured CRT basis uses fewer lanes than this capacity |
| `N5_SLICING_MAX_VARIABLES` | Maximum variables represented by the generic slicing structures |
| `N5_PRODUCTION_MAX_SLICE_VARIABLES` | Maximum variables handled by the production lambda-slicing representation |
| `N5_CRT_MAX_LIMBS` | Number of base-10^9 limbs in the fixed-size CRT integer |
| `N5_CRT_DECIMAL_BUFFER` | Output buffer capacity for decimal serialization |

These limits are implementation capacities, not mathematical answer values. The code has bounds checks and returns failure when a represented object exceeds a supported capacity; it must not silently truncate an object. The final release review should retain this failure behavior.

The CRT integer uses base `10^9`. Eight limbs provide up to 72 decimal digits of storage. The six-prime modulus product is 37 decimal digits, and the 96-byte decimal output buffer is larger than the maximum 72-digit limb representation plus its terminating NUL byte.


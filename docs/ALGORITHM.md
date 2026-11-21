# Algorithm overview

This document describes the production flow implemented by the source. It is an implementation overview, not a substitute for the mathematical definitions and proofs in the associated paper.

## 1. Size-matrix catalog

The code uses the finite spaces represented by `N5_V2_SIZE = 4` and `N5_V3_SIZE = 8`. A candidate size matrix is a `4 x 4` array whose entries are integer cell sizes.

The catalog generator enumerates six integer parameters `p, q, r, s, t, h`, each from 0 through 8, and constructs the matrix

\[
K =
\begin{pmatrix}
p+r+s+t+h-8 & h+p-q & 8+q-s-t-h-p & 8-r-h-p \\
8+q-r-t-h-p & 8-s-h-p & r+s+t+2h+p-q-8 & p \\
t+h-q & q & 8-r-t-h & r \\
8-s-t-h & s & t & h
\end{pmatrix}.
\]

A candidate is retained only if every entry lies in the permitted cell-size range and the legality predicate succeeds. The predicate checks the row sums, column sums, and XOR-diagonal sums. For each `w` in the V2 index set, the XOR-diagonal sum is

\[
\sum_{a \mathbin{\mathrm{XOR}} b = w} K_{a,b}.
\]

The required row, column, and XOR-diagonal totals are derived from the V3 dimension. These are mathematical admissibility conditions, not expected-result checks.

## 2. Symmetry and canonical orbits

The top-level cell-permutation group is generated from the V2 structure by the orbit module. The program applies the generated action to the legal size matrices, selects canonical representatives, and constructs the orbit catalog at runtime.

Orbit IDs in the terminal output are assigned by the generated canonical catalog order. They are not loaded from a frozen historical orbit table.

The catalog records each orbit's canonical matrix, member count, stabilizer size, and family classification. The family classification determines which production solver is used.

## 3. All-even production path

The all-even solver starts from a canonical all-even size matrix and performs the stages documented in `n5_all_even.h`:

1. Build the scalar-flow nullspace/reduced-flow catalog.
2. Quotient the flow collection by the implemented `GL(3,2)` and size-matrix stabilizer actions.
3. Build the second-Walsh representation.
4. Select and apply the required lambda slicing/quotient path.
5. Perform the dense modular contraction for each CRT modulus.
6. Aggregate weighted flow-orbit residues in deterministic class-index order.
7. Reconstruct the canonical CRT representative.

The implementation contains direct, symmetry-sliced, and generic-sliced contraction paths. The exact selection logic is implemented in `n5_all_even.c`, `n5_contraction.c`, `n5_production.c`, `n5_slicing.c`, and `n5_dense.c`.

## 4. Odd production path

The odd solver builds the odd kernel library and reduces the fixed-charge collection. It then processes the reduced charge classes through the runtime parallel-for scheduler.

Each class is passed to the adaptive contractor. Depending on the generated factor-graph structure and configured thresholds, the path can use ordinary contraction, the second-Walsh representation, structural slicing, or a rewritten variant. Class results are aggregated by class index, independently of worker completion order.

The implementation is primarily in `n5_odd_kernel.c`, `n5_odd_charge.c`, `n5_odd_network.c`, `n5_odd_adaptive.c`, `n5_odd_production.c`, and the shared sparse/dense modules.

## 5. Orbit weighting and CRT

For a canonical representative `K` with top-level orbit size `|O|`, the program computes the six residues of `C(K)` and multiplies each residue by `|O|` modulo its corresponding CRT modulus. Weighted residues are accumulated separately for the all-even and odd sectors and for the total.

The six pairwise-coprime CRT moduli are defined once in `n5_crt.c`; their count is defined once in `n5_crt.h`. The modulus product is constructed at runtime. The output includes:

- six residues for each `C(K)`;
- six residues for each weighted contribution `|O|C(K)`;
- aggregate all-even, odd, and total residue vectors;
- CRT-reconstructed decimal representatives and the modulus product.

**Mathematical exactness bound:** the target is a count of permutations of the 5-bit domain, which has `2^5 = 32` elements. Therefore the number of complete permutations is bounded above by the number of all bijections of this domain:

\[
0 \leq N_5 \leq 32!.
\]

The six CRT moduli in the source are pairwise-coprime primes. Their product is

\[
M = 1000292032458727685153601621373570283,
\]

whereas

\[
32! = 263130836933693530167218012160000000.
\]

Thus `32! < M`. Since each sector subtotal and each orbit-weighted contribution is a nonnegative count bounded by the total count, the canonical CRT representative is the actual integer count, not merely a residue modulo `M`, provided the mathematical orbit partition and contribution formulas are correct. The production program constructs `M` at runtime; these decimal values are recorded here only as a transparent proof of the CRT range bound.

The script `tools/verify_crt_bound.py` independently reads the prime basis from `n5_crt.c`, verifies primality and pairwise coprimality, and checks `M > 32!`. It does not inspect or compare the computed N5 answer.

## 6. Source map

| File | Role |
|---|---|
| `main.c` | CLI parsing, thread selection request, top-level invocation, final terminal summary |
| `n5_core.c/.h` | V2/V3 operations and size-matrix legality |
| `n5_catalog.c/.h` | Six-parameter candidate enumeration and runtime legal-matrix catalog |
| `n5_orbit.c/.h` | Top-level symmetry group and canonical orbit construction |
| `n5_all_even.c/.h` | All-even production solver |
| `n5_odd_kernel.c/.h` | Odd kernel library |
| `n5_odd_charge.c/.h` | Odd fixed-charge reduction |
| `n5_odd_network.c/.h` | Odd factor-network construction |
| `n5_odd_adaptive.c/.h` | Adaptive odd contraction |
| `n5_odd_production.c/.h` | Odd production aggregation and CRT residues |
| `n5_contraction.c/.h` | Flow and contraction structures |
| `n5_production.c/.h` | Lambda slicing and production contraction orchestration |
| `n5_slicing.c/.h` | Slicing plans and enumeration |
| `n5_sparse.c/.h` | Sparse factor/network representation |
| `n5_dense.c/.h` | Dense modular contraction |
| `n5_modular.c/.h` | Modular arithmetic helpers |
| `n5_crt.c/.h` | Shared CRT basis and integer reconstruction |
| `n5_parallel.c/.h` | Runtime thread detection, scheduling, and monotonic timing |
| `n5_final.c/.h` | Top-level orchestration, orbit weighting, output and aggregate CRT |


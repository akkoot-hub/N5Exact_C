# ComputeOnly scope

This project is the standalone production calculator. It is intentionally not packaged with the historical regression suite, prior result files, frozen orbit catalogs, checkpoints, or expected-answer comparisons.

## Runtime-generated objects

The executable generates the legal size-matrix collection and top-level symmetry-orbit catalog from the six-parameter definition on every run. The number of legal matrices, number of generated orbits, orbit-family distribution, and final contributions are runtime outputs.

## Excluded from production

The production source does not contain:
- a pre-entered top-level orbit count;
- a pre-entered legal size-matrix count;
- pre-entered all-even or odd orbit counts;
- frozen per-orbit contributions;
- a pre-entered final N5 total;
- a post-computation comparison against a known answer.

## Retained mathematical and implementation parameters

The program necessarily retains mathematical dimensions and finite implementation capacities, including the V2/V3 sizes, the six-prime CRT basis, sparse-factor arity/variable/factor limits, dense-modulus capacity, and decimal/big-integer storage capacity. Solver tuning parameters such as the odd adaptive target peak and maximum slice count are algorithm configuration, not known result values.

Memory-allocation checks, domain checks, arithmetic-capacity checks, and solver-failure handling remain. They prevent incomplete or invalid execution from being reported as a successful result; they do not compare results with historical answers.

## Independent review output

For every generated orbit, the program prints the canonical matrix, orbit and stabilizer sizes, solver diagnostics, the contribution `C(K)`, the orbit-weighted contribution `|O|C(K)`, and both six-prime residue vectors. It also prints aggregate residue vectors for the all-even subtotal, odd subtotal, and total.

The CRT reconstruction returns the canonical representative in the interval `[0, M)`, where `M` is the product of the six moduli. To interpret this representative as the unreduced integer count, the mathematical publication must separately establish an upper bound showing that the target count lies below `M`. The program does not encode a precomputed answer or independently prove that external bound.

# Publication audit report — 2026-11-20

## Baseline identity

- Supplied candidate: `N5Exact_C_ComputeOnly_ThreadConfig(1).zip`
- SHA-256: `db4f315c29adbeca09083dc1d0f7612bae7bb1e6b4bac696faded0a2a89f0e4f`
- Result: exact match with the ZIP digest specified by the supplied GitHub Publication Handoff.
- Original archive: 148 ZIP entries, 121 files, 20 C source files, 19 headers.

## Release-directory cleanup

Removed from the release candidate:

- `.vs/` Visual Studio cache, browse databases, IntelliSense PCH files, and solution state;
- nested `N5Exact/x64/Release/` object files, PDB, IOBJ/IPDB, logs, recipes, and TLOG files;
- root `x64/Release/` executable and PDB;
- `N5Exact.vcxproj.user`;
- development-only `THREAD_CONFIG.patch`;
- old `THREAD_CONFIG_CHANGELOG.md` and stale `THREAD_CONFIG_SHA256SUMS.txt`.

Added or updated:

- `CHANGELOG.md`;
- expanded `.gitignore` for nested Visual Studio, MSBuild, GCC, and Clang products;
- complete `N5Exact.vcxproj` header inventory (all 19 headers are listed);
- source/documentation-only `SHA256SUMS.txt`, generated after final cleanup.

## Final source-only candidate

- Files in the cleaned candidate: 59 (including the checksum manifest).
- C source files: 20.
- Header files: 19.
- Total uncompressed size: approximately 0.54 MB.
- No `.vs`, `x64`, `Debug`, `Release`, executable, PDB, object, user settings, or development patch files remain.
- `SHA256SUMS.txt` covers all other 58 release files and passes `sha256sum --check`.

## Source audit and code changes

1. Rechecked all production `.c` and `.h` files for known-answer macros and expected-result comparisons. No `N5_TOP_ORBIT_COUNT`, `N5_LEGAL_K_COUNT`, `EXPECTED`, `KNOWN_ANSWER`, or expected N5 total was found.
2. Centralized the CRT prime count as `N5_CRT_PRIME_COUNT` in `n5_crt.h`.
3. Centralized the six CRT moduli in `n5_crt.c`; all-even, odd-adaptive, odd-production, and final aggregation now use the shared basis. The duplicated sector-specific prime-count macros and redundant getter wrappers were removed.
4. Removed unused configuration macros found by whole-source reference counting: `N5_ODD_ADAPTIVE_DEFAULT_TARGET_PEAK`, `N5_ODD_ADAPTIVE_DEFAULT_MAX_SLICES`, `N5_V3_FULL_MASK`, and `N5_ODD_NETWORK_LINE_COUNT`.
5. Replaced the literal Fourier-channel capacity `70u` in `n5_odd_network.c` with the formula-derived `N5_KERNEL_MAX_CHANNELS_PER_CARDINALITY`.
6. Corrected the read-only multidimensional-domain function parameter in `n5_production.c` to a flat `const uint8_t*` view. This avoids four C17 pedantic pointer-to-array qualifier diagnostics without changing the indexing or accepted domains.
7. Hardened `n5_big_add_scaled()` so limbs beyond the target's prior `used` count are treated as zero rather than read. A GCC-only no-inline annotation prevents GCC 14's false-positive `-Wmaybe-uninitialized` diagnostic in this helper. The CRT arithmetic round-trip smoke test passes with both GCC and Clang.
8. Corrected thread terminology in CLI output and headers: requested thread count means total participating threads including the main/calling thread.
9. The Windows launcher normalizes `AUTO` to lowercase before passing it to the case-sensitive CLI parser.
10. `BUILD_RELEASE.ps1` now prefers a Visual Studio 2026 (18.x) installation with MSBuild and the x86/x64 C++ tools component, rather than accidentally selecting an older Visual Studio MSBuild found first on PATH.
11. Added per-orbit `C(K)` and `|O|C(K)` six-residue vectors and aggregate all-even, odd, and total residue vectors to terminal output. This is output-only and does not change the computed residues.
12. Added algorithm, constant/capacity, reproducibility, publication, license-selection, citation-template, and Code Availability documentation.

## Build and smoke-test evidence

Current Linux review environment:

- GCC 14.2.0
- Clang 17.0.0
- CPU affinity available to this process: 5 logical processors

Results:

| Check | Result |
|---|---|
| GCC C17 Release with `-Wall -Wextra -Wpedantic -Werror` | PASS |
| Clang C17 Release with `-Wall -Wextra -Wpedantic -Werror` | PASS |
| CLI `--help` | PASS |
| CLI `--threads auto --help` | PASS |
| CLI `--threads=auto --help` | PASS |
| CLI `-j 1 --help` | PASS |
| Reject zero threads | PASS |
| Reject request above detected available processors | PASS |
| Reject duplicate thread options | PASS |
| Reject unknown option | PASS |
| Parallel scheduler exact-once smoke test for 1, 2, 4, auto, and task-limited workers | PASS |
| Linux one-CPU affinity smoke test | PASS |
| Generic CRT residue/reconstruction round trip | PASS |
| CRT prime primality and pairwise-coprime audit | PASS |
| CRT modulus product greater than `32!` | PASS |
| Visual Studio project and filters XML parsing | PASS |
| Visual Studio project lists all 20 C files and 19 headers | PASS |

No full N5 production run was performed in this audit. No Windows/MSVC compiler was available in this environment, so the VS2026 Release x64 build remains to be run on the target Windows installation.

## CRT exactness range

The input domain has 32 elements. The number of complete permutations is therefore bounded by `32!`. The shared CRT basis product exceeds `32!`:

- `M = 1000292032458727685153601621373570283`
- `32! = 263130836933693530167218012160000000`
- `M > 32!`

This establishes the CRT range bound independently of the unknown N5 answer. The script `tools/verify_crt_bound.py` checks the prime basis and inequality directly from the source. It does not calculate or compare an expected N5 result.

## Publication blockers requiring project-owner input

The code and source-only release candidate are prepared, but the package is not yet a complete public release because the supplied materials do not establish:

- copyright holder and permission to publish;
- selected license;
- GitHub owner/repository URL;
- final paper title and author order;
- ORCID identifiers, if any;
- paper publication status / preprint URL;
- software version and release date;
- final `CITATION.cff`;
- actual Windows/MSVC build and complete production-run record tied to the frozen source SHA-256.

`LICENSE` is intentionally not fabricated. A citation template and Code Availability template are provided, but placeholders must be replaced before release.


## Windows execution evidence received after the initial audit

A user-supplied terminal log records a successful Windows build and two complete production runs.

Environment reported by the log:

- Visual Studio Developer PowerShell 18.10.3
- MSBuild 18.10.1
- Visual Studio 2026 / MSVC toolset v145
- ISO C17, Release, x64
- 16 detected logical processors
- CPU model, physical core count, and installed RAM were not included in the log

Build result:

- `BUILD_RELEASE.ps1` completed successfully and produced `x64/Release/N5Exact.exe`.
- MSVC emitted C4819 warnings because UTF-8 source/header text was interpreted under the active Windows code page (936).
- The project now explicitly passes `/utf-8` to MSVC in all four configuration/platform combinations. This is a source-encoding/build-setting correction only; no C source or mathematical algorithm was changed.
- The new `/utf-8` setting still requires one local rebuild to confirm that C4819 warnings disappear.

Two complete runs were captured from the same local working directory:

| Run | Command | Participating threads | Wall time |
|---|---|---:|---:|
| Single-thread | `N5Exact.exe --threads 1` | 1 | 973.164 s (16.219 min) |
| Automatic | `N5Exact.exe --threads auto` | up to 16 | 180.337 s (3.006 min) |

The automatic run was approximately 5.396 times faster than the single-thread run.

A line-by-line comparison of the supplied log shows all 68 per-orbit `C(K)` values, six-prime residue vectors, weighted `|O|C(K)` values, and weighted residue vectors are identical between the two runs. The aggregate all-even, odd, and total residues, reconstructed subtotals, and final total are also identical.

The two runs both report:

- generated orbit catalog: 68;
- all-even / odd partition: 11 / 57;
- all-even subtotal: `23804600145881119850496`;
- odd subtotal: `911568826006449804017664`;
- total: `935373426152330923868160`;
- CRT modulus product: `1000292032458727685153601621373570283`.

The catalog statistics printed under “Definition-derived prerequisites” are runtime values, not hard-coded expected-answer gates: `n5_catalog.c` enumerates the six parameters over their definition ranges, dynamically counts legal matrices, generates the top-level group, constructs canonical orbits, and counts the resulting families. The production source does not contain an expected N5 total comparison. Structural constants such as the V2/V3 dimensions and formula-derived group allocation remain algorithm definitions, not stored answers.

Important provenance limitation: the supplied terminal log does not contain a source-manifest hash or Git commit. Therefore, the two runs are strong execution evidence for the user's local working tree, but are not yet cryptographically bound to the exact ZIP archive. Before the final public Release, verify the source manifest in the Windows working directory and record the final archive/commit hash.

The supplied log ends with successful process completion and the two computed result summaries. This is a complete production execution record, not a proof by itself that the mathematical model and implementation are correct; independent mathematical and historical-result review remain separate.


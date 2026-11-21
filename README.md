# N5Exact_C

**N5Exact_C** is a standalone C17 program for computing the exact number of 5-bit complete permutations from the mathematical definition.

The program rebuilds the size-matrix catalog and symmetry orbits at runtime, evaluates each generated orbit through the all-even or odd production path, combines modular residues, and prints the resulting integer counts.

> **Release status:** N5Exact_C v1.0.0 release candidate. The release metadata identifies Kai Wang as the individual copyright holder and software author, and the repository is designated as [akkoot-hub/N5Exact_C](https://github.com/akkoot-hub/N5Exact_C). Publication under MIT is based on the project owner's instruction to proceed on the assumption of individual ownership. Before making the repository public, the owner should confirm that no institutional, funder, employer, co-author, or other prior agreement restricts publication.

## What the program computes

At startup, the program:

1. Enumerates the six integer parameters used by the size-matrix parametrization.
2. Constructs candidate size matrices and retains those satisfying the mathematical legality predicate.
3. Generates the top-level symmetry action and constructs the canonical orbit catalog dynamically.
4. Dispatches each generated orbit to the all-even or odd production solver.
5. Multiplies each canonical representative's contribution by its top-level orbit size.
6. Accumulates the six CRT residues and reconstructs the all-even subtotal, odd subtotal, and total.
7. Prints the generated catalog, per-orbit contributions, per-orbit residue vectors, aggregate residue vectors, subtotals, total, CRT modulus product, and elapsed time.

No expected orbit count, legal-matrix count, sector count, per-orbit contribution, or final N5 answer is compiled into the production program. The executable does not compare its computed result with a stored answer.

The mathematical legality predicates and runtime safety checks remain part of the implementation. They are not expected-answer regression checks.

## Requirements

### Windows

- Windows 10 or later (Windows 11 recommended).
- Visual Studio 2026 (18.x) with the MSVC v145 C/C++ toolset.
- MSBuild.
- ISO C17 mode, x64 Release configuration.

### Linux

- GCC or Clang with C17 support.
- POSIX threads (`pthread`).
- Linux builds use `sched_getaffinity()` where available to detect the CPUs available to the current process.

The source uses the C standard library and operating-system APIs for threading and monotonic timing. No third-party C library is required by the current source.

Linux GCC and Clang builds are included in the reproducibility instructions. Windows/MSVC must be validated on the target Windows installation before the release is tagged.

## Source encoding

The Visual Studio project passes `/utf-8` to MSVC in all configurations so UTF-8 source comments and headers are read consistently on Windows systems using a non-UTF-8 active code page.

## Build on Windows

From the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File .\BUILD_RELEASE.ps1
```

The script builds `Release|x64` using the Visual Studio project.

## Run on Windows

```powershell
# Automatic thread selection (default)
.\x64\Release\N5Exact.exe

# Explicit total participating-thread count
.\x64\Release\N5Exact.exe --threads 8

# Single-thread execution
.\x64\Release\N5Exact.exe --threads 1

# Help
.\x64\Release\N5Exact.exe --help
```

The PowerShell launcher can build and run in one step:

```powershell
powershell -ExecutionPolicy Bypass -File .\RUN_FULL_WINDOWS.ps1
powershell -ExecutionPolicy Bypass -File .\RUN_FULL_WINDOWS.ps1 -Threads 8
```

## Build on Linux

GCC:

```bash
gcc -std=c17 -O3 -pthread -I. -Wall -Wextra -Wpedantic -Werror \
    n5_*.c main.c -o N5Exact
```

Clang:

```bash
clang -std=c17 -O3 -pthread -I. -Wall -Wextra -Wpedantic -Werror \
    n5_*.c main.c -o N5Exact
```

## Run on Linux

```bash
# Automatic thread selection
./N5Exact --threads auto

# Explicit total participating-thread count
./N5Exact --threads 8

# Single-thread execution
./N5Exact --threads 1

# Help
./N5Exact --help
```

The CLI accepts `--threads N`, `--threads=N`, `-j N`, and `--threads auto`. Omitting the option is equivalent to automatic selection.

## Runtime thread selection

Thread count is selected at runtime; it is not fixed in the source or compiler settings.

- `auto` is the default and uses the detected available logical-processor count.
- An explicit value must be between 1 and the detected available logical-processor count.
- Linux prefers the current process CPU-affinity mask, respecting affinity/cpuset restrictions when exposed by the OS.
- Windows queries active logical processors across processor groups.
- Other POSIX systems fall back to `_SC_NPROCESSORS_ONLN`.
- A parallel stage uses no more threads than the number of independent tasks in that stage.
- The main thread participates in the work. Thus `--threads 8` means up to eight total participating threads, not eight additional background threads.

Automatic selection uses all detected available logical processors. Select a smaller number if the machine must remain responsive for other workloads.

## CRT basis

The all-even and odd production paths use one shared CRT basis defined in `n5_crt.c`:

```text
1000003
1000033
1000037
1000039
1000081
1000099
```

The prime count is defined once in `n5_crt.h`. The CRT modulus product is assembled at runtime from this basis; it is not stored as an expected answer.

The program prints each orbit's six residues and weighted six residues, as well as the aggregate all-even, odd, and total residue vectors. These fields support independent review of the CRT reconstruction.

## Output

A successful run prints:

- the number of six-parameter candidates examined;
- the number of legal size matrices;
- the generated top-level symmetry group size;
- the number and family distribution of generated canonical orbits;
- each canonical size matrix and its orbit/stabilizer sizes;
- per-orbit solver diagnostics, exact contribution, weighted contribution, and six-lane residue vectors;
- all-even and odd subtotals;
- the final N5 total;
- the CRT modulus product and elapsed time.

The program exits with a nonzero status if the production computation cannot complete. A partial run must not be treated as a completed result.

## Source integrity

`SHA256SUMS.txt` contains SHA-256 checksums for the release files (excluding the checksum manifest itself). Verify them from the repository root with:

```bash
sha256sum --check SHA256SUMS.txt
```

## Documentation

- [Algorithm overview](docs/ALGORITHM.md)
- [Mathematical parameters and implementation limits](docs/CONSTANTS_AND_LIMITS.md)
- [Reproducibility record template](docs/REPRODUCIBILITY.md)
- [Publication checklist](docs/PUBLISHING_CHECKLIST.md)
- [Final publication audit report](docs/PUBLICATION_AUDIT_20261120.md)
- [Code Availability statement template](docs/CODE_AVAILABILITY_TEMPLATE.md)
- [ComputeOnly scope notes](COMPUTE_ONLY_NOTES.md)

## License and citation

This release candidate uses the **MIT License** and identifies **Kai Wang** as the copyright holder. The repository URL is <https://github.com/akkoot-hub/N5Exact_C>. The root-level `CITATION.cff` identifies Kai Wang as the software author and includes the supplied ORCID: <https://orcid.org/0009-0006-9441-0466>.

The version is set to `1.0.0`. The release date is intentionally omitted until the GitHub Release is actually created. Do not cite a placeholder DOI. Once the tagged GitHub Release is archived through Zenodo, cite the DOI assigned to that specific archived version. Review [LICENSE_SELECTION.md](docs/LICENSE_SELECTION.md) for the ownership assumption and publication caveat.

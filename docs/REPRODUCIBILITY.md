# Reproducibility record

Complete this record for the exact source snapshot that is used for a paper or release. Do not copy historical timings from another source version.

## Source identity

- Git repository:
- Git commit:
- Release tag:
- Source archive SHA-256:
- `SHA256SUMS.txt` SHA-256:

## Operating system and hardware

- OS name and version:
- CPU model:
- Physical core count:
- Logical processor count:
- Available logical processors/CPU affinity:
- Installed RAM:
- Memory limit, if running in a container or job scheduler:

## Toolchain

- Compiler:
- Compiler version:
- C language standard: ISO C17
- Build configuration: Release / optimized
- Exact compiler command or MSBuild command:
- Link/runtime dependencies:

## Runtime configuration

- Command line:
- Thread selection: `auto` or explicit integer
- Requested participating threads:
- Actual per-stage participating threads, where reported:
- Environment restrictions (CPU affinity, cpuset, job limits):

## Execution

- Start timestamp and timezone:
- End timestamp and timezone:
- Total wall-clock time:
- Exit code:
- Was the complete production run finished?:
- Were there any warnings or runtime errors?:

## Captured output

Linux example:

```bash
./N5Exact --threads auto 2>&1 | tee N5Exact-run.txt
```

Windows PowerShell example:

```powershell
powershell -ExecutionPolicy Bypass -File .\\RUN_FULL_WINDOWS.ps1 -Threads auto 2>&1 |
    Tee-Object -FilePath .\\N5Exact-run.txt
```

Attach the complete unedited terminal output. It should include:

- generated catalog and family distribution;
- each generated orbit and canonical matrix;
- each per-orbit `C(K)` and weighted `|O|C(K)` value;
- each per-orbit six-prime residue vector and weighted residue vector;
- aggregate all-even, odd, and total residue vectors;
- all-even subtotal, odd subtotal, total, and CRT modulus product;
- elapsed time.

## Cross-thread reproducibility

When comparing two runs, record the source commit and all environment details for both runs. Thread count and runtime may differ. The generated catalog, per-orbit exact contributions, residue vectors, subtotals, and total must agree exactly.

## CRT exactness note

The CRT output is the canonical representative modulo the product of the six moduli. Record the external mathematical bound or theorem establishing that the target count is smaller than this modulus product. Do not infer this condition merely from the fact that CRT reconstruction succeeded.

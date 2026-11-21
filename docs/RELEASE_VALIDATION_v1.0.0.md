# Windows validation for N5Exact_C v1.0.0

The accepted source commit is `a75a80d1d25b9a73f021a275f2259f42b8db5779`. Its manifest covers 60 Git blobs. All 60 checks passed. To retrieve the manifest used for acceptance, run `git show a75a80d1d25b9a73f021a275f2259f42b8db5779:SHA256SUMS.txt`. The release documentation commit updates the manifest to cover the validation report and revised README.

## Evidence fingerprints

| Acceptance artifact | SHA-256 |
| --- | --- |
| SHA256SUMS.txt at the accepted commit | `8e9a0c4860e5b6c734528773fb87721903dd8fed4b17314c963873b71cf6232c` |
| N5Exact.exe | `d7b117b2163a19232acb181d6a2472442c192a225b1b5fadaf83eeb7ec319424` |
| build.log | `f8ac4d82e3ebe8df22aa07a91f93de309efe9af1838fa6d7d3b94b5dfbf9e2ed` |
| full-run.log | `8ef21164e83bcf184a089ce2d8ffb901b26fd801bba47559079dc1a95af3fa84` |

## Build and run

On 2026-11-21, the Windows 11 x64 Release build completed with MSBuild 18.10.1, MSVC v145 and C17. The build log contains no C4819 warning.

```powershell
powershell -ExecutionPolicy Bypass -File .\BUILD_RELEASE.ps1
powershell -ExecutionPolicy Bypass -File .\RUN_FULL_WINDOWS.ps1 -NoBuild -Threads auto
```

The program detected 16 logical processors and exited with status 0. Program wall time was 125.817 seconds. The launcher took 126.311 seconds.

## Independent comparison

All 68 canonical matrices matched the historical Windows Step10 run. For every matrix, both the exact contribution and the orbit-size-weighted contribution matched. Reducing each historical integer modulo the six CRT primes reproduced both residue vectors for every orbit. The partition and total residue vectors also passed verification.

| Sector | Generated orbits | Reconstructed count |
| --- | ---: | ---: |
| All-even | 11 | 23804600145881119850496 |
| Odd | 57 | 911568826006449804017664 |
| Total | 68 | 935373426152330923868160 |

The project owner independently reviewed the acceptance package and confirmed the build, full computation, per-orbit agreement and source-manifest verification.

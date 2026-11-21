# Changelog

## v1.0.0 (release pending)

- Provide a standalone C17 build with MSVC v145 and UTF-8 source encoding.

- Rebuild the legal size-matrix catalog and top-level orbit catalog from the mathematical definition at runtime.
- Keep the all-even and odd production solvers and the six-prime CRT reconstruction.
- Add runtime thread selection (`auto`, `--threads N`, `--threads=N`, and `-j N`).
- Detect available logical processors at runtime; respect Linux process affinity/cpuset restrictions where available.
- Centralize the CRT prime basis and prime-count definition.
- Print per-orbit and aggregate six-lane CRT residue vectors for independent review.
- Complete the Visual Studio project header inventory and expand build-artifact ignore rules.
- Keep expected-answer constants and result-comparison gates out of the production executable.

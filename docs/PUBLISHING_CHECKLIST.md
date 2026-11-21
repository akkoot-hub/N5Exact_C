# Publication checklist

This is a release-candidate checklist. Items marked **owner confirmation required** cannot be completed from source code alone.

## Source and build

- [x] Candidate ZIP identity matched the supplied handoff SHA-256.
- [x] Remove Visual Studio caches, generated binaries, PDBs, object files, and local user settings from the source release directory.
- [x] Remove the development-only thread patch.
- [x] Extend `.gitignore` for nested build products.
- [x] Add missing header entries to the Visual Studio project.
- [x] Use one shared CRT prime-count definition and one shared CRT basis.
- [x] Print per-orbit and aggregate CRT residue vectors.
- [x] Generate and verify `SHA256SUMS.txt` for all release files (manifest excludes itself).
- [x] GCC C17 strict build (`-Wall -Wextra -Wpedantic -Werror`) passes in the current Linux review environment.
- [x] Clang C17 strict build (`-Wall -Wextra -Wpedantic -Werror`) passes in the current Linux review environment.
- [x] Windows Visual Studio 2026 / MSVC v145 Release x64 build completed on the target Windows machine (user-supplied log); initial build emitted C4819 encoding warnings.
- [x] Two complete production calculations completed: `--threads 1` and `--threads auto`; all per-orbit values and aggregate residues/results match exactly.
- [x] Preserve complete unedited production output (user-supplied terminal log).
- [ ] Rebuild after adding `/utf-8` to the Visual Studio project and confirm C4819 warnings are gone.
- [ ] Bind the complete run log to the final frozen source SHA-256 / Git commit.

## Mathematical documentation

- [x] Document the runtime size-matrix and orbit-generation flow.
- [x] Document the all-even and odd production paths at implementation level.
- [x] Document the shared CRT basis and residue output.
- [x] Include the trivial permutation bound `N5 <= 32!` and show that it is below the six-prime CRT modulus product.
- [x] Add `tools/verify_crt_bound.py` to independently verify the CRT basis and range inequality.

## Publication metadata — owner confirmation required

- [x] GitHub owner and repository name: `akkoot-hub/N5Exact_C` (repository creation/public visibility still pending).
- [x] Copyright holder recorded as Kai Wang, based on the owner's instruction to proceed on the individual-ownership assumption; owner must still confirm there are no conflicting institutional, funder, employer, or co-author rights.
- [x] Selected MIT License and populated root `LICENSE` file.
- [ ] Paper title.
- [x] Software author recorded as Kai Wang; final paper author order remains to be confirmed separately.
- [x] ORCID recorded: `0009-0006-9441-0466`.
- [ ] Publication status and paper URL/preprint DOI, if any.
- [x] Software version set to `1.0.0`; release date intentionally deferred until the GitHub Release is created.
- [x] Root `CITATION.cff` populated with author, ORCID, version, and repository URL.
- [ ] Replace the DOI placeholder in `CODE_AVAILABILITY_TEMPLATE.md` after Zenodo archive creation.

## Archive and DOI

- [ ] Create a clean initial Git commit.
- [ ] Inspect `git status` and the complete tracked-file list.
- [ ] Create a version tag and GitHub Release.
- [ ] Archive that exact Release with Zenodo.
- [ ] Obtain the actual DOI and update the paper's Code Availability statement.


## Owner metadata recorded for the release candidate

- Copyright holder: `Kai Wang` (individual-ownership assumption, as instructed by the project owner).
- GitHub repository: `https://github.com/akkoot-hub/N5Exact_C`.
- Software author: `Kai Wang`.
- ORCID: `https://orcid.org/0009-0006-9441-0466`.
- License: MIT.
- Software version: `1.0.0`.
- Release date: intentionally not set until the GitHub Release exists.

Before making the repository public, the owner should verify the individual-ownership assumption against any institutional, funder, employer, or co-author agreements. The repository itself, GitHub Release, Zenodo archive, DOI, final paper metadata, and DOI-backed Code Availability statement remain to be created or confirmed.

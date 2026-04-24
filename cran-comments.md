# smfa

## Resubmission

This is a resubmission addressing comments from the CRAN team.

### Major Change: Package Renaming
- **Package Name**: Renamed the package from `metafrontieR` to `smfa` to resolve the conflict with the existing CRAN package `metafrontier`. All internal references, documentation, and vignettes have been updated accordingly.

### Other Changes:
- Removed `\dontrun{}` from all fast-executing examples in `smfa.Rd`.
- Added a lightweight toy example for the Sample Selection model to allow for automatic CRAN testing of this feature within the 5-second limit.
- Enabled code execution in all vignettes by changing `eval = FALSE` to `eval = TRUE` in the setup chunks. The models in the vignettes now provide executable demonstrations for users and automated tests.
- Wrapped informational `cat()` messages in `smfa()` with `if (printInfo)` checks. This allows users to suppress progress messages, ensuring the package follows R best practices for silent execution by default.
- Fixed `ic.smfa` to return a data frame (instead of only printing to console), resolving a row-mismatch error in the `efficiency-extraction` vignette.
- Removed executable installation commands (`devtools::install_github`) from `getting-started.Rmd` to comply with CRAN policies.
- Synchronized `smfa.Rd` and `smfa.R` and fixed a minor syntax error in the documentation.

---

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  - "New submission" is expected for a first-time CRAN submission.
  - "Possibly misspelled words in DESCRIPTION" — all domain-specific terms and author names are documented in `inst/WORDLIST`.
  - URL 404s for `smfa` are expected as the package is being renamed and the repository transition is in progress.
  - One note regarding "future file timestamps" was observed locally but is expected to disappear on CRAN servers.

## win-builder (R-devel, x86_64-w64-mingw32)

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  - "New submission" — first-time submission, expected.
  - "Possibly misspelled words in DESCRIPTION" — the flagged words (`Dakpo`, `Metafrontier`, `al`, `et`, `metafrontier`) are correct authors' names and domain-specific terminology.
  - All other sub-items in the NOTE have been addressed.

## Platform

* macOS Ventura 13.7.8, R 4.5.2 (x86_64-apple-darwin20)
* Windows Server 2022, R-devel (win-builder)

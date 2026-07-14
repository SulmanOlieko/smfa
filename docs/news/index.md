# Changelog

## smfa v1.0.1

- Updated the reference for the generalised exponential distribution to
  Papadopoulos (2021) in the documentation.
- Fixed duplicated citation of Papadopoulos (2015) in documentation.
- Added comprehensive, CRAN-compliant unit tests via `testthat`.
- Fixed missing alt text for badges in README files.
- Refactored `README.Rmd` code examples to execute dynamically via
  `knitr` rather than using static mock outputs.

## smfa v1.0.0: Imitation is the sincerest form of flattery

CRAN release: 2026-04-28

The definitive, original implementation of stochastic metafrontier
analysis for R is officially stable.

After rigorous development, extensive methodological testing, and
refinement for CRAN submission, `smfa` v1.0.0 provides a robust,
production-ready environment for productivity and performance
benchmarking across firms operating under different technologies. This
release establishes the standard for stochastic metafrontier analysis in
the R ecosystem.

### Initial CRAN release

- First public release of `smfa`.
- Implements stochastic metafrontier analysis for productivity and
  performance benchmarking across firms operating under different
  technologies.
- Supports three group-frontier types via ‘sfaR’:
  - Standard SFA (`sfacross`)
  - Latent class SFA (`sfalcmcross`)
  - Sample-selection-corrected SFA (`sfaselectioncross`)
- Three metafrontier estimation methods:
  - Linear programming (LP) deterministic envelope
  - Quadratic programming (QP) deterministic envelope
  - Second-stage stochastic frontier (Huang et al. 2014; O’Donnell et
    al. 2008)
- Full efficiency outputs via
  [`efficiencies()`](https://SulmanOlieko.github.io/smfa/reference/efficiencies.md):
  group TE (JLMS and BC), metafrontier TE, and metatechnology ratios
  (MTR).
- Five vignettes illustrating each major use case.

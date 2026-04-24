# Extract information criteria of stochastic metafrontier models

`ic` returns information criterion from stochastic metafrontier models
estimated with
[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

## Usage

``` r
# S3 method for class 'smfa'
ic(object, IC = NULL, ...)
```

## Arguments

- object:

  A stochastic metafrontier model returned by
  [`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

- IC:

  Character string. Information criterion measure. Three criteria are
  available:

  - `'AIC'` for Akaike information criterion (default)

  - `'BIC'` for Bayesian information criterion

  - `'HQIC'` for Hannan-Quinn information criterion

  .

- ...:

  Currently ignored.

## Value

`ic` returns a data frame with the values of the information criteria
(AIC, BIC and HQIC) of the maximum likelihood coefficients. If the `IC`
argument is provided, it returns only the requested criterion as a
numeric value.

## Details

The different information criteria are computed as follows:

- AIC: \\-2 \log{LL} + 2 \* K\\

- BIC: \\-2 \log{LL} + \log{N} \* K\\

- HQIC: \\-2 \log{LL} + 2 \log{\left\[\log{N}\right\]} \* K\\

where \\LL\\ is the maximum likelihood value, \\K\\ the number of
parameters estimated and \\N\\ the number of observations.

## See also

[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), for the
stochastic metafrontier analysis model fitting function using
cross-sectional or pooled data.

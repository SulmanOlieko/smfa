# Extract log-likelihood value of stochastic metafrontier models

`logLik` extracts the log-likelihood value(s) from stochastic
metafrontier models estimated with
[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

## Usage

``` r
# S3 method for class 'smfa'
logLik(object, individual = FALSE, ...)
```

## Arguments

- object:

  A stochastic metafrontier model returned by
  [`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

- individual:

  Logical. If `FALSE` (default), the sum of all observations'
  log-likelihood values is returned. If `TRUE`, a vector of each
  observation's log-likelihood value is returned.

- ...:

  Currently ignored.

## Value

`logLik` returns either an object of class `'logLik'`, which is the
log-likelihood value with the total number of observations (`nobs`) and
the number of free parameters (`df`) as attributes, when
`individual = FALSE`, or a list of elements, containing the
log-likelihood of each observation (`logLik`), the total number of
observations (`Nobs`) and the number of free parameters (`df`), when
`individual = TRUE`.

## See also

[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), for the
stochastic metafrontier analysis model fitting function using
cross-sectional or pooled data.

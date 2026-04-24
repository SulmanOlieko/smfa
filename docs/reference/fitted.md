# Extract fitted values of stochastic metafrontier models

`fitted` returns the fitted frontier values from stochastic metafrontier
models estimated with
[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

## Usage

``` r
# S3 method for class 'smfa'
fitted(object, ...)
```

## Arguments

- object:

  A stochastic metafrontier model returned by
  [`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md).

- ...:

  Currently ignored.

## Value

A vector of fitted values is returned.

## Note

The fitted values are ordered in the same way as the corresponding
observations in the dataset used for the estimation.

## See also

[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), for the
stochastic metafrontier analysis model fitting function using
cross-sectional or pooled data.

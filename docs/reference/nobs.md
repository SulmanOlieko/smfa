# Extract total number of observations used in frontier models

This function extracts the total number of 'observations' from a fitted
point frontier model.

## Usage

``` r
# S3 method for class 'smfa'
nobs(object, ...)
```

## Arguments

- object:

  a `smfa` object for which the number of total observations is to be
  extracted.

- ...:

  Currently ignored.

## Value

A single number, normally an integer.

## Details

`nobs` gives the number of observations actually used by the estimation
procedure.

## See also

[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), for the
stochastic metafrontier analysis model fitting function using
cross-sectional or pooled data

# Extract coefficients of stochastic metafrontier models

From an object of class `'summary.smfa'`, `coef` extracts the
coefficients, their standard errors, z-values, and (asymptotic)
P-values.

From on object of class `'smfa'`, it extracts only the estimated
coefficients.

## Usage

``` r
# S3 method for class 'smfa'
coef(object, ...)

# S3 method for class 'summary.smfa'
coef(object, ...)
```

## Arguments

- object:

  A stochastic metafrontier model returned by
  [`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), or an
  object of class `'summary.smfa'`.

- ...:

  Currently ignored.

## Value

For objects of class `'summary.smfa'`, `coef` returns a matrix with four
columns. Namely, the estimated coefficients, their standard errors,
z-values, and (asymptotic) P-values.

For objects of class `'smfa'`, `coef` returns a numeric vector of the
estimated coefficients.

## See also

[`smfa`](https://SulmanOlieko.github.io/smfa/reference/smfa.md), for the
stochastic metafrontier analysis model fitting function using
cross-sectional or pooled data.

# Direct standardization over an arbitrary standard population

For each provider effect \\\gamma_i\\, the expected number of events
\\E_i = \sum_k N_k\\ \mathrm{expit}(\gamma_i + \eta_k)\\ over the rows
of a standard population with linear predictors \\\eta_k\\ and trials
\\N_k\\, and \\V_i = \sum_k N_k\\ p_k (1 - p_k)\\. The direct
standardized rate is \\E_i / \sum_k N_k\\; its delta-method standard
error with \\\hat\beta\\ held fixed is \\V_i / \sum_k N_k \times
se(\hat\gamma_i)\\.

## Usage

``` r
direct_expectation(gamma, linear_predictor, n = 1, threads = 1)
```

## Arguments

- gamma:

  provider effects.

- linear_predictor:

  \\Z\hat\beta\\ of the standard population's rows.

- n:

  number of trials per row of the standard population (default 1).

- threads:

  number of OpenMP threads.

## Value

A data frame with `expected` and `variance`, one row per element of
`gamma`.

## Examples

``` r
direct_expectation(c(-1, 0, 1), linear_predictor = rnorm(100))
#>   expected variance
#> 1 28.08335 17.26296
#> 2 47.48502 20.64803
#> 3 67.36212 18.19297
```

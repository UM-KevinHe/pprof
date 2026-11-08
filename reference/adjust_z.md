# Adjust statistics with a null model

Adjust statistics with a null model

## Usage

``` r
adjust_z(z, null)
```

## Arguments

- z:

  statistics (a numeric vector, optionally named by provider).

- null:

  a `pp_null` object.

## Value

\\(z - \mathrm{null\\mean}) / \mathrm{null\\sd}\\.

## Examples

``` r
adjust_z(c(a = 2, b = -1), null_model("fixed", sd = 1.5))
#>          a          b 
#>  1.3333333 -0.6666667 
```

# Prepare standardized table for heatmap

Receives a data.frame with Score and creates/ensures the columns: miRNA,
target, gene, utr, Score. Returns a data.frame ready for heatmap
generation.

## Usage

``` r
prepare_for_heatmap(df)
```

## Arguments

- df:

  Data.frame with at least miRNA, target/gene, and Score.

## Value

Data.frame with ordered columns: miRNA, target, gene, utr, Score, and
optional extra columns.

## Examples

``` r
df <- data.frame(miRNA = c("a","b"), target = c("G1_UTR_1","G2"), Score = c(-10, -20))
prepare_for_heatmap(df)
#>   miRNA   target gene  utr Score
#> 1     a G1_UTR_1   G1    1   -10
#> 2     b       G2   G2 <NA>   -20
```

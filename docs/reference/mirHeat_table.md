# Generate a standardized interaction table for mirHeat

This function isolates the tabular output of the mirHeat pipeline. It
returns a clean, standardized miRNA–target interaction table suitable
for inspection or export (e.g. CSV), without generating a heatmap.

## Usage

``` r
mirHeat_table(df)
```

## Arguments

- df:

  A data.frame produced by parse_file(), select_score(), and
  apply_numeric_filters().

## Value

A tibble with standardized columns: miRNA, gene, utr, target, Score.

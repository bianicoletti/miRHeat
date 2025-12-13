# Select which numeric column will be used as Score

Automatically detects numeric columns and defines the selected column as
a new \`Score\` column in the returned data.frame. Interactivity is
avoided: if multiple numeric columns are detected, \`score_column\` must
be explicitly provided.

## Usage

``` r
select_score(df, score_column = NULL)
```

## Arguments

- df:

  Data.frame (usually the output of parse_file()).

- score_column:

  Character or NULL. Name of the numeric column to be used as Score. If
  NULL and only one numeric column exists, it will be used
  automatically.

## Value

A data.frame with a new \`Score\` column

## Examples

``` r
df <- data.frame(miRNA = c("a","b"), target = c("T1","T2"), energy = c(-10, -20))
select_score(df) # uses energy automatically
#> Automatically using the only numeric column detected: energy
#>   miRNA target energy Score
#> 1     a     T1    -10   -10
#> 2     b     T2    -20   -20
```

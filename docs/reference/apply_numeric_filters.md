# Apply numeric filters to Score

Filters the data.frame to keep only rows with Score values within the
specified limits.

## Usage

``` r
apply_numeric_filters(df, min_value = NULL, max_value = NULL, remove_na = TRUE)
```

## Arguments

- df:

  Data.frame containing a \`Score\` column.

- min_value:

  Numeric or NULL. Keeps Score \>= min_value if provided.

- max_value:

  Numeric or NULL. Keeps Score \<= max_value if provided.

- remove_na:

  Logical. If TRUE, removes rows with NA Score values.

## Value

A filtered data.frame.

## Examples

``` r
df <- data.frame(miRNA = c("a","b","c"), Score = c(-10, NA, -30))
apply_numeric_filters(df, min_value = -25)
#>   miRNA Score
#> 1     a   -10
```

# Aplica filtros numericos ao Score

Filtra o data.frame para manter apenas linhas com Score dentro dos
limites.

## Usage

``` r
apply_numeric_filters(df, min_value = NULL, max_value = NULL, remove_na = TRUE)
```

## Arguments

- df:

  Data.frame com coluna \`Score\`.

- min_value:

  Numeric or NULL. Mantem Score \>= min_value se fornecido.

- max_value:

  Numeric or NULL. Mantem Score \<= max_value se fornecido.

- remove_na:

  Logical. Se TRUE remove linhas com Score NA.

## Value

Data.frame filtrado.

## Examples

``` r
df <- data.frame(miRNA = c("a","b","c"), Score = c(-10, NA, -30))
apply_numeric_filters(df, min_value = -25)
#>   miRNA Score
#> 1     a   -10
```

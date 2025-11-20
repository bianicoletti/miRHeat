# Seleciona qual coluna numerica sera usada como Score

Detecta automaticamente colunas numericas e define a coluna selecionada
como nova coluna \`Score\` no data.frame retornado. Evita
interatividade: se houver multiplas colunas numericas requer que
\`score_column\` seja fornecido.

## Usage

``` r
select_score(df, score_column = NULL)
```

## Arguments

- df:

  Data.frame (geralmente resultado de parse_file()).

- score_column:

  Character or NULL. Nome da coluna numerica a usar como Score. Se NULL
  e apenas uma coluna numerica existir, ela sera usada automaticamente.

## Value

Data.frame com nova coluna \`Score\`.

## Examples

``` r
df <- data.frame(miRNA = c("a","b"), target = c("T1","T2"), energy = c(-10, -20))
select_score(df) # usa energy automaticamente
#> Usando automaticamente a unica coluna numerica encontrada: energy
#>   miRNA target energy Score
#> 1     a     T1    -10   -10
#> 2     b     T2    -20   -20
```

# Prepara tabela padronizada para heatmap

Recebe o data.frame com Score e cria/garante as colunas: miRNA, target,
gene, utr, Score. Retorna um data.frame pronto para o heatmap.

## Usage

``` r
prepare_for_heatmap(df)
```

## Arguments

- df:

  Data.frame com pelo menos miRNA, target/gene e Score.

## Value

Data.frame com colunas ordenadas: miRNA, target, gene, utr, Score e
outras colunas opcionais.

## Examples

``` r
df <- data.frame(miRNA = c("a","b"), target = c("G1_UTR_1","G2"), Score = c(-10, -20))
prepare_for_heatmap(df)
#>   miRNA   target gene  utr Score
#> 1     a G1_UTR_1   G1    1   -10
#> 2     b       G2   G2 <NA>   -20
```

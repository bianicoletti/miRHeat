# Gera Heatmap de Interacoes miRNA x UTR

Gera um heatmap dos scores de interacao entre miRNA e UTR. Determina
automaticamente o numero ideal de clusters usando corte dinamico de
arvores (dynamic tree cutting).

## Usage

``` r
plot_mirheat(df, output_file = NULL, width = 2000, height = 1800)
```

## Arguments

- df:

  Data frame contendo as colunas miRNA, UTR (ou utr) e Score.

- output_file:

  Caminho opcional para exportar o PNG (padrao = NULL).

- width:

  Largura do PNG ao exportar (em pixels).

- height:

  Altura do PNG ao exportar (em pixels).

## Value

Retorna invisivelmente o objeto ComplexHeatmap.

## Examples

``` r
if (FALSE) { # \dontrun{
file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")
df <- parse_file(file)
df <- select_score(df, "interaction_energy")
df <- prepare_for_heatmap(df)
plot_mirheat(df)
} # }
```

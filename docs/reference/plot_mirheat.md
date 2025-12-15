# Generate miRNA vs UTR interaction heatmap

Generates a heatmap of interaction scores between miRNAs and UTRs.
Automatically determines the optimal number of clusters using dynamic
tree cutting.

## Usage

``` r
plot_mirheat(df, output_file = NULL, width = 2000, height = 1800)
```

## Arguments

- df:

  Data.frame containing the columns miRNA, UTR (or utr), and Score.

- output_file:

  Optional path to export the PNG file (default = NULL).

- width:

  Width of the PNG when exporting (in pixels).

- height:

  Height of the PNG when exporting (in pixels).

## Value

Invisibly returns the ComplexHeatmap object.

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

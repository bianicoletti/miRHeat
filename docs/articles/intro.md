# Introduction to miRHeat: from input files to heatmap

## Introduction

The **miRHeat** package provides a complete pipeline to:

- Read output files from **IntaRNA** and **RNAHybrid** for exemple
- Extract and select numeric columns as **Score**
- Standardize the data structure for visualization purposes
- Generate a **miRNA × UTR interaction heatmap** with automatic
  clustering

It was designed to facilitate the visual inspection of predicted
interactions between miRNAs and gene 3’UTR regions.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("bianicoletti/miRHeat")
```

Load the package::

``` r
library(miRHeat)
```

## Example file

The package includes an example file that can be used for demonstration
purposes:

``` r
file <- system.file("extdata", "example_blocks.txt", package = "miRHeat")
file
```

## 1. Reading and structuring the data

The first step is to use
[`parse_file()`](https://github.com/bianicoletti/miRHeat/reference/parse_file.md):

``` r
df <- parse_file(arquivo)
head(df)
```

## 2. Choosing the Score metric

Predictive softwares can produce multiple numeric columns (energy, ΔG,
etc.).

The
[`select_score()`](https://github.com/bianicoletti/miRHeat/reference/select_score.md)
function automatically identifies the only numeric column or allows you
to specify which one should be used:

``` r
df <- select_score(df, score_column = "interaction_energy")
head(df$Score)
```

## 3. Preparing the data for the heatmap

The
[`prepare_for_heatmap()`](https://github.com/bianicoletti/miRHeat/reference/prepare_for_heatmap.md)
function ensures that all required columns are present:

- `miRNA`
- `target`
- `gene`
- `utr`
- `Score`

``` r
df <- prepare_for_heatmap(df)
head(df)
```

## 4. Optional filtering

Numeric filters can be applied naturally:

``` r
df_filt <- apply_numeric_filters(df, Score < -10)
```

## 5. Generating the final heatmap

The main visualization function is:

``` r
plot_mirheat(df)
```

The plot:

- automatically builds the miRNA × UTR matrix
- handles missing values
- computes distances and clusters
- determines the optimal number of groups
- uses`ComplexHeatmap` for rendering

## 6. Exporting the heatmap to a file

``` r
plot_mirheat(df, output_file = "heatmap_miRHeat.png")
```

## PComplete Pipeline

``` r
file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")

df <- parse_file(file)
df <- select_score(df, "interaction_energy")
df <- prepare_for_heatmap(df)

plot_mirheat(df)
```

## Conclusions

\*miRHeat\*\* provides a simple and automated workflow to analyze
miRNA–3’UTR interactions predicted by IntaRNA — from raw block parsing
to the final heatmap visualization.

``` r
library(miRHeat)
```

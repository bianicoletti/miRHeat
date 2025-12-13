# Parse interaction files (block-formatted TXT or CSV/TSV)

Reads an interaction file either in block format (as produced by
RNAhybrid/IntaRNA) or as a tabular CSV/TSV file. For block-formatted TXT
files, the function assumes the following structure: a line with an
identifier such as \`GENE_UTR_123\`, followed by alignment lines, a line
containing the miRNA name (e.g. \`hsa-miR-7160-5p\`), and at least one
line with a \`key = value\` pair (e.g. \`interaction energy = -24.54
kcal/mol\`).

## Usage

``` r
parse_file(
  file,
  pattern = "^.*_UTR_.*",
  guess_tabular = TRUE,
  mirna_col = NULL,
  target_col = NULL
)
```

## Arguments

- file:

  Path to the input file.

- pattern:

  Character string. Regular expression pattern used to detect the
  beginning of a block (default: '^.\*\_UTR\_.\*').

- guess_tabular:

  Logical. If TRUE, the function attempts to automatically detect
  whether the file is tabular.

- mirna_col:

  Optional character. Name of the column containing miRNAs (applies to
  tabular files).

- target_col:

  Optional character. Name of the column containing targets/genes
  (applies to tabular files).

## Value

A data.frame containing the extracted columns (miRNA, target, gene, utr,
and any detected numeric values).

## Details

For tabular files, the function attempts to automatically detect miRNA
and target columns.

## Examples

``` r
if (FALSE) { # \dontrun{
file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")
df <- parse_file(file)
head(df)
} # }
```

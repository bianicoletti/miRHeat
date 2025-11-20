# Parse de arquivo de interacoes (TXT em blocos ou CSV/TSV)

Le um arquivo de interacoes (formato em blocos tipo RNAhybrid/IntaRNA ou
arquivo tabular CSV/TSV). Para arquivos TXT em blocos assume o padrao
informado: uma linha com identificador do tipo \`GENE_UTR_123\`, seguida
por linhas de alinhamento, uma linha com o miRNA (ex:
\`hsa-miR-7160-5p\`) e pelo menos uma linha com \`chave = valor\` (ex:
\`interaction energy = -24.54 kcal/mol\`).

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

  Path para o arquivo de entrada.

- pattern:

  String. Padrao regex para detectar inicio de bloco (default
  '^.\*\_UTR\_.\*').

- guess_tabular:

  Logical. Se TRUE, tenta detectar automaticamente se o arquivo e
  tabular.

- mirna_col:

  Optional character. Nome da coluna que contem miRNAs (aplica-se a
  arquivos tabulares).

- target_col:

  Optional character. Nome da coluna que contem targets/genes (aplica-se
  a arquivos tabulares).

## Value

Um data.frame com colunas extraidas (miRNA, target, gene, utr e
quaisquer valores numericos detectados).

## Details

Para arquivos tabulares, tenta detectar colunas de miRNA e target
automaticamente.

## Examples

``` r
if (FALSE) { # \dontrun{
file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")
df <- parse_file(file)
head(df)
} # }
```

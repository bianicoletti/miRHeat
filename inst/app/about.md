
# Sobre o App miRHeat

Este aplicativo é uma interface interativa do pacote **miRHeat**, que
permite processar arquivos de interações miRNA–UTR usando a função
`parse_file()`.

## Funcionalidades

- Carregar arquivos tabulares ou texto (.txt, .csv, .tsv, etc.)
- Aplicar **expressões regulares** (pattern) para filtrar linhas/colunas
- Detectar automaticamente o formato tabular (opção “Detectar formato
  automaticamente”)
- Visualizar os dados processados em uma tabela interativa
- Mensagens informativas sobre o processamento (sucesso, erro ou sem
  resultados)

## Como usar

1.  Faça upload do seu arquivo na barra lateral.
2.  Ajuste o **Pattern** para filtrar os dados.
3.  Marque ou desmarque a opção de detecção automática do formato.
4.  Veja a prévia da tabela processada no painel principal.

## Sobre o pacote

O app utiliza o pacote **miRHeat**, que contém funções para:

- Parsing de arquivos (`parse_file`)
- Seleção de scores (`select_score`)
- Aplicação de filtros (`apply_numeric_filters`)
- Preparação de dados para heatmaps (`prepare_for_heatmap`)
- Plotagem de heatmaps com clusters (`plot_mirheat`)

## Contato

Desenvolvido por **Maria Nicoletti**.

#' Parse de arquivo de interações (TXT em blocos ou CSV/TSV)
#'
#' Lê um arquivo de interações (formato em blocos tipo RNAhybrid/IntaRNA ou
#' arquivo tabular CSV/TSV). Para arquivos TXT em blocos assume o padrão informado:
#' uma linha com identificador do tipo `GENE_UTR_123`, seguida por linhas de alinhamento,
#' uma linha com o miRNA (ex: `hsa-miR-7160-5p`) e pelo menos uma linha com `chave = valor`
#' (ex: `interaction energy = -24.54 kcal/mol`).
#'
#' Para arquivos tabulares, tenta detectar colunas de miRNA e target automaticamente.
#'
#' @param file Path para o arquivo de entrada.
#' @param pattern String. Padrão regex para detectar início de bloco (default '^.*_UTR_.*').
#' @param guess_tabular Logical. Se TRUE, tenta detectar automaticamente se o arquivo é tabular.
#' @return Um data.frame com colunas extraídas (miRNA, target_raw, e quaisquer valores numéricos detectados).
#' @examples
#' \dontrun{
#' df <- parse_file("exemplo_blocos.txt")
#' }
#' @export
parse_file <- function(file, pattern = "^.*_UTR_.*", guess_tabular = TRUE) {
  if (!file.exists(file)) stop("Arquivo não encontrado: ", file)

  ext <- tolower(tools::file_ext(file))
  # Helper: normaliza nomes de colunas
  clean_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    tolower(trimws(x))
  }

  # ----- Caso tabular (csv/tsv) -----
  if (guess_tabular && ext %in% c("csv", "tsv", "txt")) {
    # tenta inferir se há separador na primeira linha
    first_line <- readLines(file, n = 1, warn = FALSE)
    if (grepl(",", first_line) && ext == "csv") {
      df_tab <- tryCatch(utils::read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else if (grepl("\t", first_line) && ext %in% c("tsv", "txt")) {
      df_tab <- tryCatch(utils::read.delim(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else {
      df_tab <- NULL
    }

    if (!is.null(df_tab)) {
      # detecta colunas candidatas
      names(df_tab) <- clean_name(names(df_tab))
      nms <- names(df_tab)
      # tenta encontrar miRNA column
      mirna_col <- nms[grepl("mirna|mi-rna|mir-|mir", nms)]
      target_col <- nms[grepl("gene|target|utr", nms)]
      # detecta colunas numéricas que pareçam scores
      score_col <- nms[vapply(df_tab, is.numeric, logical(1))]

      # se não achar explicitamente miRNA/target, pergunta interativamente
      if (length(mirna_col) == 0 || length(target_col) == 0) {
        message("Colunas detectadas no arquivo tabular:")
        print(nms)
        if (length(mirna_col) == 0) {
          mirna_col <- readline(prompt = "Digite o nome da coluna que contém os miRNAs (ou deixe em branco para sair): ")
          if (mirna_col == "") stop("miRNA column não fornecida.")
          mirna_col <- clean_name(mirna_col)
        }
        if (length(target_col) == 0) {
          target_col <- readline(prompt = "Digite o nome da coluna que contém os targets (gene/UTR) (ou deixe em branco para sair): ")
          if (target_col == "") stop("target column não fornecida.")
          target_col <- clean_name(target_col)
        }
      } else {
        mirna_col <- mirna_col[1]
        target_col <- target_col[1]
      }

      df_out <- data.frame(
        miRNA = df_tab[[mirna_col]],
        target_raw = df_tab[[target_col]],
        stringsAsFactors = FALSE
      )

      # anexa colunas numéricas detectadas (scores)
      for (sc in score_col) {
        # evita sobrescrever colunas base
        if (!(sc %in% c(mirna_col, target_col))) {
          df_out[[clean_name(sc)]] <- df_tab[[sc]]
        }
      }
      return(df_out)
    }
  }

  # ----- Caso TXT em blocos -----
  linhas <- readLines(file, warn = FALSE)
  indices_inicio <- grep(pattern, linhas)
  if (length(indices_inicio) == 0) {
    stop("Nenhum bloco encontrado com o padrão: ", pattern)
  }
  indices_fim <- c((indices_inicio[-1] - 1), length(linhas))
  blocos <- Map(function(i, f) linhas[i:f], indices_inicio, indices_fim)

  # funcoes auxiliares de extração
  extrai_miRNA <- function(bloco) {
    # procura padrão tipo hsa-miR-... ou miR...
    pat <- "(hsa-)?miR[-A-Za-z0-9_]+"
    hit <- unique(unlist(regmatches(bloco, gregexpr(pat, bloco, ignore.case = TRUE))))
    if (length(hit) >= 1) {
      return(hit[1])
    }
    # fallback: procura linha que contenha 'miR' ignorando case
    idx <- grep("mir", bloco, ignore.case = TRUE)
    if (length(idx) >= 1) {
      cand <- trimws(bloco[idx[1]])
      return(cand)
    }
    return(NA_character_)
  }
  extrai_target <- function(first_line) {
    txt <- trimws(first_line)
    # tenta extrair gene e utr: GENE_UTR_123 ou variantes
    m <- regexec("^([A-Za-z0-9]+)[_\\.-]*UTR[_\\.-]*([0-9]+)?", txt, ignore.case = TRUE)
    mm <- regmatches(txt, m)[[1]]
    if (length(mm) >= 2) {
      gene <- mm[2]
      utrid <- ifelse(length(mm) >= 3 && nzchar(mm[3]), mm[3], NA_character_)
      return(list(gene = toupper(gene), utr = utrid, target = ifelse(is.na(utrid), gene, paste0(gene, "_UTR_", utrid))))
    } else {
      # se não bater, usa toda a primeira linha como target bruto
      return(list(gene = txt, utr = NA_character_, target = txt))
    }
  }
  extrai_pairs_numeric <- function(bloco) {
    # captura pares "nome = valor" onde valor é numérico (com sinal opcional e decimal)
    pat <- "([A-Za-z0-9\\+\\-\\_\\s\\+\\+\\.\\%]+?)\\s*=\\s*(-?[0-9]+\\.?[0-9]*)"
    hits <- regmatches(bloco, gregexpr(pat, bloco, perl = TRUE, ignore.case = TRUE))
    pairs <- unlist(hits)
    if (length(pairs) == 0) return(list())
    # processa cada match: queremos nome (limpo) e valor numérico
    res <- list()
    for (p in pairs) {
      m <- regexec(pat, p, perl = TRUE, ignore.case = TRUE)
      mm <- regmatches(p, m)[[1]]
      if (length(mm) >= 3) {
        name_raw <- trimws(mm[2])
        # limpa o nome pra usar como coluna
        name_col <- gsub("[^A-Za-z0-9_]", "_", tolower(name_raw))
        name_col <- gsub("_+", "_", name_col)
        val <- as.numeric(mm[3])
        res[[name_col]] <- val
      }
    }
    return(res)
  }

  # itera blocos
  rows <- vector("list", length(blocos))
  for (i in seq_along(blocos)) {
    b <- blocos[[i]]
    mi <- extrai_miRNA(b)
    target_info <- extrai_target(b[1])
    pairs <- extrai_pairs_numeric(b)
    rows[[i]] <- c(list(miRNA = mi, target_raw = target_info$target, gene = target_info$gene, utr = target_info$utr), pairs)
  }

  # une em data.frame (preenchendo NA onde necessário)
  df <- do.call(rbind, lapply(rows, function(x) {
    # garante que todos os nomes se alinhem
    as.data.frame(lapply(x, function(y) if (is.null(y)) NA else y), stringsAsFactors = FALSE)
  }))
  # corrige nomes e tipos numéricos
  # converter colunas que parecem numéricas
  for (nm in names(df)) {
    if (all(is.na(df[[nm]]) | grepl("^-?[0-9]+\\.?[0-9]*$", as.character(df[[nm]])))) {
      df[[nm]] <- as.numeric(df[[nm]])
    }
  }
  # renomeia target_raw para target se existir
  if ("target_raw" %in% names(df)) names(df)[names(df) == "target_raw"] <- "target"
  # organiza colunas: miRNA, target, gene, utr, outros...
  core <- c("miRNA", "target", "gene", "utr")
  others <- setdiff(names(df), core)
  df <- df[c(core, others)]
  rownames(df) <- NULL
  return(df)
}


#' Seleciona qual coluna numérica será usada como Score
select_score <- function(df, score_column = NULL) {
  if (!"data.frame" %in% class(df)) {
    stop("df deve ser um data.frame produzido por parse_file().")
  }

  # Descobre automaticamente colunas numéricas
  numeric_cols <- names(df)[sapply(df, is.numeric)]

  if (length(numeric_cols) == 0) {
    stop("Nenhuma coluna numérica encontrada para ser usada como Score.")
  }

  # Se o usuário não informou a coluna, mas só existe uma → usar ela
  if (is.null(score_column)) {
    if (length(numeric_cols) == 1) {
      score_column <- numeric_cols[1]
      message("Usando automaticamente a única coluna numérica encontrada: ", score_column)
    } else {
      stop(
        paste0(
          "Existem múltiplas colunas numéricas: ",
          paste(numeric_cols, collapse = ", "),
          "\nEspecifique qual deseja usar com score_column = \"nome_da_coluna\"."
        )
      )
    }
  }

  # Verifica se o nome existe
  if (!score_column %in% names(df)) {
    stop("Coluna selecionada não encontrada no data.frame: ", score_column)
  }

  df$Score <- df[[score_column]]
  df
}



#' Aplica filtros numéricos ao Score
#'
#' Filtra o data.frame para manter apenas linhas com Score dentro dos limites.
#'
#' @param df Data.frame com coluna `Score`.
#' @param min_value Numeric or NULL. Mantém Score >= min_value se fornecido.
#' @param max_value Numeric or NULL. Mantém Score <= max_value se fornecido.
#' @param remove_na Logical. Se TRUE remove linhas com Score NA.
#' @return Data.frame filtrado.
#' @examples
#' \dontrun{
#' df_f <- apply_numeric_filters(df2, min_value = -20)
#' }
#' @export
apply_numeric_filters <- function(df, min_value = NULL, max_value = NULL, remove_na = TRUE) {
  if (!"Score" %in% names(df)) stop("Data.frame deve conter a coluna 'Score'. Use select_score() antes.")
  keep <- rep(TRUE, nrow(df))
  if (!is.null(min_value)) keep <- keep & (!is.na(df$Score) & df$Score >= min_value)
  if (!is.null(max_value)) keep <- keep & (!is.na(df$Score) & df$Score <= max_value)
  if (remove_na) keep <- keep & !is.na(df$Score)
  df2 <- df[keep, , drop = FALSE]
  rownames(df2) <- NULL
  return(df2)
}


#' Prepara tabela padronizada para heatmap
#'
#' Recebe o data.frame com Score e cria/garante as colunas:
#' miRNA, Target, Gene, UTR, Score. Retorna um data.frame pronto para o heatmap.
#'
#' @param df Data.frame com pelo menos miRNA, target/gene e Score.
#' @return Data.frame com colunas ordenadas: miRNA, Target, Gene, UTR, Score, e outras colunas opcionais.
#' @examples
#' \dontrun{
#' df_ready <- prepare_for_heatmap(df_filtered)
#' }
#' @export
prepare_for_heatmap <- function(df) {
  if (!is.data.frame(df)) stop("df deve ser um data.frame.")
  # checagens
  if (!("miRNA" %in% names(df))) stop("Coluna 'miRNA' não encontrada.")
  if (!("Score" %in% names(df))) stop("Coluna 'Score' não encontrada. Use select_score().")
  # target: se não existir, tenta construir a partir de gene/utr
  if (!("target" %in% names(df))) {
    if ("gene" %in% names(df) && "utr" %in% names(df)) {
      df$target <- ifelse(is.na(df$utr), df$gene, paste0(df$gene, "_UTR_", df$utr))
    } else if ("gene" %in% names(df)) {
      df$target <- df$gene
    } else {
      stop("Impossível construir 'target' a partir dos dados existentes.")
    }
  }
  # Garante colunas gene/utr
  if (!("gene" %in% names(df))) {
    # tenta extrair gene de target (antes do _UTR_)
    df$gene <- sub("(_UTR_.*$)|([_\\.-]UTR[_\\.-].*$)", "", df$target)
  }
  if (!("utr" %in% names(df))) {
    # tenta extrair UTR id se houver
    utrs <- regmatches(df$target, regexec("_UTR_([0-9]+)$", df$target))
    df$utr <- sapply(utrs, function(x) if (length(x) >= 2) x[2] else NA_character_)
  }
  # organiza colunas
  core <- c("miRNA", "target", "gene", "utr", "Score")
  rest <- setdiff(names(df), core)
  df_out <- df[c(core, rest)]
  rownames(df_out) <- NULL
  return(df_out)
}

#' Plot Heatmap for miRNA–UTR Interactions
#'
#' @description
#' Generates a heatmap of miRNA × UTR interaction scores.
#' Automatically determines the optimal number of clusters using
#' dynamic tree cutting.
#'
#' @param df Data frame containing miRNA, UTR, and Score columns.
#' @param output_file Optional path to export PNG (default = NULL).
#' @param width Width of PNG if exporting (in pixels).
#' @param height Height of PNG if exporting (in pixels).
#'
#' @return Returns the ComplexHeatmap object invisibly.
#'
#' @export
plot_mirheat <- function(df,
                         output_file = NULL,
                         width = 2000,
                         height = 1800) {

  # Required packages
  require(reshape2)
  require(ComplexHeatmap)
  require(circlize)
  require(dynamicTreeCut)

  message("→ Construindo matriz miRNA × UTR...")

  # Reformat
  matriz <- reshape2::dcast(df, miRNA ~ utr, value.var = "Score")
  rownames(matriz) <- matriz$miRNA
  matriz$miRNA <- NULL

  # Replace NAs with mean of row
  matriz <- t(apply(matriz, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }))

  # Distance + clustering
  d_miRNA <- dist(matriz)
  hc_miRNA <- hclust(d_miRNA, method = "ward.D2")

  message("→ Detectando automaticamente número de clusters...")

  clusters <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_miRNA,
    distM = as.matrix(d_miRNA),
    deepSplit = 2,
    pamRespectsDendro = TRUE
  )

  message(paste("→ Número detectado:", length(unique(clusters)), "clusters."))

  annotation <- HeatmapAnnotation(
    Cluster = clusters,
    col = list(Cluster = structure(
      circlize::rand_color(length(unique(clusters))),
      names = unique(clusters)
    ))
  )

  # Palette for scores
  paleta <- circlize::colorRamp2(
    seq(min(matriz), max(matriz), length.out = 5),
    c("#b2182b", "#f4a582", "#f7f7f7", "#92c5de", "#2166ac")
  )

  # Heatmap
  ht <- Heatmap(
    matriz,
    name = "Score",
    col = paleta,
    cluster_rows = hc_miRNA,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    left_annotation = rowAnnotation(Cluster = clusters)
  )

  # Export PNG if requested
  if (!is.null(output_file)) {
    message("→ Exportando heatmap para: ", output_file)
    png(output_file, width = width, height = height, res = 200)
    draw(ht)
    dev.off()
  }

  message("→ Heatmap pronto!")
  invisible(ht)
}

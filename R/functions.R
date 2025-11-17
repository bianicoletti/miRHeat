#' Parse de arquivo de interacoes (TXT em blocos ou CSV/TSV)
#'
#' Le um arquivo de interacoes (formato em blocos tipo RNAhybrid/IntaRNA ou
#' arquivo tabular CSV/TSV). Para arquivos TXT em blocos assume o padrao informado:
#' uma linha com identificador do tipo `GENE_UTR_123`, seguida por linhas de alinhamento,
#' uma linha com o miRNA (ex: `hsa-miR-7160-5p`) e pelo menos uma linha com `chave = valor`
#' (ex: `interaction energy = -24.54 kcal/mol`).
#'
#' Para arquivos tabulares, tenta detectar colunas de miRNA e target automaticamente.
#'
#' @param file Path para o arquivo de entrada.
#' @param pattern String. Padrao regex para detectar inicio de bloco (default '^.*_UTR_.*').
#' @param guess_tabular Logical. Se TRUE, tenta detectar automaticamente se o arquivo e tabular.
#' @param mirna_col Optional character. Nome da coluna que contem miRNAs (aplica-se a arquivos tabulares).
#' @param target_col Optional character. Nome da coluna que contem targets/genes (aplica-se a arquivos tabulares).
#' @return Um data.frame com colunas extraidas (miRNA, target, gene, utr e quaisquer valores numericos detectados).
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")
#' df <- parse_file(file)
#' head(df)
#' }
#' @export
parse_file <- function(file,
                       pattern = "^.*_UTR_.*",
                       guess_tabular = TRUE,
                       mirna_col = NULL,
                       target_col = NULL) {
  if (!file.exists(file)) stop("Arquivo não encontrado: ", file)

  ext <- tolower(tools::file_ext(file))

  # Helper: normaliza nomes de colunas
  clean_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    tolower(trimws(x))
  }

  # ----- Caso tabular (csv/tsv/txt) -----
  if (guess_tabular && ext %in% c("csv", "tsv", "txt")) {
    # tenta inferir separador pela primeira linha
    first_line <- readLines(file, n = 1, warn = FALSE)
    df_tab <- NULL
    if (grepl(",", first_line) && ext %in% c("csv", "txt")) {
      df_tab <- tryCatch(utils::read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else if (grepl("\t", first_line) && ext %in% c("tsv", "txt")) {
      df_tab <- tryCatch(utils::read.delim(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else if (grepl(";", first_line) && ext %in% c("csv","txt")) {
      # possível CSV com ponto-e-vírgula
      df_tab <- tryCatch(utils::read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
    }

    if (!is.null(df_tab) && is.data.frame(df_tab) && ncol(df_tab) > 0) {
      names(df_tab) <- clean_name(names(df_tab))
      nms <- names(df_tab)

      # Se o usuário forneceu nomes, usa-os (após limpeza)
      if (!is.null(mirna_col)) mirna_col <- clean_name(mirna_col)
      if (!is.null(target_col)) target_col <- clean_name(target_col)

      # tenta encontrar miRNA/target automaticamente se não fornecido
      if (is.null(mirna_col)) {
        mirna_candidates <- nms[grepl("mirna|mirna_|mi-rna|mi-rna|\\bmir\\b|mi_r|mi-r", nms, ignore.case = TRUE)]
        if (length(mirna_candidates) >= 1) mirna_col <- mirna_candidates[1]
      }
      if (is.null(target_col)) {
        target_candidates <- nms[grepl("gene|target|utr", nms, ignore.case = TRUE)]
        if (length(target_candidates) >= 1) target_col <- target_candidates[1]
      }

      # Se ainda faltam colunas, aborta com mensagem informativa (sem interatividade)
      missing <- character(0)
      if (is.null(mirna_col)) missing <- c(missing, "miRNA")
      if (is.null(target_col)) missing <- c(missing, "target/gene/utr")
      if (length(missing) > 0) {
        stop(
          "Não foi possível identificar automaticamente as colunas tabulares necessárias: ",
          paste(missing, collapse = ", "),
          ". Colunas detectadas: ", paste(nms, collapse = ", "),
          ".\nSolução: reexporte o arquivo com cabeçalho claro ou chame parse_file(..., mirna_col = 'nome', target_col = 'nome')."
        )
      }

      # detecta colunas numéricas que pareçam scores (exceto colunas chave)
      score_col <- setdiff(nms[vapply(df_tab, is.numeric, logical(1))], c(mirna_col, target_col))

      df_out <- data.frame(
        miRNA = df_tab[[mirna_col]],
        target = df_tab[[target_col]],
        stringsAsFactors = FALSE
      )

      # anexa colunas numéricas detectadas (scores)
      for (sc in score_col) {
        out_name <- clean_name(sc)
        # evita sobrescrever colunas base
        if (!(out_name %in% c("mirna", "target"))) {
          df_out[[out_name]] <- df_tab[[sc]]
        }
      }

      return(df_out)
    }
    # se df_tab não foi detectado, continua para modo blocos
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
    # procura padrões como hsa-miR-..., miR-..., miRNA-..., miR...
    pat <- "(?:hsa-)?(?:miR|miRNA)[-A-Za-z0-9_]+"
    hit <- unique(unlist(regmatches(bloco, gregexpr(pat, bloco, ignore.case = TRUE))))
    if (length(hit) >= 1 && nzchar(hit[1])) {
      return(hit[1])
    }
    # fallback: procura linha que contenha 'miR' ou 'mirna' ignorando case
    idx <- grep("mirna|mir", bloco, ignore.case = TRUE)
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
    # captura pares "nome = valor" onde valor é numérico (com sinal opcional e decimal),
    # aceita percentuais em nome, ignora unidades
    pat <- "([A-Za-z0-9\\+\\-\\_\\s\\.\\%]+?)\\s*=\\s*(-?[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?)"
    hits <- regmatches(bloco, gregexpr(pat, bloco, perl = TRUE, ignore.case = TRUE))
    pairs <- unlist(hits)
    if (length(pairs) == 0) return(list())
    res <- list()
    for (p in pairs) {
      m <- regexec(pat, p, perl = TRUE, ignore.case = TRUE)
      mm <- regmatches(p, m)[[1]]
      if (length(mm) >= 3) {
        name_raw <- trimws(mm[2])
        name_col <- gsub("[^A-Za-z0-9_]", "_", tolower(name_raw))
        name_col <- gsub("_+", "_", name_col)
        # remove trailing/leading underscores
        name_col <- gsub("^_|_$", "", name_col)
        val <- as.numeric(mm[3])
        # se já existir, torna em vetor (mantém último valor se necessário)
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
    rows[[i]] <- c(list(miRNA = mi, target = target_info$target, gene = target_info$gene, utr = target_info$utr), pairs)
  }

  # une em data.frame (preenchendo NA onde necessário)
  # coletar todos nomes possíveis
  all_names <- unique(unlist(lapply(rows, names)))
  df_list <- lapply(rows, function(x) {
    # garante que todos os nomes existam e sejam NA se ausentes
    res <- setNames(vector("list", length(all_names)), all_names)
    for (nm in all_names) {
      res[[nm]] <- if (!is.null(x[[nm]])) x[[nm]] else NA
    }
    as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  })
  df <- do.call(rbind, df_list)

  # Evitar conversão de colunas chave para numérico
  core_cols <- c("miRNA", "target", "gene", "utr")
  other_cols <- setdiff(names(df), core_cols)

  # Converter colunas "other" que parecem numéricas
  for (nm in other_cols) {
    # se toda a coluna for NA ou números no formato, converte
    vals_chr <- as.character(df[[nm]])
    if (all(is.na(vals_chr) | grepl("^\\s*-?[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?\\s*$", vals_chr))) {
      df[[nm]] <- as.numeric(vals_chr)
    }
  }

  # renomeia target para target (já OK) e organiza colunas
  # padroniza nomes: se existir 'target' ok; se estiver 'target_raw', renomeia
  if ("target_raw" %in% names(df) && !("target" %in% names(df))) {
    names(df)[names(df) == "target_raw"] <- "target"
  }

  # organiza colunas: miRNA, target, gene, utr, outros...
  core <- c("miRNA", "target", "gene", "utr")
  others <- setdiff(names(df), core)
  df <- df[c(core, others)]
  rownames(df) <- NULL
  return(df)
}


#' Seleciona qual coluna numerica sera usada como Score
#'
#' Detecta automaticamente colunas numericas e define a coluna selecionada como
#' nova coluna `Score` no data.frame retornado. Evita interatividade: se houver
#' multiplas colunas numericas requer que `score_column` seja fornecido.
#'
#' @param df Data.frame (geralmente resultado de parse_file()).
#' @param score_column Character or NULL. Nome da coluna numerica a usar como Score.
#'                     Se NULL e apenas uma coluna numerica existir, ela sera usada automaticamente.
#' @return Data.frame com nova coluna `Score`.
#' @examples
#' df <- data.frame(miRNA = c("a","b"), target = c("T1","T2"), energy = c(-10, -20))
#' select_score(df) # usa energy automaticamente
#' @export


select_score <- function(df, score_column = NULL) {
  if (!is.data.frame(df)) {
    stop("df deve ser um data.frame produzido por parse_file().")
  }

  # Descobre automaticamente colunas numéricas, exceto colunas chave
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("utr", "target", "miRNA", "gene"))

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
          ". Especifique qual deseja usar com score_column = \"nome_da_coluna\"."
        )
      )
    }
  }

  # Verifica se o nome existe e é numérico
  if (!score_column %in% names(df)) {
    stop("Coluna selecionada não encontrada no data.frame: '", score_column, "'.")
  }
  if (!is.numeric(df[[score_column]])) {
    stop("A coluna selecionada não é numérica: '", score_column, "'.")
  }

  df$Score <- df[[score_column]]
  return(df)
}


#' Aplica filtros numericos ao Score
#'
#' Filtra o data.frame para manter apenas linhas com Score dentro dos limites.
#'
#' @param df Data.frame com coluna `Score`.
#' @param min_value Numeric or NULL. Mantem Score >= min_value se fornecido.
#' @param max_value Numeric or NULL. Mantem Score <= max_value se fornecido.
#' @param remove_na Logical. Se TRUE remove linhas com Score NA.
#' @return Data.frame filtrado.
#' @examples
#' df <- data.frame(miRNA = c("a","b","c"), Score = c(-10, NA, -30))
#' apply_numeric_filters(df, min_value = -25)
#' @export
apply_numeric_filters <- function(df, min_value = NULL, max_value = NULL, remove_na = TRUE) {
  if (!is.data.frame(df)) stop("df deve ser um data.frame.")
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
#' miRNA, target, gene, utr, Score. Retorna um data.frame pronto para o heatmap.
#'
#' @param df Data.frame com pelo menos miRNA, target/gene e Score.
#' @return Data.frame com colunas ordenadas: miRNA, target, gene, utr, Score e outras colunas opcionais.
#' @examples
#' df <- data.frame(miRNA = c("a","b"), target = c("G1_UTR_1","G2"), Score = c(-10, -20))
#' prepare_for_heatmap(df)
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
    df$gene <- sub("(_UTR_.*$)|([_\\.-]UTR[_\\.-].*$)", "", df$target)
  }
  if (!("utr" %in% names(df))) {
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


#' Plot Heatmap for miRNA-UTR Interactions
#'
#' @description
#' Generates a heatmap of miRNA x UTR interaction scores.
#' Automatically determines the optimal number of clusters using
#' dynamic tree cutting.
#'
#' @param df Data frame containing miRNA, UTR (or utr column) and Score columns.
#' @param output_file Optional path to export PNG (default = NULL).
#' @param width Width of PNG if exporting (in pixels).
#' @param height Height of PNG if exporting (in pixels).
#' @return Returns the ComplexHeatmap object invisibly.
#' @examples
#' \dontrun{
#' file <- system.file("extdata", "exemplo_blocos.txt", package = "miRHeat")
#' df <- parse_file(file)
#' df <- select_score(df, "interaction_energy")
#' df <- prepare_for_heatmap(df)
#' plot_mirheat(df)
#' }
#' @export
plot_mirheat <- function(df,
                         output_file = NULL,
                         width = 2000,
                         height = 1800) {

  # valida entradas
  if (!is.data.frame(df)) stop("df deve ser um data.frame.")
  if (!("miRNA" %in% names(df))) stop("Coluna 'miRNA' não encontrada.")
  if (!("Score" %in% names(df))) stop("Coluna 'Score' não encontrada. Use select_score().")
  if (!("utr" %in% names(df))) {
    # tenta criar via prepare_for_heatmap
    df <- prepare_for_heatmap(df)
  }

  message("→ Construindo matriz miRNA × UTR...")

  # Reformat usando namespace explícito
  matriz <- reshape2::dcast(df, miRNA ~ utr, value.var = "Score")
  # dcast cria miRNA como coluna; movemos para rownames
  rownames(matriz) <- matriz$miRNA
  matriz$miRNA <- NULL

  # Replace NAs with mean of row (se toda linha NA, manter NA -> substituir por 0)
  matriz <- t(apply(matriz, 1, function(x) {
    if (all(is.na(x))) {
      x[is.na(x)] <- 0
    } else {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
    x
  }))

  # Distance + clustering
  d_miRNA <- stats::dist(matriz)
  hc_miRNA <- stats::hclust(d_miRNA, method = "ward.D2")

  message("→ Detectando automaticamente número de clusters...")

  clusters <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_miRNA,
    distM = as.matrix(d_miRNA),
    deepSplit = 2,
    pamRespectsDendro = TRUE
  )

  n_clusters <- length(unique(clusters))
  message(sprintf("→ Número detectado: %d clusters.", n_clusters))

  # Annotations
  heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = clusters,
    col = list(Cluster = structure(
      circlize::rand_color(n_clusters),
      names = as.character(unique(clusters))
    ))
  )

  # Palette for scores
  paleta <- circlize::colorRamp2(
    seq(min(matriz, na.rm = TRUE), max(matriz, na.rm = TRUE), length.out = 5),
    c("#b2182b", "#f4a582", "#f7f7f7", "#92c5de", "#2166ac")
  )

  # Heatmap
  ht <- ComplexHeatmap::Heatmap(
    matriz,
    name = "Score",
    col = paleta,
    cluster_rows = hc_miRNA,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    top_annotation = heatmap_annotation
  )

  # Export PNG if requested
  if (!is.null(output_file)) {
    message("→ Exportando heatmap para: ", output_file)
    grDevices::png(output_file, width = width, height = height, res = 200)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }

  message("→ Heatmap pronto!")
  invisible(ht)
}

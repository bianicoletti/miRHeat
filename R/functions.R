#' Parse interaction files (block-formatted TXT or CSV/TSV)
#'
#' Reads an interaction file either in block format (as produced by RNAhybrid/IntaRNA)
#' or as a tabular CSV/TSV file. For block-formatted TXT files, the function assumes
#' the following structure: a line with an identifier such as `GENE_UTR_123`, followed
#' by alignment lines, a line containing the miRNA name (e.g. `hsa-miR-7160-5p`), and
#' at least one line with a `key = value` pair (e.g. `interaction energy = -24.54 kcal/mol`).
#'
#' For tabular files, the function attempts to automatically detect miRNA and target
#' columns.
#'
#' @param file Path to the input file.
#' @param pattern Character string. Regular expression pattern used to detect the
#'   beginning of a block (default: '^.*_UTR_.*').
#' @param guess_tabular Logical. If TRUE, the function attempts to automatically detect
#'   whether the file is tabular.
#' @param mirna_col Optional character. Name of the column containing miRNAs
#'   (applies to tabular files).
#' @param target_col Optional character. Name of the column containing targets/genes
#'   (applies to tabular files).
#' @return A data.frame containing the extracted columns (miRNA, target, gene, utr,
#'   and any detected numeric values).
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
  if (!file.exists(file)) stop("Input file not found. Please check the file path.", file)

  ext <- tolower(tools::file_ext(file))

  # Helper: normalize the columns names
  clean_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    tolower(trimws(x))
  }

  # ----- (csv/tsv/txt) -----
  if (guess_tabular && ext %in% c("csv", "tsv", "txt")) {
    # attempt to infer the delimiter from the first line

    first_line <- readLines(file, n = 1, warn = FALSE)
    df_tab <- NULL
    if (grepl(",", first_line) && ext %in% c("csv", "txt")) {
      df_tab <- tryCatch(utils::read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else if (grepl("\t", first_line) && ext %in% c("tsv", "txt")) {
      df_tab <- tryCatch(utils::read.delim(file, stringsAsFactors = FALSE), error = function(e) NULL)
    } else if (grepl(";", first_line) && ext %in% c("csv","txt")) {
      # possible semicolon-separated CSV file
      df_tab <- tryCatch(utils::read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
    }

    if (!is.null(df_tab) && is.data.frame(df_tab) && ncol(df_tab) > 0) {
      names(df_tab) <- clean_name(names(df_tab))
      nms <- names(df_tab)

      # If the user provided column names, use them (after cleaning)
      if (!is.null(mirna_col)) mirna_col <- clean_name(mirna_col)
      if (!is.null(target_col)) target_col <- clean_name(target_col)

      # Attempt to automatically detect miRNA/target columns if not provided
      if (is.null(mirna_col)) {
        mirna_candidates <- nms[grepl("mirna|mirna_|mi-rna|mi-rna|\\bmir\\b|mi_r|mi-r", nms, ignore.case = TRUE)]
        if (length(mirna_candidates) >= 1) mirna_col <- mirna_candidates[1]
      }
      if (is.null(target_col)) {
        target_candidates <- nms[grepl("gene|target|utr", nms, ignore.case = TRUE)]
        if (length(target_candidates) >= 1) target_col <- target_candidates[1]
      }

      # If required columns are still missing, abort with an informative message (no interactivity)
      missing <- character(0)
      if (is.null(mirna_col)) missing <- c(missing, "miRNA")
      if (is.null(target_col)) missing <- c(missing, "target/gene/utr")
      if (length(missing) > 0) {
        stop(
          "Unable to automatically identify the required tabular columns: ",
          paste(missing, collapse = ", "),
          ". Colunas detectadas: ", paste(nms, collapse = ", "),
          ".\nSolucao: re-export the file with a clear header or specify mirna_col and target_col in parse_file()."
        )
      }

      # detect numeric columns that resemble scores (excluding key columns)
      score_col <- setdiff(nms[vapply(df_tab, is.numeric, logical(1))], c(mirna_col, target_col))

      df_out <- data.frame(
        miRNA = df_tab[[mirna_col]],
        target = df_tab[[target_col]],
        stringsAsFactors = FALSE
      )

      # append detected numeric columns (scores)
      for (sc in score_col) {
        out_name <- clean_name(sc)
        # evita sobrescrever colunas base
        if (!(out_name %in% c("mirna", "target"))) {
          df_out[[out_name]] <- df_tab[[sc]]
        }
      }

      return(df_out)
    }
    # if df_tab was not detected, continue with block mode
  }

  # ----- Block-formatted TXT case -----
  linhas <- readLines(file, warn = FALSE)
  indices_inicio <- grep(pattern, linhas)
  if (length(indices_inicio) == 0) {
    stop("Nenhum bloco encontrado com o padrao: ", pattern)
  }
  indices_fim <- c((indices_inicio[-1] - 1), length(linhas))
  blocos <- Map(function(i, f) linhas[i:f], indices_inicio, indices_fim)

  # helper functions for extraction
  extrai_miRNA <- function(bloco) {
    # procura padroes como hsa-miR-..., miR-..., miRNA-..., miR...
    pat <- "(?:hsa-)?(?:miR|miRNA)[-A-Za-z0-9_]+"
    hit <- unique(unlist(regmatches(bloco, gregexpr(pat, bloco, ignore.case = TRUE))))
    if (length(hit) >= 1 && nzchar(hit[1])) {
      return(hit[1])
    }
    # fallback: search for a line containing 'miR' or 'mirna', ignoring case
    idx <- grep("mirna|mir", bloco, ignore.case = TRUE)
    if (length(idx) >= 1) {
      cand <- trimws(bloco[idx[1]])
      return(cand)
    }
    return(NA_character_)
  }

  extrai_target <- function(first_line) {
    txt <- trimws(first_line)
    # attempt to extract gene and UTR: GENE_UTR_123 or variants
    m <- regexec("^([A-Za-z0-9]+)[_\\.-]*UTR[_\\.-]*([0-9]+)?", txt, ignore.case = TRUE)
    mm <- regmatches(txt, m)[[1]]
    if (length(mm) >= 2) {
      gene <- mm[2]
      utrid <- ifelse(length(mm) >= 3 && nzchar(mm[3]), mm[3], NA_character_)
      return(list(gene = toupper(gene), utr = utrid, target = ifelse(is.na(utrid), gene, paste0(gene, "_UTR_", utrid))))
    } else {
      # if it does not match, use the entire first line as the raw target
      return(list(gene = txt, utr = NA_character_, target = txt))
    }
  }

  extrai_pairs_numeric <- function(bloco) {
    # capture "name = value" pairs where the value is numeric (optional sign and decimal allowed),
    # accept percentages in the name and ignore units
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
        # if it already exists, coerce to a vector (keeping the last value if necessary)
        res[[name_col]] <- val
      }
    }
    return(res)
  }

  # iterate over blocks
  rows <- vector("list", length(blocos))
  for (i in seq_along(blocos)) {
    b <- blocos[[i]]
    mi <- extrai_miRNA(b)
    target_info <- extrai_target(b[1])
    pairs <- extrai_pairs_numeric(b)
    rows[[i]] <- c(list(miRNA = mi, target = target_info$target, gene = target_info$gene, utr = target_info$utr), pairs)
  }

  # combine into a data.frame (filling NA where necessary)
  # collect all possible column names
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

  # avoid converting key columns to numeric
  core_cols <- c("miRNA", "target", "gene", "utr")
  other_cols <- setdiff(names(df), core_cols)

  # convert "other" columns that appear to be numeric
  for (nm in other_cols) {
    # if the entire column is NA or numeric-like, convert it
    vals_chr <- as.character(df[[nm]])
    if (all(is.na(vals_chr) | grepl("^\\s*-?[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?\\s*$", vals_chr))) {
      df[[nm]] <- as.numeric(vals_chr)
    }
  }

  # rename target to 'target' and organize columns
  # standardize names: if 'target' exists, keep it; if 'target_raw' exists, rename it
  if ("target_raw" %in% names(df) && !("target" %in% names(df))) {
    names(df)[names(df) == "target_raw"] <- "target"
  }

  # organize columns: miRNA, target, gene, utr, others...
  core <- c("miRNA", "target", "gene", "utr")
  others <- setdiff(names(df), core)
  df <- df[c(core, others)]
  rownames(df) <- NULL
  return(df)
}


#' Select which numeric column will be used as Score
#'
#' Automatically detects numeric columns and defines the selected column as a new
#' `Score` column in the returned data.frame. Interactivity is avoided: if multiple
#' numeric columns are detected, `score_column` must be explicitly provided.
#'
#' @param df Data.frame (usually the output of parse_file()).
#' @param score_column Character or NULL. Name of the numeric column to be used as Score.
#'   If NULL and only one numeric column exists, it will be used automatically.
#' @return A data.frame with a new `Score` column
#' @examples
#' df <- data.frame(miRNA = c("a","b"), target = c("T1","T2"), energy = c(-10, -20))
#' select_score(df) # uses energy automatically
#' @export


select_score <- function(df, score_column = NULL) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame returned by parse_file().")
  }

  # If the user explicitly provides the column
  if (!is.null(score_column)) {

    # 1. Check whether it exists
    if (!score_column %in% names(df)) {
      stop("Specified column not found: '", score_column, "'.")
    }

    # 2. Check whether it is numeric
    if (!is.numeric(df[[score_column]])) {
      stop("The specified column is not numeric: '", score_column, "'.")
    }

    df$Score <- df[[score_column]]
    return(df)
  }

  # ---- If score_column is NULL — automatic detection ----

  # identify numeric columns
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("utr", "target", "miRNA", "gene"))

  # no numeric columns
  if (length(numeric_cols) == 0) {
    stop("No numeric column found to be used as Score.")
  }

  # single numeric column, use it
  if (length(numeric_cols) == 1) {
    score_column <- numeric_cols[1]
    message("Automatically using the only numeric column detected: ", score_column)
    df$Score <- df[[score_column]]
    return(df)
  }

  # multiple numeric columns
  stop(
    paste0(
      "Multiple numeric columns detected: ",
      paste(numeric_cols, collapse = ", "),
      ". Please specify the column to use with score_column = \"column_name\"."
    )
  )
}

#' Apply numeric filters to Score
#'
#' Filters the data.frame to keep only rows with Score values within the specified limits.
#'
#' @param df Data.frame containing a `Score` column.
#' @param min_value Numeric or NULL. Keeps Score >= min_value if provided.
#' @param max_value Numeric or NULL. Keeps Score <= max_value if provided.
#' @param remove_na Logical. If TRUE, removes rows with NA Score values.
#' @return A filtered data.frame.
#' @examples
#' df <- data.frame(miRNA = c("a","b","c"), Score = c(-10, NA, -30))
#' apply_numeric_filters(df, min_value = -25)
#' @export
apply_numeric_filters <- function(df, min_value = NULL, max_value = NULL, remove_na = TRUE) {
  if (!is.data.frame(df)) stop("df must be a data.frame.")
  if (!"Score" %in% names(df)) stop("Data.frame must contain a 'Score' column. Use select_score() first.")

  keep <- rep(TRUE, nrow(df))
  if (!is.null(min_value)) keep <- keep & (!is.na(df$Score) & df$Score >= min_value)
  if (!is.null(max_value)) keep <- keep & (!is.na(df$Score) & df$Score <= max_value)
  if (remove_na) keep <- keep & !is.na(df$Score)
  df2 <- df[keep, , drop = FALSE]
  rownames(df2) <- NULL
  return(df2)
}


#' Prepare standardized table for heatmap
#'
#' Receives a data.frame with Score and creates/ensures the columns:
#' miRNA, target, gene, utr, Score. Returns a data.frame ready for heatmap generation.
#'
#' @param df Data.frame with at least miRNA, target/gene, and Score.
#' @return Data.frame with ordered columns: miRNA, target, gene, utr, Score, and optional extra columns.
#' @examples
#' df <- data.frame(miRNA = c("a","b"), target = c("G1_UTR_1","G2"), Score = c(-10, -20))
#' prepare_for_heatmap(df)
#' @export
prepare_for_heatmap <- function(df) {
  if (!is.data.frame(df)) stop("df must be a data.frame.")
  # checks
  if (!("miRNA" %in% names(df))) stop("Column 'miRNA' not found.")
  if (!("Score" %in% names(df))) stop("Column 'Score' not found. Use select_score().")

  # target: if not present, try to construct it from gene/utr
  if (!("target" %in% names(df))) {
    if ("gene" %in% names(df) && "utr" %in% names(df)) {
      df$target <- ifelse(is.na(df$utr), df$gene, paste0(df$gene, "_UTR_", df$utr))
    } else if ("gene" %in% names(df)) {
      df$target <- df$gene
    } else {
      stop("Unable to construct 'target' from the available data.")
    }
  }

  # Ensure gene/utr columns
  if (!("gene" %in% names(df))) {
    df$gene <- sub("(_UTR_.*$)|([_\\.-]UTR[_\\.-].*$)", "", df$target)
  }
  if (!("utr" %in% names(df))) {
    utrs <- regmatches(df$target, regexec("_UTR_([0-9]+)$", df$target))
    df$utr <- sapply(utrs, function(x) if (length(x) >= 2) x[2] else NA_character_)
  }

  # organize columns
  core <- c("miRNA", "target", "gene", "utr", "Score")
  rest <- setdiff(names(df), core)
  df_out <- df[c(core, rest)]
  rownames(df_out) <- NULL
  return(df_out)
}


#' Generate miRNA × UTR interaction heatmap
#'
#' @description
#' Generates a heatmap of interaction scores between miRNAs and UTRs.
#' Automatically determines the optimal number of clusters using
#' dynamic tree cutting.
#'
#' @param df Data.frame containing the columns miRNA, UTR (or utr), and Score.
#' @param output_file Optional path to export the PNG file (default = NULL).
#' @param width Width of the PNG when exporting (in pixels).
#' @param height Height of the PNG when exporting (in pixels).
#' @return Invisibly returns the ComplexHeatmap object.
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

  # validate inputs
  if (!is.data.frame(df)) stop("df must be a data.frame.")
  if (!("miRNA" %in% names(df))) stop("Column 'miRNA' not found.")
  if (!("Score" %in% names(df))) stop("Column 'Score' not found. Use select_score().")
  if (!("utr" %in% names(df))) {
    # try to create it via prepare_for_heatmap
    df <- prepare_for_heatmap(df)
  }

  message("Building miRNA × UTR matrix...")

  # Reformat  using explicit namespace
  matriz <- reshape2::dcast(df, miRNA ~ utr, value.var = "Score")
  # dcast creates miRNA as a column; move it to rownames
  rownames(matriz) <- matriz$miRNA
  matriz$miRNA <- NULL

  # Replace NAs with row mean (if entire row is NA, replace with 0)
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

  message("Automatically detecting number of clusters...")

  clusters <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_miRNA,
    distM = as.matrix(d_miRNA),
    deepSplit = 2,
    pamRespectsDendro = TRUE
  )

  n_clusters <- length(unique(clusters))
  message(sprintf("Detected: %d clusters.", n_clusters))

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
    message("Exportando heatmap para: ", output_file)
    grDevices::png(output_file, width = width, height = height, res = 200)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }

  message("Heatmap ready!")
  invisible(ht)
}

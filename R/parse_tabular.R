#' Parse tabular interaction file
#'
#' Reads interaction data from tabular files such as CSV, TSV, or TXT.
#' The function attempts to automatically detect columns corresponding to
#' miRNA, target, and interaction score.
#'
#' @param file Character. Path to the input file.
#' @param mirna_col Optional character. Name of the column containing miRNA identifiers.
#' @param target_col Optional character. Name of the column containing target identifiers.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{miRNA}
#'     \item \code{target}
#'     \item \code{Score}
#'   }
#'
#' @details
#' The function automatically detects delimiters (comma, semicolon, or tab)
#' and attempts to identify relevant columns based on column names.
#'
#' If automatic detection fails, users should provide \code{mirna_col} and
#' \code{target_col}.
#'
#' @examples
#' \dontrun{
#' df <- parse_tabular("interactions.csv")
#' df <- parse_tabular("interactions.tsv", mirna_col = "miRNA", target_col = "gene")
#' }

parse_tabular <- function(file,
                          mirna_col = NULL,
                          target_col = NULL) {

  if (!file.exists(file)) {
    stop("Input file not found. Please check the file path.")
  }

  # helper
  clean_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    tolower(trimws(x))
  }

  # detectar separador de forma segura
  first_line <- readLines(file, n = 1, warn = FALSE)

  sep <- if (grepl("\t", first_line)) {
    "\t"
  } else if (grepl(";", first_line)) {
    ";"
  } else {
    ","  # default
  }

  df_tab <- tryCatch(
    utils::read.table(
      file,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) {
      stop("Failed to read tabular file.")
    }
  )

  if (!is.data.frame(df_tab) || ncol(df_tab) == 0) {
    stop("Invalid or empty tabular file.")
  }

  names(df_tab) <- clean_name(names(df_tab))
  nms <- names(df_tab)

  # limpar nomes fornecidos
  if (!is.null(mirna_col)) mirna_col <- clean_name(mirna_col)
  if (!is.null(target_col)) target_col <- clean_name(target_col)

  # detectar colunas automaticamente
  if (is.null(mirna_col)) {
    mirna_candidates <- nms[grepl("mirna|\\bmir\\b", nms, ignore.case = TRUE)]
    if (length(mirna_candidates) >= 1) mirna_col <- mirna_candidates[1]
  }

  if (is.null(target_col)) {
    target_candidates <- nms[grepl("target|gene|utr", nms, ignore.case = TRUE)]
    if (length(target_candidates) >= 1) target_col <- target_candidates[1]
  }

  # erro se não encontrar
  if (is.null(mirna_col) || is.null(target_col)) {
    stop(
      "Could not detect required columns.\n",
      "Detected columns: ", paste(nms, collapse = ", "), "\n",
      "Provide 'mirna_col' and 'target_col'."
    )
  }

  # montar output
  df_out <- data.frame(
    miRNA = df_tab[[mirna_col]],
    target = df_tab[[target_col]],
    Score = NA_real_,
    stringsAsFactors = FALSE
  )

  # detectar score
  score_candidates <- nms[grepl("energy|mfe|score", nms, ignore.case = TRUE)]

  if (length(score_candidates) >= 1) {
    df_out$Score <- suppressWarnings(as.numeric(df_tab[[score_candidates[1]]]))
  }

  return(.standardize_output(df_out))
}

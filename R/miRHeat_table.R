#' Generate a standardized interaction table for mirHeat
#'
#' This function isolates the tabular output of the mirHeat pipeline.
#' It returns a clean, standardized miRNA–target interaction table
#' suitable for inspection or export (e.g. CSV), without generating
#' a heatmap.
#'
#' @param df A data.frame produced by parse_file(), select_score(), and apply_numeric_filters().
#'
#' @return A tibble with standardized columns: miRNA, gene, utr, target, Score.
#'
#' @export
#'
mirHeat_table <- function(df) {

  # -------------------------------------------------
  # 1) Basic validation
  # -------------------------------------------------
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }

  if (!"miRNA" %in% colnames(df)) {
    stop("Required column 'miRNA' not found.")
  }

  if (!"Score" %in% colnames(df)) {
    stop("Required column 'Score' not found.")
  }

  # -------------------------------------------------
  # 2) Ensure target column exists
  # -------------------------------------------------
  if (!"target" %in% colnames(df)) {

    if (!all(c("gene", "utr") %in% colnames(df))) {
      stop(
        "Column 'target' not found. ",
        "Provide either 'target' or both 'gene' and 'utr' columns."
      )
    }

    df$target <- paste(df$gene, "UTR", df$utr, sep = "_")
  }

  # -------------------------------------------------
  # 3) Extract gene and utr from target if needed
  # -------------------------------------------------
  if (!"gene" %in% colnames(df) || !"utr" %in% colnames(df)) {

    parts <- strsplit(df$target, "_UTR_")

    df$gene <- vapply(parts, `[`, character(1), 1)
    df$utr  <- vapply(parts, `[`, character(1), 2)
  }

  # -------------------------------------------------
  # 4) Build standardized output table
  # -------------------------------------------------
  out <- df[, c("miRNA", "gene", "utr", "target", "Score"), drop = FALSE]

  # -------------------------------------------------
  # 5) Return as tibble
  # -------------------------------------------------
  tibble::as_tibble(out)
}

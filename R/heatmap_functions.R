#' Select or redefine Score column
#'
#' Allows the user to define or override the `Score` column in a data.frame.
#' This is mainly useful for custom or external tabular data where the Score
#' was not automatically detected.
#'
#' @param df A data.frame containing at least `miRNA` and `target` columns.
#' @param score_column Character or NULL. Name of the numeric column to be used as Score.
#'
#' @return A data.frame with a `Score` column.
#'
#' @examples
#' df <- data.frame(miRNA = c("a","b"), target = c("T1","T2"), energy = c(-10, -20))
#' select_score(df, "energy")
#'
#' @export
#'
select_score <- function(df, score_column = NULL) {

  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }


  if ("Score" %in% names(df) && is.null(score_column)) {
    return(df)
  }


  if (!is.null(score_column)) {

    if (!score_column %in% names(df)) {
      stop("Column not found: '", score_column, "'.")
    }

    if (!is.numeric(df[[score_column]])) {
      stop("Column is not numeric: '", score_column, "'.")
    }

    df$Score <- df[[score_column]]
    return(df)
  }


  numeric_cols <- names(df)[sapply(df, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, "Score")

  if (length(numeric_cols) == 1) {
    df$Score <- df[[numeric_cols]]
    return(df)
  }

  stop(
    "Could not determine Score column automatically. ",
    "Please provide 'score_column'."
  )
}

#' Apply numeric filters to Score
#'
#' Filters a data.frame to keep only rows with Score values within specified limits.
#'
#' @param df Data.frame containing a numeric `Score` column.
#' @param min_value Numeric or NULL. Keeps Score >= min_value if provided.
#' @param max_value Numeric or NULL. Keeps Score <= max_value if provided.
#' @param remove_na Logical. If TRUE, removes rows with NA Score values.
#'
#' @return A filtered data.frame.
#'
#' @examples
#' df <- data.frame(miRNA = c("a","b","c"), Score = c(-10, NA, -30))
#' apply_numeric_filters(df, min_value = -25)
#'
#' @export
apply_numeric_filters <- function(df,
                                  min_value = NULL,
                                  max_value = NULL,
                                  remove_na = TRUE) {

  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }

  if (!"Score" %in% names(df)) {
    stop("Data.frame must contain a 'Score' column.")
  }

  if (!is.numeric(df$Score)) {
    stop("'Score' column must be numeric.")
  }

  keep <- rep(TRUE, nrow(df))

  # remover NA primeiro (mais limpo)
  if (remove_na) {
    keep <- keep & !is.na(df$Score)
  }

  if (!is.null(min_value)) {
    keep <- keep & df$Score >= min_value
  }

  if (!is.null(max_value)) {
    keep <- keep & df$Score <= max_value
  }

  df2 <- df[keep, , drop = FALSE]
  rownames(df2) <- NULL

  return(df2)
}

#' Prepare standardized table for heatmap
#'
#' Ensures the input data.frame has the required structure for heatmap generation.
#' Optionally aggregates duplicated miRNA-target pairs.
#'
#' @param df A data.frame containing the columns \code{miRNA}, \code{target}, and \code{Score}.
#' @param remove_na Logical. If \code{TRUE}, rows with \code{NA} values are removed. Default is \code{TRUE}.
#' @param aggregate_fun Function or NULL. Function used to aggregate duplicate miRNA-target pairs
#'   (e.g., \code{mean}, \code{min}, \code{max}). If \code{NULL} (default), no aggregation is performed.
#'
#' @return Data.frame com colunas: miRNA, target, Score.
#'
#' @examples
#' df <- data.frame(
#'   miRNA = c("a","a"),
#'   target = c("T1","T1"),
#'   Score = c(-10, -20)
#' )
#' prepare_for_heatmap(df, aggregate_fun = mean)
#'
#' @export
prepare_for_heatmap <- function(df,
                                remove_na = TRUE,
                                aggregate_fun = NULL) {

  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }

  required <- c("miRNA", "target", "Score")

  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  if (!is.numeric(df$Score)) {
    stop("'Score' column must be numeric.")
  }

  df_out <- df[required]

  if (!is.null(aggregate_fun)) {
    df_out <- aggregate(
      Score ~ miRNA + target,
      data = df_out,
      FUN = aggregate_fun
    )
  }


  if (remove_na) {
    df_out <- df_out[complete.cases(df_out), , drop = FALSE]
  }

  rownames(df_out) <- NULL

  return(df_out)
}

#' Generate miRNA vs target interaction heatmap
#'
#' @description
#' Generates a heatmap of interaction scores between miRNAs and targets using
#' ComplexHeatmap. Automatically detects the number of clusters via dynamic tree
#' cutting, with optional manual override. Supports built-in and custom color
#' palettes. When the number of unique interactions is large, only the top
#' strongest interactions are plotted.
#'
#' @param df Data.frame containing columns \code{miRNA}, \code{target}, and
#'   \code{Score}. Scores are expected to be negative (e.g., binding energy):
#'   lower (more negative) values indicate stronger interactions.
#' @param output_file Character or NULL. Optional path to export the plot as a
#'   PNG file. Default is \code{NULL} (no export).
#' @param width Integer. Width of the PNG in pixels when exporting. Default is
#'   \code{2400}.
#' @param height Integer. Height of the PNG in pixels when exporting. Default
#'   is \code{2000}.
#' @param top_n Integer or NULL. Maximum number of interactions (miRNA-target
#'   pairs) to display, selected by lowest (strongest) Score. Default is
#'   \code{50}. Set to \code{NULL} to display all interactions.
#' @param n_clusters Integer or NULL. Number of row clusters to use. If
#'   \code{NULL} (default), the number is detected automatically via
#'   \code{dynamicTreeCut}. If a value is provided, it overrides the automatic
#'   detection and a warning is issued.
#' @param palette Character. Built-in color palette to use. One of
#'   \code{"RdBu"} (default), \code{"RdYlBu"}, or \code{"viridis"}. Ignored
#'   if \code{custom_palette} is provided.
#' @param custom_palette A vector of colors to build a custom palette via
#'   \code{circlize::colorRamp2}. Must have at least 2 colors. The breakpoints
#'   are distributed evenly across the Score range. Default is \code{NULL}.
#' @param cluster_colors Character vector or NULL. Colors for cluster annotation.
#'   If NULL (default), colors are assigned automatically.
#'
#' @return Invisibly returns the \code{ComplexHeatmap} object.
#'
#' @examples
#' \dontrun{
#' df <- parse_interactions("file.txt", type = "rnahybrid")
#' df <- apply_numeric_filters(df, min_value = -30)
#' df <- prepare_for_heatmap(df)
#'
#' # Default: top 50, automatic clusters, RdBu palette
#' plot_mirheat(df)
#'
#' # Show all interactions, force 5 clusters, use viridis
#' plot_mirheat(df, top_n = NULL, n_clusters = 5, palette = "viridis")
#'
#' # Custom palette
#' plot_mirheat(df, custom_palette = c("navy", "white", "firebrick"))
#'
#' # Export to PNG
#' plot_mirheat(df, output_file = "heatmap.png")
#' }
#'
#' @export
plot_mirheat <- function(df,
                         output_file  = NULL,
                         width        = 2400,
                         height       = 2000,
                         top_n        = 200,
                         n_clusters   = NULL,
                         palette      = c("RdBu", "RdYlBu", "viridis"),
                         custom_palette = NULL,
                         cluster_colors = NULL)
  {

  # ---- input validation -------------------------------------------------------
  if (!is.data.frame(df)) stop("'df' must be a data.frame.")

  required <- c("miRNA", "target", "Score")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(df$Score)) stop("'Score' column must be numeric.")

  palette <- match.arg(palette)

  if (!is.null(n_clusters)) {
    if (!is.numeric(n_clusters) || length(n_clusters) != 1 || n_clusters < 1) {
      stop("'n_clusters' must be a single positive integer or NULL.")
    }
    n_clusters <- as.integer(n_clusters)
  }

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
      stop("'top_n' must be a single positive integer or NULL.")
    }
    top_n <- as.integer(top_n)
  }

  # ---- filter top interactions ------------------------------------------------
  n_total <- nrow(df)

  if (!is.null(top_n) && n_total > top_n) {
    message(
      sprintf(
        "Large dataset detected: %d interactions found. Displaying top %d targets based on strongest interaction (lowest Score).\nTo display all targets, set top_n = NULL.",
        n_total, top_n
      )
    )
    if (!is.null(top_n)) {


      target_score <- aggregate(
        Score ~ target,
        data = df,
        FUN = min
      )


      target_score <- target_score[order(target_score$Score), ]


      top_targets <- head(target_score$target, top_n)


      df <- df[df$target %in% top_targets, ]
    }
  }

  # ---- build matrix -----------------------------------------------------------
  message("Building miRNA x target Score matrix...")

  safe_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    min(x)
  }

  mat <- reshape2::dcast(
    df,
    miRNA ~ target,
    value.var     = "Score",
    fun.aggregate = safe_min
  )

  rownames(mat) <- mat$miRNA
  mat$miRNA     <- NULL
  mat           <- as.matrix(mat)



  # fill structural NAs (missing pairs) with 0
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message(sprintf(
      "Note: %d missing miRNA-target pairs filled with 0 (no interaction data).",
      n_na
    ))
    mat[is.na(mat)] <- 0
  }

  if (nrow(mat) < 2) {
    stop(
      "At least 2 miRNAs are required to build a heatmap. ",
      "Consider relaxing your filters."
    )
  }

  # ---- clustering -------------------------------------------------------------
  d_mirna  <- stats::dist(mat)
  hc_mirna <- stats::hclust(d_mirna, method = "ward.D2")

  message("Automatically detecting number of clusters...")

  auto_clusters <- tryCatch(
    dynamicTreeCut::cutreeDynamic(
      dendro            = hc_mirna,
      distM             = as.matrix(d_mirna),
      deepSplit         = 2,
      pamRespectsDendro = TRUE
    ),
    error = function(e) NULL
  )

  n_auto <- if (!is.null(auto_clusters)) length(unique(auto_clusters)) else 1

  # fallback: if dynamicTreeCut detects only 1 cluster (flat dendrogram or
  # sparse matrix), use cutree with k = min(3, nrow)
  if (n_auto <= 1) {
    k_fallback  <- min(3L, nrow(mat))
    auto_clusters <- stats::cutree(hc_mirna, k = k_fallback)
    n_auto      <- k_fallback
    message(sprintf(
      "dynamicTreeCut could not detect structure. Falling back to k = %d clusters.",
      n_auto
    ))
  }

  message(sprintf("Detected: %d cluster(s).", n_auto))

  if (!is.null(n_clusters)) {
    if (n_clusters != n_auto) {
      warning(sprintf(
        "Overriding automatic cluster detection: using %d cluster(s) instead of %d. ",
        n_clusters, n_auto
      ))
    }
    cluster_vec <- stats::cutree(hc_mirna, k = n_clusters)
  } else {
    cluster_vec <- auto_clusters
    n_clusters  <- n_auto
  }

  # ---- color palette ----------------------------------------------------------
  score_range <- range(mat, na.rm = TRUE)
  breaks      <- seq(score_range[1], score_range[2], length.out = 5)

  if (!is.null(custom_palette)) {
    if (length(custom_palette) < 2) {
      stop("'custom_palette' must contain at least 2 colors.")
    }
    pal_breaks <- seq(score_range[1], score_range[2],
                      length.out = length(custom_palette))
    col_fun <- circlize::colorRamp2(pal_breaks, custom_palette)

  } else {
    col_fun <- switch(
      palette,
      "RdBu" = circlize::colorRamp2(
        breaks,
        c("#b2182b", "#f4a582", "#f7f7f7", "#92c5de", "#2166ac")
      ),
      "RdYlBu" = circlize::colorRamp2(
        breaks,
        c("#d73027", "#fdae61", "#ffffbf", "#abd9e9", "#4575b4")
      ),
      "viridis" = circlize::colorRamp2(
        breaks,
        c("#440154", "#3b528b", "#21908c", "#5dc963", "#fde725")
      )
    )
  }

  # ---- cluster annotation colors ----------------------------------------------
  cluster_ids <- as.character(sort(unique(cluster_vec)))

  if (is.null(cluster_colors)) {

    cluster_colors <- structure(
      circlize::rand_color(length(cluster_ids), luminosity = "bright"),
      names = cluster_ids
    )
  } else {


    if (!is.vector(cluster_colors)) {
      stop("'cluster_colors' must be a vector of colors.")
    }

    if (length(cluster_colors) != length(cluster_ids)) {
      stop(
        "Length of 'cluster_colors' must match the number of clusters (",
        length(cluster_ids), ")."
      )
    }

    cluster_colors <- structure(cluster_colors, names = cluster_ids)
  }

  row_anno <- ComplexHeatmap::rowAnnotation(
    Cluster = as.character(cluster_vec),
    col     = list(Cluster = cluster_colors),
    annotation_legend_param = list(
      Cluster = list(title = "Cluster")
    )
  )

  # ---- heatmap ----------------------------------------------------------------
  # Sizing guidance:
  #   - unit-based cell sizing scales better than fixed px dimensions
  #   - row/col name font size tuned for typical publication figures
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  if (n_rows > 50) {
    message("Large number of miRNAs detected. Consider increasing output size or reducing top_n.")
  }
  if (n_cols > 50) {
    message("Large number of targets detected. Consider increasing output size or reducing top_n.")
  }

  if (is.null(row_fontsize)) {
    row_fontsize <- max(6, min(12, 180 / n_rows))
  }

  if (is.null(col_fontsize)) {
    col_fontsize <- max(6, min(12, 180 / n_cols))
  }

  ht <- ComplexHeatmap::Heatmap(
    mat,
    name                  = "Score",
    col                   = col_fun,
    cluster_rows          = hc_mirna,
    cluster_columns       = TRUE,
    show_row_names        = TRUE,
    show_column_names     = TRUE,
    row_names_gp          = grid::gpar(fontsize = row_fontsize),
    column_names_gp       = grid::gpar(fontsize = col_fontsize),
    cell_fun              = NULL,
    left_annotation       = row_anno,
    row_title             = "miRNA",
    column_title          = "Target",
    row_title_gp          = grid::gpar(fontsize = 12, fontface = "bold"),
    column_title_gp       = grid::gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param  = list(
      title          = "Score",
      legend_height  = grid::unit(4, "cm"),
      title_gp       = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp      = grid::gpar(fontsize = 9)
    )
  )

  ComplexHeatmap::draw(ht)

  # ---- export -----------------------------------------------------------------
  if (!is.null(output_file)) {
    message("Exporting heatmap to: ", output_file)

    # auto-scale export dimensions if not overridden by user
    export_w <- max(width,  n_cols * 40 + 400)
    export_h <- max(height, n_rows * 40 + 400)

    grDevices::png(output_file, width = export_w, height = export_h, res = 200)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()

    message(sprintf(
      "Saved: %s (%d x %d px, 200 dpi).", output_file, export_w, export_h
    ))
  }

  message("Done!")
  invisible(ht)
}

#' Standardize parsed interaction output
#'
#' Internal helper function to ensure consistent output format across
#' different parsers.
#'
#' @param df A data.frame containing parsed interaction data.
#'
#' @return A data.frame with standardized columns:
#' \code{miRNA}, \code{target}, and \code{Score}.
#'
#' @keywords internal


.standardize_output <- function(df) {

  required <- c("miRNA", "target", "Score")

  for (col in required) {
    if (!col %in% names(df)) {
      df[[col]] <- NA
    }
  }

  df <- df[required]

  df$Score <- as.numeric(df$Score)

  return(df)
}

#' Parse RNAhybrid output file
#'
#' Parses interaction data from RNAhybrid block-formatted output files.
#' Extracts miRNA identifiers, target names, and minimum free energy (mfe)
#' values as interaction scores.
#'
#' @param file Character. Path to the RNAhybrid output file.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{miRNA}
#'     \item \code{target}
#'     \item \code{Score} (mfe values)
#'   }
#'
#' @details
#' The function expects blocks starting with \code{"target:"} and containing
#' lines with \code{"miRNA :"} and \code{"mfe:"}.
#'
#' @examples
#' \dontrun{
#' df <- parse_rnahybrid("rnahybrid_output.txt")
#' }
parse_rnahybrid <- function(file) {
  linhas <- readLines(file, warn = FALSE)

  blocks <- split(linhas, cumsum(grepl("^target:", linhas)))

  rows <- lapply(blocks, function(b) {

    miRNA <- sub(".*miRNA\\s*:\\s*", "", grep("miRNA", b, value = TRUE)[1])

    target <- sub("target:\\s*", "", grep("^target:", b, value = TRUE)[1])

    score <- as.numeric(sub("mfe:\\s*", "", sub(" kcal/mol", "", grep("mfe:", b, value = TRUE)[1])))

    data.frame(
      miRNA = miRNA,
      target = target,
      Score = score,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  .standardize_output(df)
}

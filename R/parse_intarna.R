#' Parse IntaRNA output file
#'
#' Parses interaction data from IntaRNA text output files.
#' Extracts miRNA identifiers, target names, and interaction energy values.
#'
#' @param file Character. Path to the IntaRNA output file.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{miRNA}
#'     \item \code{target}
#'     \item \code{Score} (energy values)
#'   }
#'
#' @details
#' The function expects blocks starting with lines beginning with \code{">"},
#' where the first header corresponds to the target and the second to the miRNA.
#'
#' @examples
#' \dontrun{
#' df <- parse_intarna("intarna_output.txt")
#' }

parse_intarna <- function(file) {
  linhas <- readLines(file, warn = FALSE)

  sep_idx <- grep("^=+", linhas)

  start_idx <- c(1, sep_idx + 1)
  end_idx <- c(sep_idx - 1, length(linhas))

  blocks <- Map(function(i, j) linhas[i:j], start_idx, end_idx)

  rows <- lapply(blocks, function(b) {

    b <- trimws(b)
    b <- b[b != ""]

    headers <- grep("^>", b, value = TRUE)
    headers <- sub("^>", "", headers)

    if (length(headers) < 2) return(NULL)


    mirna_idx <- grep("miR", headers, ignore.case = TRUE)

    if (length(mirna_idx) == 0) return(NULL)

    miRNA <- headers[mirna_idx[1]]
    target <- headers[-mirna_idx[1]][1]

    score_line <- grep("energy:", b, value = TRUE)
    if (length(score_line) == 0) return(NULL)

    score <- as.numeric(
      sub("energy:\\s*", "",
          sub(" kcal/mol", "", score_line[1]))
    )

    data.frame(
      miRNA = miRNA,
      target = target,
      Score = score,
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)

  df <- do.call(rbind, rows)

  df <- unique(df)

  .standardize_output(df)
}

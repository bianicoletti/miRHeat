#' Parse interaction files
#'
#' Main interface to parse interaction data from different input formats.
#' This function dispatches the parsing to the appropriate method based on
#' the specified type.
#'
#' @param file Character. Path to the input file.
#' @param type Character. Type of input file. Must be one of:
#'   \code{"rnahybrid"}, \code{"intarna"}, or \code{"tabular"}.
#'
#' @return A data.frame with standardized columns:
#'   \itemize{
#'     \item \code{miRNA} - miRNA identifier
#'     \item \code{target} - target identifier
#'     \item \code{Score} - interaction score (numeric)
#'   }
#'
#' @examples
#' \dontrun{
#' parse_interactions("file.txt", type = "rnahybrid")
#' parse_interactions("file.txt", type = "intarna")
#' parse_interactions("file.csv", type = "tabular")
#' }
#'
#' @export

parse_interactions <- function(file,
                               type = c("rnahybrid", "intarna", "tabular")) {

  type <- match.arg(type)

  if (type == "rnahybrid") {
    return(parse_rnahybrid(file))
  }

  if (type == "intarna") {
    return(parse_intarna(file))
  }

  if (type == "tabular") {
    return(parse_tabular(file))
  }
}


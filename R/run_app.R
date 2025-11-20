#' Executa o aplicativo Shiny do pacote
#'
#' @return Inicia o aplicativo Shiny
#' @export
run_app <- function() {
  app_dir <- system.file("app", package = "miRHeat") # troque por seu pacote
  if (app_dir == "") stop("Nao foi possivel encontrar a pasta do app dentro do pacote.")
  shiny::runApp(app_dir)
}


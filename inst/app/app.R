library(shiny)
library(DT)
library(bslib)
library(miRHeat)

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "minty"),

  titlePanel("App do miRHEat - parse_file()"),

  sidebarLayout(

    sidebarPanel(
      h4("Configurações de entrada"),

      fileInput(
        "file",
        "Escolha um arquivo",
        placeholder = "Selecione um arquivo .txt, .csv, .tsv, etc."
      ),
      helpText("Carregue o arquivo que deseja processar. O app tentará identificar automaticamente seu formato."),

      textInput(
        "pattern",
        "Pattern (regex)",
        value = "^.*_UTR_.*",
        placeholder = "^.*_UTR_.*"
      ),
      helpText("Padrão usado para filtrar colunas/linhas via expressão regular (regex)."),

      checkboxInput("guess", "Detectar formato automaticamente", TRUE),
      helpText("Se ativado, o app tentará reconhecer automaticamente o formato tabular."),

      hr(),

      h4("Sobre o app"),
      helpText("Este aplicativo demonstra o uso da função parse_file() do pacote."),
      helpText("Ele permite:"),
      tags$ul(
        tags$li("Carregar arquivos tabulares ou texto"),
        tags$li("Aplicar uma expressão regular (pattern) para filtrar dados"),
        tags$li("Visualizar a tabela resultante")
      ),
      helpText("Ideal para testes rápidos, análises simples e inspeção de dados.")
    ),

    mainPanel(
      h4("Prévia dos dados processados"),
      helpText("Os dados aparecerão aqui após o upload."),
      DTOutput("table"),
      br(),
      uiOutput("msg")  # área para mensagens adicionais
    )
  )
)

server <- function(input, output, session) {

  # Processamento reativo
  parsed <- reactive({
    req(input$file)

    # Tentativa de leitura com tryCatch para erros no parse
    out <- tryCatch({
      parse_file(
        file = input$file$datapath,
        pattern = input$pattern,
        guess_tabular = input$guess
      )
    }, error = function(e) {
      return(NULL)
    })

    out
  })

  # Validações antes de mostrar a tabela
  output$table <- renderDT({
    validate(
      need(input$file, "Nenhum arquivo enviado ainda."),
      need(parsed(), "Não foi possível processar o arquivo. Verifique o formato.")
    )

    # Se o resultado for vazio, mostra mensagem
    if (nrow(parsed()) == 0) {
      validate("Nenhuma linha corresponde ao pattern informado. Tente outro pattern.")
    }

    datatable(parsed(), options = list(pageLength = 10))
  })

  # Mensagens explicativas adicionais
  output$msg <- renderUI({
    if (is.null(input$file)) {
      return(NULL)
    }

    if (is.null(parsed())) {
      return(tags$div(
        class = "alert alert-danger",
        "O arquivo foi carregado, mas ocorreu um erro ao processá-lo."
      ))
    }

    if (nrow(parsed()) == 0) {
      return(tags$div(
        class = "alert alert-warning",
        "Nenhum resultado encontrado para o pattern informado."
      ))
    }

    # Se tudo deu certo
    tags$div(
      class = "alert alert-success",
      paste("Arquivo processado com sucesso! Linhas:", nrow(parsed()))
    )
  })
}

shinyApp(ui, server)


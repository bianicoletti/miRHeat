#' Filtra interações por energia em saída do RNAhybrid / IntaRNA
#'
#' Lê um arquivo de saída do RNAhybrid (ou IntaRNA), identifica blocos de interação
#' (cada bloco iniciando por uma linha com o identificador da sequência, por exemplo
#' "GENE_UTR_..."), extrai a linha que contém "interaction energy = <valor>" e
#' filtra blocos com energia menor (mais negativa) que um corte definido.
#'
#' @details
#' - O arquivo de entrada deve conter blocos que iniciem por uma linha com o
#'   identificador da sequência (por exemplo "MICA_UTR_1", "TP53_UTR_2" etc.).
#' - Cada bloco deve incluir ao menos uma linha com o texto "interaction energy = <valor>"
#'   (o sufixo "kcal/mol" é ignorado na extração).
#' - O argumento `pattern` serve para localizar as linhas iniciais dos blocos.
#'
#' @param input_file String. Caminho para o arquivo de saída do RNAhybrid/IntaRNA.
#' @param output_file String. Caminho do arquivo `.txt` onde serão gravados os blocos filtrados.
#' @param pattern String. Expressão regular usada para identificar o início de cada bloco
#'   (por padrão: `"^MICA_UTR_"`). Deve identificar a linha inicial do bloco.
#' @param energy_cutoff Numeric. Corte de energia; blocos com energia < energy_cutoff
#'   serão mantidos (padrão = -8).
#' @param verbose Logical. Se TRUE, imprime mensagens informativas (padrão = TRUE).
#'
#' @return Invisivelmente, retorna uma lista com:
#'   - blocks: lista dos blocos filtrados (cada elemento é um vetor de linhas);
#'   - summary: data.frame com colunas block_id, energy, block_index.
#'
#' @examples
#' \dontrun{
#' res <- filter_energy("resultado_pred_intaRNA.txt",
#'                      output_file = "interacoes_filtradas.txt",
#'                      pattern = "^MICA_UTR_",
#'                      energy_cutoff = -8)
#' str(res$summary)
#' }
#'
#' @export
filter_energy <- function(input_file,
                          output_file = "interacoes_filtradas.txt",
                          pattern = "^MICA_UTR_",
                          energy_cutoff = -8,
                          verbose = TRUE) {

  # checagens iniciais
  if (!file.exists(input_file)) {
    stop("Arquivo não encontrado: ", input_file)
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("`pattern` deve ser um único string com expressão regular.")
  }

  # lê linhas do arquivo (preserva encoding e evita warning em arquivos grandes)
  linhas <- readLines(input_file, warn = FALSE)

  # identifica inícios de bloco
  indices_inicio <- grep(pattern, linhas)
  if (length(indices_inicio) == 0) {
    stop("Nenhum bloco encontrado com o padrão: ", pattern,
         ". Verifique o formato do arquivo ou o argumento `pattern`.")
  }
  indices_fim <- c((indices_inicio[-1] - 1), length(linhas))

  # separa blocos
  blocos <- Map(function(i, f) linhas[i:f], indices_inicio, indices_fim)

  # função interna para extrair energia (procura última ocorrência dentro do bloco)
  extrai_energia <- function(bloco) {
    # procura linhas que contenham 'interaction energy'
    linhas_energy <- grep("interaction energy", bloco, value = TRUE, ignore.case = TRUE)
    if (length(linhas_energy) == 0) return(NA_real_)
    # pega a última ocorrência para maior robustez
    linha <- tail(linhas_energy, 1)
    # regex que captura número com ou sem decimal, com sinal opcional
    match <- regmatches(linha, regexec("interaction energy\\s*=\\s*(-?[0-9]+\\.?[0-9]*)", linha, ignore.case = TRUE))[[1]]
    if (length(match) >= 2) {
      return(as.numeric(match[2]))
    } else {
      return(NA_real_)
    }
  }

  # extrai identificador da primeira linha de cada bloco (trim)
  extrai_id <- function(bloco) {
    first_line <- bloco[1]
    id <- trimws(first_line)
    # retorna a primeira linha como identificador bruto (o usuário pode pós-processar)
    return(id)
  }

  # filtra blocos por energia
  summary_list <- vector("list", length(blocos))
  blocos_filtrados <- list()
  for (i in seq_along(blocos)) {
    energia <- extrai_energia(blocos[[i]])
    id <- extrai_id(blocos[[i]])
    summary_list[[i]] <- list(block_index = i, id = id, energy = energia)
    if (!is.na(energia) && energia < energy_cutoff) {
      blocos_filtrados[[length(blocos_filtrados) + 1]] <- blocos[[i]]
    }
  }

  # cria data.frame resumo
  summary_df <- do.call(rbind, lapply(summary_list, function(x) {
    data.frame(block_index = x$block_index,
               block_id = x$id,
               energy = x$energy,
               stringsAsFactors = FALSE)
  }))

  # mensagem
  nkept <- length(blocos_filtrados)
  if (verbose) {
    message(nkept, " bloco(s) encontrados com energia < ", energy_cutoff, ".")
  }

  # salva blocos filtrados em arquivo .txt (cada bloco separado por linha vazia)
  if (nkept > 0) {
    writeLines(unlist(lapply(blocos_filtrados, function(b) c(b, ""))),
               con = output_file)
    if (verbose) message("Blocos salvos em: ", output_file)
  } else {
    warning("Nenhum bloco com energia < ", energy_cutoff, " encontrado. Nenhum arquivo gravado.")
  }

  # retorno invisível: blocos + resumo (útil para pipeline em R)
  out <- list(blocks = blocos_filtrados, summary = summary_df)
  invisible(out)
}

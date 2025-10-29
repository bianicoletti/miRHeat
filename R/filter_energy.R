#' Filtra interações por energia em saída do RNAhybrid / IntaRNA
#'
#' Lê um arquivo de saída do RNAhybrid (ou IntaRNA), identifica blocos de interação
#' (cada bloco iniciando por uma linha com identificador do tipo "GENE_UTR_*"),
#' extrai a linha que contém "interaction energy = <valor>" e
#' filtra blocos com energia menor (mais negativa) que um corte definido.
#'
#' @details
#' - Por default, a função detecta automaticamente qualquer bloco que tenha
#'   uma linha inicial terminando com "_UTR_*", independentemente do gene.
#' - Cada bloco precisa conter pelo menos uma linha com "interaction energy = <valor>".
#'
#' @param input_file String. Caminho para o arquivo de saída do RNAhybrid/IntaRNA.
#' @param output_file String. Caminho do arquivo `.txt` onde serão gravados os blocos filtrados.
#' @param pattern String. Expressão regular usada para identificar o início de cada bloco.
#'   Por padrão: `'^.*_UTR_.*'` (qualquer gene/UTR). Pode ser ajustado se necessário.
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
#'                      energy_cutoff = -8)
#' str(res$summary)
#' }
#'
#' @export
filter_energy <- function(input_file,
                          output_file = "interacoes_filtradas.txt",
                          pattern = "^.*_UTR_.*",  # padrão genérico
                          energy_cutoff = -8,
                          verbose = TRUE) {

  # checagens iniciais
  if (!file.exists(input_file)) stop("Arquivo não encontrado: ", input_file)
  if (!is.character(pattern) || length(pattern) != 1) stop("`pattern` deve ser um único string.")

  # lê linhas do arquivo
  linhas <- readLines(input_file, warn = FALSE)

  # identifica inícios de bloco
  indices_inicio <- grep(pattern, linhas)
  if (length(indices_inicio) == 0) stop("Nenhum bloco encontrado com o padrão: ", pattern)
  indices_fim <- c((indices_inicio[-1] - 1), length(linhas))

  # separa blocos
  blocos <- Map(function(i, f) linhas[i:f], indices_inicio, indices_fim)

  # função interna para extrair energia
  extrai_energia <- function(bloco) {
    linhas_energy <- grep("interaction energy", bloco, value = TRUE, ignore.case = TRUE)
    if (length(linhas_energy) == 0) return(NA_real_)
    linha <- tail(linhas_energy, 1)
    match <- regmatches(linha, regexec("interaction energy\\s*=\\s*(-?[0-9]+\\.?[0-9]*)", linha, ignore.case = TRUE))[[1]]
    if (length(match) >= 2) return(as.numeric(match[2])) else return(NA_real_)
  }

  # extrai identificador do bloco (primeira linha)
  extrai_id <- function(bloco) {
    first_line <- bloco[1]
    id <- trimws(first_line)
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

  # mensagens
  nkept <- length(blocos_filtrados)
  if (verbose) message(nkept, " bloco(s) encontrados com energia < ", energy_cutoff, ".")

  # salva blocos filtrados em arquivo .txt
  if (nkept > 0) {
    writeLines(unlist(lapply(blocos_filtrados, function(b) c(b, ""))), con = output_file)
    if (verbose) message("Blocos salvos em: ", output_file)
  } else {
    warning("Nenhum bloco com energia < ", energy_cutoff, " encontrado. Nenhum arquivo gravado.")
  }

  # retorno invisível
  invisible(list(blocks = blocos_filtrados, summary = summary_df))
}

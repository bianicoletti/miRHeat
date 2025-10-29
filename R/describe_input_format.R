#' Mostra exemplo do formato de entrada esperado
#'
#' Imprime no console um exemplo simples do formato aceitável para ser lido
#' por `filter_energy()` e funções relacionadas.
#'
#' @return Invisivelmente NULL. Apenas imprime no console.
#' @examples
#' \dontrun{
#' describe_input_format()
#' }
#' @export
describe_input_format <- function() {
  cat(
    "Exemplo de formato esperado (saída do RNAhybrid / IntaRNA):\n\n",
    "GENE_UTR_1\n",
    "alignment info ...\n",
    "more alignment info ...\n",
    "interaction energy = -12.3 kcal/mol\n",
    "\n",
    "GENE_UTR_2\n",
    "alignment info ...\n",
    "interaction energy = -7.5 kcal/mol\n",
    "\n",
    "Notas:\n",
    "- Cada bloco deve começar com o identificador (linha única) do tipo 'GENE_UTR_...'.\n",
    "- Cada bloco precisa ter uma linha contendo 'interaction energy = <valor>'.\n",
    "- Ajuste o argumento `pattern` da função `filter_energy()` para localizar o início dos blocos.\n"
  )
  invisible(NULL)
}

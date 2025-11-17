test_that("parse_file() lê corretamente arquivos em blocos (formato RNAhybrid)", {
  txt <- "
MICA_UTR_123
hsa-miR-125a-5p
interaction energy = -28.4 kcal/mol

MICA_UTR_123
hsa-miR-7160-5p
minimum free energy = -31.55
p-value = 0.0041
"
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  df <- parse_file(tmp)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
  expect_equal(df$miRNA[1], "hsa-miR-125a-5p")
  expect_equal(df$miRNA[2], "hsa-miR-7160-5p")
  expect_equal(df$gene, c("MICA", "MICA"))
  expect_equal(df$utr,  c("123","123"))
  expect_equal(df$interaction_energy[1], -28.4)
  expect_equal(df$minimum_free_energy[2], -31.55)
  expect_equal(df$p_value[2], 0.0041)
})



test_that("parse_file() lê corretamente arquivos CSV com colunas padrão", {
  txt <- "miRNA,target,score
hsa-miR-34a-5p,MICA_UTR_555,-25.5
hsa-miR-21-5p,MICA_UTR_555,-18.2"

  tmp <- tempfile(fileext = ".csv")
  writeLines(txt, tmp)

  df <- parse_file(tmp)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
  expect_equal(df$target, c("MICA_UTR_555","MICA_UTR_555"))
  expect_true("score" %in% names(df))
  expect_equal(df$score, c(-25.5, -18.2))
})



test_that("parse_file() identifica separador automaticamente", {
  txt <- "miRNA;target;score
miR-A;GENE_UTR_1;-10"
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  df <- parse_file(tmp)

  expect_s3_class(df, "data.frame")
  expect_equal(df$miRNA, "miR-A")
  expect_equal(df$target, "GENE_UTR_1")
  expect_equal(df$score, -10)
})



test_that("parse_file() falha adequadamente para arquivos inexistentes", {
  expect_error(parse_file("arquivo_que_nao_existe.txt"))
})



test_that("parse_file() retorna erro quando nem tabular nem blocos são detectados", {
  txt <- "isto não deveria ser reconhecido"
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  expect_error(parse_file(tmp))
})



test_that("parse_file() lida com valores numéricos diferentes por linha", {
  txt <- "
GENE_UTR_9
miR-1
score1 = 10
score2 = -2.5

GENE_UTR_9
miR-2
score1 = 50
score2 = -100.3
"

  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  df <- parse_file(tmp)

  expect_equal(nrow(df), 2)
  expect_equal(df$score1, c(10, 50))
  expect_equal(df$score2, c(-2.5, -100.3))
})



test_that("parse_file() não quebra quando uma chave numérica está ausente em alguns blocos", {
  txt <- "
GENE_UTR_7
miR-1
score1 = 10

GENE_UTR_7
miR-2
score2 = 5
"
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  df <- parse_file(tmp)

  expect_true(all(c("score1", "score2") %in% names(df)))
  expect_true(is.na(df$score1[2]))
  expect_true(is.na(df$score2[1]))
})

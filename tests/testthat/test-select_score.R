test_that("select_score usa automaticamente a única coluna numérica", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20)
  )

  expect_message(
    out <- select_score(df),
    "automaticamente"
  )

  expect_equal(out$Score, df$energy)
})

test_that("select_score requer score_column quando há múltiplas colunas numéricas", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  expect_error(
    select_score(df),
    "múltiplas colunas numéricas"
  )
})

test_that("select_score funciona quando score_column é fornecido", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  out <- select_score(df, score_column = "pvalue")
  expect_equal(out$Score, df$pvalue)
})

test_that("select_score lança erro se score_column não existir", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20)
  )

  expect_error(
    select_score(df, score_column = "nao_existe"),
    "não encontrada"
  )
})

test_that("select_score lança erro se a coluna não é numérica", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c("x", "y") # não numérico
  )

  expect_error(
    select_score(df, score_column = "energy"),
    "não é numérica"
  )
})

test_that("select_score lança erro quando não há colunas numéricas", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    texto = c("x", "y")
  )

  expect_error(
    select_score(df),
    "Nenhuma coluna numérica"
  )
})

test_that("select_score lança erro se df não é data.frame", {
  expect_error(
    select_score(NULL),
    "df deve ser um data.frame"
  )
})

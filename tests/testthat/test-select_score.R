test_that("select_score usa automaticamente a unica coluna numerica", {
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

test_that("select_score requer score_column quando ha multiplas colunas numericas", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  expect_error(
    select_score(df),
    "multiplas colunas numericas"
  )
})

test_that("select_score funciona quando score_column e fornecido", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  out <- select_score(df, score_column = "pvalue")
  expect_equal(out$Score, df$pvalue)
})

test_that("select_score lanca erro se score_column nao existir", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20)
  )

  expect_error(
    select_score(df, score_column = "nao_existe"),
    "nao encontrada"
  )
})

test_that("select_score lanca erro se a coluna nao e numerica", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c("x", "y")
  )

  expect_error(
    select_score(df, score_column = "energy"),
    "nao e numerica"
  )
})

test_that("select_score lanca erro quando nao ha colunas numericas", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    texto = c("x", "y")
  )

  expect_error(
    select_score(df),
    "Nenhuma coluna numerica"
  )
})

test_that("select_score lanca erro se df nao e data.frame", {
  expect_error(
    select_score(NULL),
    "df deve ser um data.frame"
  )
})

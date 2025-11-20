test_that("apply_numeric_filters lança erro se df nao e data.frame", {
  expect_error(apply_numeric_filters(NULL), "data.frame")
})

test_that("apply_numeric_filters lança erro se falta a coluna Score", {
  df <- data.frame(a = 1:3)
  expect_error(apply_numeric_filters(df), "Score")
})

test_that("apply_numeric_filters remove NA por padrao", {
  df <- data.frame(
    miRNA = c("a", "b", "c"),
    Score = c(-10, NA, -30)
  )

  out <- apply_numeric_filters(df)

  expect_equal(nrow(out), 2)
  expect_false(any(is.na(out$Score)))
})

test_that("apply_numeric_filters mantem NA quando remove_na = FALSE", {
  df <- data.frame(
    miRNA = c("a", "b", "c"),
    Score = c(-10, NA, -30)
  )

  out <- apply_numeric_filters(df, remove_na = FALSE)

  expect_equal(nrow(out), 3)
  expect_true(any(is.na(out$Score)))
})

test_that("apply_numeric_filters aplica min_value corretamente", {
  df <- data.frame(
    miRNA = c("a", "b", "c"),
    Score = c(-10, -20, -30)
  )

  out <- apply_numeric_filters(df, min_value = -25)

  expect_equal(nrow(out), 2)
  expect_equal(out$Score, c(-10, -20))
})

test_that("apply_numeric_filters aplica max_value corretamente", {
  df <- data.frame(
    miRNA = c("a", "b", "c"),
    Score = c(-10, -20, -30)
  )

  out <- apply_numeric_filters(df, max_value = -20)

  expect_equal(nrow(out), 2)
  expect_equal(out$Score, c(-20, -30))
})

test_that("apply_numeric_filters aplica min_value e max_value juntos", {
  df <- data.frame(
    miRNA = c("a", "b", "c", "d"),
    Score = c(-5, -15, -25, -35)
  )

  out <- apply_numeric_filters(df, min_value = -30, max_value = -10)

  expect_equal(nrow(out), 2)
  expect_equal(out$Score, c(-15, -25))
})

test_that("apply_numeric_filters mantem as demais colunas intactas", {
  df <- data.frame(
    miRNA = c("a", "b", "c"),
    target = c("T1", "T2", "T3"),
    Score = c(-10, -20, -30)
  )

  out <- apply_numeric_filters(df, min_value = -25)

  expect_equal(names(out), names(df))
  expect_equal(out$target, c("T1", "T2"))
})

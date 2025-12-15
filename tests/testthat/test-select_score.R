test_that("select_score automatically uses the single numeric column", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20)
  )

  out <- select_score(df)
  expect_true("Score" %in% colnames(out))


  expect_equal(out$Score, df$energy)
})

test_that("select_score requires score_column when multiple numeric columns exist", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  expect_error(
    select_score(df),
    "numeric columns"
  )

})

test_that("select_score works when score_column is provided", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20),
    pvalue = c(0.01, 0.05)
  )

  out <- select_score(df, score_column = "pvalue")
  expect_equal(out$Score, df$pvalue)
})

test_that("select_score throws an error if score_column does not exist", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c(-10, -20)
  )

  expect_error(
    select_score(df, score_column = "does_not_exist"),
    "not found"
  )
})

test_that("select_score throws an error if the column is not numeric", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    energy = c("x", "y")
  )

  expect_error(
    select_score(df, score_column = "energy"),
    "not numeric"
  )
})

test_that("select_score throws an error when no numeric columns exist", {
  df <- data.frame(
    miRNA = c("a", "b"),
    target = c("T1", "T2"),
    text = c("x", "y")
  )

  expect_error(
    select_score(df),
    "No numeric column"
  )
})

test_that("select_score throws an error if df is not a data.frame", {
  expect_error(
    select_score(NULL),
    "df must be a data.frame"
  )
})

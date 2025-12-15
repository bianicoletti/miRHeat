test_that("mirHeat_table() errors if input is not a data.frame", {
  expect_error(
    mirHeat_table("not_a_df"),
    "df must be a data.frame"
  )
})
test_that("mirHeat_table() errors if miRNA column is missing", {
  df <- data.frame(
    gene = "MICA",
    utr = "1",
    Score = -15
  )

  expect_error(
    mirHeat_table(df),
    "Required column 'miRNA' not found"
  )
})
test_that("mirHeat_table() errors if Score column is missing", {
  df <- data.frame(
    miRNA = "miR-21",
    gene = "MICA",
    utr = "1"
  )

  expect_error(
    mirHeat_table(df),
    "Required column 'Score' not found"
  )
})
test_that("mirHeat_table() creates target column from gene and utr when missing", {
  df <- data.frame(
    miRNA = "miR-21",
    gene = "MICA",
    utr = "1",
    Score = -12
  )

  out <- mirHeat_table(df)

  expect_true("target" %in% colnames(out))
  expect_equal(out$target, "MICA_UTR_1")
})
test_that("mirHeat_table() extracts gene and utr from target when missing", {
  df <- data.frame(
    miRNA = "miR-21",
    target = "MICA_UTR_2",
    Score = -18
  )

  out <- mirHeat_table(df)

  expect_true(all(c("gene", "utr") %in% colnames(out)))
  expect_equal(out$gene, "MICA")
  expect_equal(out$utr, "2")
})

test_that("mirHeat_table() returns standardized columns in correct order", {
  df <- data.frame(
    miRNA = "miR-21",
    gene = "MICA",
    utr = "3",
    Score = -20
  )

  out <- mirHeat_table(df)

  expect_equal(
    colnames(out),
    c("miRNA", "gene", "utr", "target", "Score")
  )
})
test_that("mirHeat_table() returns a tibble", {
  df <- data.frame(
    miRNA = "miR-21",
    gene = "MICA",
    utr = "1",
    Score = -10
  )

  out <- mirHeat_table(df)

  expect_s3_class(out, "tbl_df")
})

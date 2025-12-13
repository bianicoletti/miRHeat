test_that("prepare_for_heatmap() requires a data.frame", {
  expect_error(prepare_for_heatmap(123), "data.frame")
  expect_error(prepare_for_heatmap("abc"), "data.frame")
})

test_that("prepare_for_heatmap() requires the miRNA column", {
  df <- data.frame(target = "G1", Score = -10)
  expect_error(prepare_for_heatmap(df), "miRNA")
})

test_that("prepare_for_heatmap() requires the Score column", {
  df <- data.frame(miRNA = "a", target = "G1")
  expect_error(prepare_for_heatmap(df), "Score")
})

# ------------------------------------------------------------------
# Existing target
# ------------------------------------------------------------------
test_that("prepare_for_heatmap() keeps target when it already exists", {
  df <- data.frame(
    miRNA = "a",
    target = "G1_UTR_1",
    Score = -10
  )

  out <- prepare_for_heatmap(df)

  expect_equal(out$target, "G1_UTR_1")
  expect_equal(out$gene, "G1")
  expect_equal(out$utr, "1")
})


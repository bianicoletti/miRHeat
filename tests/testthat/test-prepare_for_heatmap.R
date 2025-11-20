test_that("prepare_for_heatmap() exige um data.frame", {
  expect_error(prepare_for_heatmap(123), "data.frame")
  expect_error(prepare_for_heatmap("abc"), "data.frame")
})

test_that("prepare_for_heatmap() exige coluna miRNA", {
  df <- data.frame(target = "G1", Score = -10)
  expect_error(prepare_for_heatmap(df), "miRNA")
})

test_that("prepare_for_heatmap() exige coluna Score", {
  df <- data.frame(miRNA = "a", target = "G1")
  expect_error(prepare_for_heatmap(df), "Score")
})

# ------------------------------------------------------------------
# Target já existente
# ------------------------------------------------------------------
test_that("prepare_for_heatmap() mantem target se ele ja existe", {
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


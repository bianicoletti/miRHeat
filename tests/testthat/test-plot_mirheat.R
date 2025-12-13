test_that("plot_mirheat() validates inputs correctly", {
  # df must be a data.frame
  expect_error(plot_mirheat("a"), "data.frame")

  # missing miRNA
  df1 <- data.frame(utr = "1", Score = -10)
  expect_error(plot_mirheat(df1), "miRNA")

  # missing Score
  df2 <- data.frame(miRNA = "a", utr = "1")
  expect_error(plot_mirheat(df2), "Score")
})

test_that("plot_mirheat() works with minimal input", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  # should emit specific messages
  expect_message(plot_mirheat(df), "Building matrix")
  expect_message(plot_mirheat(df), "clusters")

  # should invisibly return a ComplexHeatmap object
  ht <- plot_mirheat(df)
  expect_true(inherits(ht, "Heatmap"))
})

test_that("plot_mirheat() automatically creates target via prepare_for_heatmap()", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    gene = c("G1", "G2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  expect_true(inherits(plot_mirheat(df), "Heatmap"))
})

test_that("plot_mirheat() generates a PNG file when output_file is provided", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  tmp <- tempfile(fileext = ".png")

  expect_message(plot_mirheat(df, output_file = tmp), "Exporting...")

  # file should exist
  expect_true(file.exists(tmp))

  # size should be > 0
  expect_gt(file.size(tmp), 0)

  unlink(tmp)
})

test_that("plot_mirheat() replaces NAs correctly", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(NA, -30)
  )

  ht <- plot_mirheat(df)

  expect_true(inherits(ht, "Heatmap"))
})

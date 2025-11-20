test_that("plot_mirheat() valida entradas corretamente", {
  # df precisa ser data.frame
  expect_error(plot_mirheat("a"), "data.frame")

  # falta miRNA
  df1 <- data.frame(utr = "1", Score = -10)
  expect_error(plot_mirheat(df1), "miRNA")

  # falta Score
  df2 <- data.frame(miRNA = "a", utr = "1")
  expect_error(plot_mirheat(df2), "Score")
})

test_that("plot_mirheat() funciona com entrada minima", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  # deve emitir mensagens específicas
  expect_message(plot_mirheat(df), "Construindo matriz")
  expect_message(plot_mirheat(df), "clusters")

  # deve retornar invisivelmente um objeto ComplexHeatmap
  ht <- plot_mirheat(df)
  expect_true(inherits(ht, "Heatmap"))
})

test_that("plot_mirheat() cria target automaticamente via prepare_for_heatmap", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    gene = c("G1", "G2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  expect_true(inherits(plot_mirheat(df), "Heatmap"))
})

test_that("plot_mirheat() gera arquivo PNG quando output_file e fornecido", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(-10, -20)
  )

  tmp <- tempfile(fileext = ".png")

  expect_message(plot_mirheat(df, output_file = tmp), "Exportando")

  # arquivo deve existir
  expect_true(file.exists(tmp))

  # tamanho > 0
  expect_gt(file.size(tmp), 0)

  unlink(tmp)
})

test_that("plot_mirheat() substitui NAs corretamente", {
  df <- data.frame(
    miRNA = c("mir1", "mir2"),
    utr = c("1", "2"),
    Score = c(NA, -30)
  )

  ht <- plot_mirheat(df)

  expect_true(inherits(ht, "Heatmap"))
})

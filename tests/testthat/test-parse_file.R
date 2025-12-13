test_that("parse_file() correctly reads block-formatted files", {
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



test_that("parse_file() correctly reads CSV files with standard columns", {
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



test_that("parse_file() automatically detects the field separator", {
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



test_that("parse_file() fails appropriately for non-existent files", {
  expect_error(parse_file("file_that_does_not_exist.txt"))
})



test_that("parse_file() returns an error when neither tabular nor block format is detected", {
  txt <- "this should not be recognized"
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)

  expect_error(parse_file(tmp))
})



test_that("parse_file() handles different numeric values per block", {
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



test_that("parse_file() does not break when a numeric key is missing in some blocks", {
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

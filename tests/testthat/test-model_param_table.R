context("test-model_param_table")

test_that("Test dimensions of model_param_table", {

  expect_equal(dim(model_param_table(reduced = FALSE)), c(1248, 33))

  expect_equal(dim(model_param_table(reduced = TRUE)), c(78, 33))

})

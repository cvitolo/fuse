context("test-fuserouting.ranges")

test_that("Get correct parameter ranges for fuserouting", {

  expect_equal(fuserouting.ranges()$timedelay, c(0.01, 5.00))

})

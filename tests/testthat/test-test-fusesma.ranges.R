context("test-fusesma.ranges")

test_that("Get correct parameter ranges for fusesma", {

  expect_equal(fusesma.ranges()$rferr_add, 0)

  expect_equal(fusesma.ranges()$rferr_mlt, 1)

  expect_equal(fusesma.ranges()$maxwatr_1, c(25, 500))

  expect_equal(fusesma.ranges()$maxwatr_2, c(50, 5000))

  expect_equal(fusesma.ranges()$fracten, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$frchzne, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$fprimqb, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$rtfrac1, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$percrte, c(1e-02, 1e+03))

  expect_equal(fusesma.ranges()$percexp, c(1, 20))

  expect_equal(fusesma.ranges()$sacpmlt, c(1, 250))

  expect_equal(fusesma.ranges()$sacpexp, c(1, 5))

  expect_equal(fusesma.ranges()$percfrac, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$iflwrte, c(1e-02, 1e+03))

  expect_equal(fusesma.ranges()$baserte, c(1e-03, 1e+03))

  expect_equal(fusesma.ranges()$qb_powr, c(1, 10))

  expect_equal(fusesma.ranges()$qb_prms, c(0.001, 0.250))

  expect_equal(fusesma.ranges()$qbrate_2a, c(0.001, 0.250))

  expect_equal(fusesma.ranges()$qbrate_2b, c(0.001, 0.250))

  expect_equal(fusesma.ranges()$sareamax, c(0.05, 0.95))

  expect_equal(fusesma.ranges()$axv_bexp, c(0.001, 3.000))

  expect_equal(fusesma.ranges()$loglamb, c(5, 10))

  expect_equal(fusesma.ranges()$tishape, c(2, 5))

})

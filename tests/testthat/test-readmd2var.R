context("test-readmd2var")

test_that("multiplication works", {
  expect_equal(readmd2var(5),
               structure(c("additive_e", "tension1_1", "unlimfrc_2",
                           "arno_x_vic", "perc_f2sat", "sequential",
                           "intflwnone", "rout_gamma"),
                         .Names = c("rferr", "arch1", "arch2", "qsurf",
                                    "qperc", "esoil", "qintf", "q_tdh")))
})

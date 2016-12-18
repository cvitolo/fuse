context("test-generateParameters")

test_that("Example parameter generation", {

  x <- structure(list(rferr_add = 0,
                      rferr_mlt = 1,
                      maxwatr_1 = 220.076697208005,
                      maxwatr_2 = 3385.10251865608,
                      fracten = 0.571165350806239,
                      frchzne = 0.0638556442721411,
                      fprimqb = 0.0651783785818094,
                      rtfrac1 = 0.80924498260935,
                      percrte = 41.0264674696365,
                      percexp = 16.7094968114303,
                      sacpmlt = 6.90925490667795,
                      sacpexp = 2.40113263144184,
                      percfrac = 0.709296633433387,
                      iflwrte = 194.805087148865,
                      baserte = 94.623821242617,
                      qb_powr = 3.531803924947,
                      qb_prms = 0.214332920810052,
                      qbrate_2a = 0.0443311457108546,
                      qbrate_2b = 0.101732776064363,
                      sareamax = 0.317480958245574,
                      axv_bexp = 2.83883580924427,
                      loglamb = 5.54628534604476,
                      tishape = 2.19522694409714,
                      timedelay = 2.76455549789233),
                 .Names = c("rferr_add",
                            "rferr_mlt",
                            "maxwatr_1",
                            "maxwatr_2",
                            "fracten",
                            "frchzne",
                            "fprimqb",
                            "rtfrac1",
                            "percrte",
                            "percexp",
                            "sacpmlt",
                            "sacpexp",
                            "percfrac",
                            "iflwrte",
                            "baserte",
                            "qb_powr",
                            "qb_prms",
                            "qbrate_2a",
                            "qbrate_2b",
                            "sareamax",
                            "axv_bexp",
                            "loglamb",
                            "tishape",
                            "timedelay"),
                 class = "data.frame", row.names = c(NA, -1L))

  set.seed(123)

  expect_equal(generateParameters(1), x)

})

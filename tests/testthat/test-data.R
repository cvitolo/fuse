context("Data")

test_that("Single run", {

  data(fuse_hydrological_timeseries)
  myDELTIM <- 1
  myMID <- 60

  set.seed(1)
  parameters <- generateParameters(1)

  x <- round(coredata(fuse(DATA = fuse_hydrological_timeseries, mid = myMID, deltim = myDELTIM,
            ParameterSet = parameters)), 3)

  # dput(x, 'fuse/inst/tests/testthat/example01')
  y <- dget(system.file(package = 'fuse', 'inst/tests/testthat/example01'))

  expect_that(all(x==y), equals(TRUE))

})

test_that("Ensemble run", {

  data(fuse_hydrological_timeseries)
  myDELTIM <- 1
  mids <- c(60, 230, 342, 426)

  set.seed(1)
  parameters <- generateParameters(10)
  numberOfRuns <- 10

  discharges <- matrix(NA,ncol=4*numberOfRuns,nrow=dim(fuse_hydrological_timeseries)[1])
  kCounter <- 0

  for (m in 1:4){

    myMID <- mids[m]

    for (pid in 1:numberOfRuns){

      kCounter <- kCounter + 1
      ParameterSet <- as.list(parameters[pid,])

      discharges[,kCounter] <- round(coredata(fuse(fuse_hydrological_timeseries, myMID, myDELTIM, parameters[pid,])), 3)

    }
  }

  # dput(discharges, 'fuse/inst/tests/testthat/example02')
  y <- dget(system.file(package = 'fuse', 'inst/tests/testthat/example02'))

  expect_that(all(discharges==y), equals(TRUE))

})

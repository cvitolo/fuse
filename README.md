
-   [fuse: Framework for Understanding Structural Errors (Hydrological Modelling)](#fuse-framework-for-understanding-structural-errors-hydrological-modelling)
-   [Installation](#installation)
-   [Sample data and parameters](#sample-data-and-parameters)
-   [Applications](#applications)
    -   [Example usage with 1 model structure](#example-usage-with-1-model-structure)
    -   [Ensemble modelling](#ensemble-modelling)
-   [Use fuse with hydromad](#use-fuse-with-hydromad)
-   [Meta](#meta)

<!-- Edit the README.Rmd only!!! The README.md is generated automatically from README.Rmd. -->
fuse: Framework for Understanding Structural Errors (Hydrological Modelling)
----------------------------------------------------------------------------

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.14005.svg)](http://dx.doi.org/10.5281/zenodo.14005) [![status](http://joss.theoj.org/papers/392a55daada04a86f95eaa8da134a28d/status.svg)](http://joss.theoj.org/papers/392a55daada04a86f95eaa8da134a28d)

[![Build Status](https://travis-ci.org/cvitolo/fuse.svg)](https://travis-ci.org/cvitolo/fuse.svg?branch=master) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/cvitolo/fuse?branch=master&svg=true)](https://ci.appveyor.com/project/cvitolo/fuse) [![codecov.io](https://codecov.io/github/cvitolo/fuse/coverage.svg?branch=master)](https://codecov.io/github/cvitolo/fuse?branch=master)

Implementation of the framework for hydrological modelling FUSE described in Clark et al. (2008) and based on the Fortran code provided by M. Clark in 2011. The package consists of two modules: Soil Moisture Accounting module (fusesma.sim) and Gamma routing module (fuserouting.sim). It also contains default parameter ranges (fusesma.ranges and fuserouting.ranges) and three data objects: DATA (sample input dataset), parameters (sample parameters) and modlist (list of FUSE model structures).

Installation
------------

If you have already installed [R](https://cran.r-project.org/) (and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows users), you can install the dependencies for the examples below from CRAN:

``` r
packs <- c("dplyr", "zoo", "tgp", "stats", "utils", "BH", "Rcpp", "testthat",
           "devtools")
new.packages <- packs[!(packs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

and the package itself from Github with [devtools](https://github.com/hadley/devtools):

``` r
devtools::install_github("cvitolo/fuse")
```

Load the package:

``` r
library(fuse)
```

Sample data and parameters
--------------------------

Load sample data (daily time step)

``` r
data(DATA)
myDELTIM <- 1
```

Define parameter ranges

``` r
DefaultRanges <- data.frame(t(data.frame(c(fusesma.ranges(),
                                           fuserouting.ranges()))))
names(DefaultRanges) <- c("Min","Max")
```

Sample parameter set using Latin Hypercube method

``` r
numberOfRuns <- 100
parameters <- tgp::lhs( numberOfRuns, as.matrix(DefaultRanges) )
parameters <- data.frame(parameters)
names(parameters) <- row.names(DefaultRanges)
```

Alternatively, sample parameter set using built-in function

``` r
parameters <- generateParameters(100)
```

Applications
------------

### Example usage with 1 model structure

Define the model to use, e.g. TOPMODEL (MID = 60)

``` r
myMID <- 60
```

Use the built-in function to run FUSE for the 1st sampled parameter set

``` r
x <- fuse(DATA, myMID, myDELTIM, parameters[1,])

plot(x,xlab="",ylab="Streamflow [mm/day]")
```

Run FUSE for all the sampled parameter sets

``` r
plot(DATA$Q,type="l",xlab="",ylab="Streamflow [mm/day]")
allQ <- data.frame(matrix(NA,ncol=numberOfRuns,nrow=dim(DATA)[1]))
for (i in 1:numberOfRuns){
  allQ[,i] <- fuse(DATA, myMID, myDELTIM, parameters[i,])
  lines(zoo::zoo(allQ[,i],order.by=zoo::index(DATA)),col="gray",lwd=0.1)
}
lines(DATA$Q,col="black")
```

### Ensemble modelling

Define a group of model structures to use

``` r
mids <- c(60, 230, 342, 426)
```

Run a multi-model calibration using the Nash-Sutcliffe efficiency as objective function

``` r
indices <- rep(NA,4*numberOfRuns)
discharges <- matrix(NA,ncol=4*numberOfRuns,nrow=dim(DATA)[1])
kCounter <- 0

for (m in 1:4){

  myMID <- mids[m]

  for (pid in 1:numberOfRuns){

    kCounter <- kCounter + 1
    ParameterSet <- as.list(parameters[pid,])
    
    Qrout <- fuse(DATA, myMID, myDELTIM, parameters[pid,])
 
    indices[kCounter] <- qualV::EF(DATA$Q,Qrout)  
    discharges[,kCounter] <- Qrout
    
    }
}
```

Compare results

``` r
bestRun <- which(indices == max(indices))
 
bestModel <- function(runNumber){
 if (runNumber<(numberOfRuns+1)) myBestModel <- "TOPMODEL"
 if (runNumber>(numberOfRuns+1) & runNumber<(2*numberOfRuns+1)) myBestModel <- "ARNOXVIC"
 if (runNumber>(2*numberOfRuns+1) & runNumber<(3*numberOfRuns+1)) myBestModel <- "PRMS"
 if (runNumber>(3*numberOfRuns+1) & runNumber<(4*numberOfRuns+1)) myBestModel <- "SACRAMENTO"
 return(myBestModel)
}
bestModel(bestRun)
 
plot(zoo::coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=0.5)
 
for(pid in 1:(4*numberOfRuns)){
 lines(discharges[,pid], col="gray", lwd=3)
}
 
lines(zoo::coredata(DATA$Q),col="black", lwd=1)
lines(discharges[,bestRun],col="red", lwd=1)
```

How the best simulation of each model structure compare to each other?

``` r
bestRun0060 <- which(indices[1:numberOfRuns] == max(indices[1:numberOfRuns]))
bestRun0230 <- numberOfRuns + which(indices[(numberOfRuns+1):(2*numberOfRuns)] == max(indices[(numberOfRuns+1):(2*numberOfRuns)]))
bestRun0342 <- 2*numberOfRuns + which(indices[(2*numberOfRuns+1):(3*numberOfRuns)] == max(indices[(2*numberOfRuns+1):(3*numberOfRuns)]))
bestRun0426 <- 3*numberOfRuns + which(indices[(3*numberOfRuns+1):(4*numberOfRuns)] == max(indices[(3*numberOfRuns+1):(4*numberOfRuns)]))
 
plot(zoo::coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=1)
lines(discharges[,bestRun0060], col="green", lwd=1)
lines(discharges[,bestRun0230], col="blue", lwd=1)
lines(discharges[,bestRun0342], col="pink", lwd=1)
lines(discharges[,bestRun0426], col="orange", lwd=1)
 
legend("top", 
        c("TOPMODEL", "ARNOXVIC", "PRMS","SACRAMENTO"), 
        col = c("green", "blue", "pink", "orange"),
        lty = c(1, 1, 1, 1))
```

Use fuse with hydromad
----------------------

[Hydromad](http://hydromad.catchment.org/) is an excellent framework for hydrological modelling, optimization, sensitivity analysis and assessment of results. It contains a large set of soil moisture accounting modules and routing functions. Thanks to Joseph Guillaume (hydromad’s maintainer), fuse is now compatible with hydromad and below are some examples Joseph provided to use fuse within the hydromad environment.

``` r
# Install and load hydromad
install.packages(c("latticeExtra", "polynom", "car", 
                   "Hmisc", "reshape", "DEoptim", "coda"))
install.packages("dream", repos="http://hydromad.catchment.org")
install.packages("hydromad", repos="http://hydromad.catchment.org")
library(hydromad) 

# Set the parameter ranges using hydromad.options
hydromad.options(fusesma = fusesma.ranges(),
                 fuserouting = fuserouting.ranges())

# Set up the model using the example dataset (DATA) loaded in the previous examples
modspec <- hydromad(DATA,
                    sma = "fusesma", 
                    routing = "fuserouting", 
                    mid = 1:1248, 
                    deltim = 1)

# Randomly generate 1 parameter set
myNewParameterSet <- parameterSets(coef(modspec, warn=FALSE),1,method="random")

# Run a single simulation using the parameter set generated above
modx <- update(modspec, newpars = myNewParameterSet)

# Generate a summary of the result
summary(modx)

# The instantaneous runoff is
U <- modx$U

# The routed discharge is
Qrout <- modx$fitted.values

# Plot the Observed vs Simulated value
hydromad:::xyplot.hydromad(modx)

# Add the precipitation to the above plot
hydromad:::xyplot.hydromad(modx, with.P=TRUE)

# Calibrate FUSE using hydromad's fitBy method and the Shuffled Complex Evolution algorithm
modfit <- fitBySCE(modspec)

# Get a summary of the result
summary(modfit)
```

Meta
----

-   Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
-   Please [report any issues or bugs](https://github.com/cvitolo/fuse/issues).
-   License: [GPL-3](https://opensource.org/licenses/GPL-3.0)
-   Get citation information for the `fuse` package in R doing `citation(package = 'fuse')`

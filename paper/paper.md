---
title: 'fuse: An R package for ensemble Hydrological Modelling'
bibliography: paper.bib
date: 6 September 2016
tags:
- hydrological modelling
- R
authors:
- affiliation: Imperial College London
  name: Claudia Vitolo
  orcid: 0000-0002-4252-1176
- affiliation: Lutra Consulting
  name: Peter Wells
- affiliation: Lutra Consulting
  name: Martin Dobias
- affiliation: Imperial College London
  name: Wouter Buytaert
  orcid: 0000-0001-6994-4454
---

# Summary

Fuse [@fuse-archive] is an R package [@R-base] that implements the framework for hydrological modelling FUSE [@Clark2008] and based on the Fortran code kindly provided by Martyn Clark in 2011. The package consists of two modules: the soil moisture accounting module (fusesma.sim) and the gamma routing module (fuserouting.sim). 

The fuse framework takes as input rainfall and potential evapotranspiration time series (areal averages over the river catchment area) and returns a simulated time series of river discharges. It can be used to understand the variability of expected hydrological responses based on model structures. 

The package contains default parameter ranges (fusesma.ranges and fuserouting.ranges) and three data objects: fuse_hydrological_timeseries (sample input dataset), parameters (sample parameters) and modlist (list of FUSE model structures).

An early version of this package was used as modelling backend of web applications for the Environmental Virtual Observatory pilot project[@Vitolo2015185,@Wilkinson201538].

The fuse package could in future be submitted to CRAN and included in the Task View dedicated to Analysis of Ecological and Environmental Data (Envirometrics). This already includes a number of packages for hydrological modelling such as [@topmodel], [@dynatopmodel] and [@wasim]. These packages only implement a single model structures,while fuse would complement them providing a framework for ensemble modelling. The package hydromad [@Andrews20111171] also implements ensemble modelling but it is not currently on CRAN.

A new version of the fuse Fortran code was recently released on [GitHub](https://github.com/naddor/fuse). Fortran users are advised to refer to this latest version of fuse. This package is not an interface for the latest Fortran code but any contribution in this direction is welcome.

# References

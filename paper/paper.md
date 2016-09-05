---
title: 'fuse: An R package for ensemble Hydrological Modelling'
bibliography: paper.bib
date: "8 August 2016"
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

The package contains default parameter ranges (fusesma.ranges and fuserouting.ranges) and three data objects: DATA (sample input dataset), parameters (sample parameters) and modlist (list of FUSE model structures).

# References

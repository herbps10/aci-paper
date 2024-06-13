# AdaptiveConformal: An `R` Package for Adaptive Conformal Inference

[![build and publish](https://github.com/computorg/template-computo-R/actions/workflows/build.yml/badge.svg)](https://github.com/computorg/template-computo-R/actions/workflows/build.yml)

Authors: 

- Herbert Susmann -- CEREMADE (UMR 7534), Université Paris-Dauphine PSL, Place du Maréchal de Lattre de Tassigny, Paris, 75016, France
- Antoine Chambaz -- Université Paris Cité, CNRS, MAP5, F-75006 Paris, France
- Julie Josse -- Inria PreMeDICaL team, Idesp, Université de Montpellier

Conformal Inference (CI) is a popular approach for generating finite sample prediction intervals based on the output of any point prediction method when data are exchangeable. Adaptive Conformal Inference (ACI) algorithms extend CI to the case of sequentially observed data, such as time series, and exhibit strong theoretical guarantees without having to assume exchangeability of the observed data. The common thread that unites algorithms in the ACI family is that they adaptively adjust the width of the generated prediction intervals in response to the observed data. We provide a detailed description of five ACI algorithms and their theoretical guarantees, and test their performance in simulation studies. We then present a case study of producing prediction intervals for influenza incidence in the United States based on black-box point forecasts. Implementations of all the algorithms are released as an open-source `R` package, `AdaptiveConformal`, which also includes tools for visualizing and summarizing conformal prediction intervals.

## Build instructions
The main paper can be found in `paper.qmd`. Several helper functions, such as functions that execute the simulation studies, are defined in `helpers.R`. 
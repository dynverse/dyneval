
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://api.travis-ci.org/dynverse/dyneval.svg)](https://travis-ci.org/dynverse/dyneval)
[![codecov](https://codecov.io/gh/dynverse/dynguidelines/branch/master/graph/badge.svg)](https://codecov.io/gh/dynverse/dynguidelines)
<img src="man/figures/logo.png" align="right" />

# Metrics to compare two trajectories

This R package implements several metrics for comparing two single-cell
trajectories.

These include:

  - **Specific metrics**, metrics which look at the similarity of a
    specific part of the trajectory, such as the topology or the
    cellular ordering

| metric\_id     | long\_name                         | category      |
| :------------- | :--------------------------------- | :------------ |
| correlation    | Geodesic distance correlation      | ordering      |
| rf\_nmse       | Random Forest MSE                  | neighbourhood |
| rf\_mse        | Random Forest Normalised MSE       | neighbourhood |
| rf\_rsq        | Random Forest R²                   | neighbourhood |
| lm\_mse        | Linear regression MSE              | neighbourhood |
| lm\_rsq        | Linear regression R²               | neighbourhood |
| edge\_flip     | Edge flip                          | topology      |
| him            | Hamming-Ipsen-Mikhailov similarity | topology      |
| isomorphic     | Isomorphic                         | topology      |
| F1\_branches   | Overlap between the branches       | clustering    |
| F1\_milestones | Overlap between the milestones     | clustering    |

  - **Application metrics**, which assess the accuracy of some
    downstream analyses of trajectories

| metric\_id      | long\_name                     |
| :-------------- | :----------------------------- |
| featureimp\_cor | Feature importance correlation |
| fimp\_ks        | Feature importance enrichment  |

  - **Overall metrics**, which combine several specific and/or
    application metrics to analyse the overall similarity of two
    trajectories

| metric\_id | long\_name    |
| :--------- | :------------ |
| harm\_mean | Harmonic mean |


<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://api.travis-ci.org/dynverse/dyneval.svg)](https://travis-ci.org/dynverse/dyneval) [![codecov](https://codecov.io/gh/dynverse/dynguidelines/branch/master/graph/badge.svg)](https://codecov.io/gh/dynverse/dynguidelines) <img src="man/figures/logo.png" align="right" />

Metrics to compare two trajectories
===================================

This R package implements several metrics for comparing two single-cell trajectories.

These include:

-   **Specific metrics**, metrics which look at the similarity of a specific part of the trajectory, such as the topology or the cellular ordering

| Name                       | Long name                          | Category      |
|:---------------------------|:-----------------------------------|:--------------|
| cor<sub>dist</sub>         | Geodesic distance correlation      | ordering      |
| NMSE<sub>rf</sub>          | Random Forest MSE                  | neighbourhood |
| MSE<sub>rf</sub>           | Random Forest Normalised MSE       | neighbourhood |
| R<sup>2</sup><sub>rf</sub> | Random Forest R²                   | neighbourhood |
| MSE<sub>lm</sub>           | Linear regression MSE              | neighbourhood |
| R<sup>2</sup><sub>lm</sub> | Linear regression R²               | neighbourhood |
| edgeflip                   | Edge flip                          | topology      |
| HIM                        | Hamming-Ipsen-Mikhailov similarity | topology      |
| Isomorphic                 | isomorphic                         | topology      |
| F1<sub>branches</sub>      | Overlap between the branches       | clustering    |
| F1<sub>milestones</sub>    | Overlap between the milestones     | clustering    |

-   **Application metrics**, which assess the accuracy of some downstream analyses of trajectories

| Name                     | Long name                               |
|:-------------------------|:----------------------------------------|
| cor<sub>features</sub>   | Feature importance correlation          |
| wcor<sub>features</sub>  | Feature importance weighted correlation |
| ks<sub>feature</sub>     | Feature importance enrichment ks        |
| wilcox<sub>feature</sub> | Feature importance enrichment wilcox    |

-   **Overall metrics**, which combine several specific and/or application metrics to analyse the overall similarity of two trajectories

| Name                      | Long name       |
|:--------------------------|:----------------|
| mean<sub>harmonic</sub>   | Harmonic mean   |
| mean<sub>geometric</sub>  | Geometric mean  |
| mean<sub>arithmetic</sub> | Arithmetic mean |

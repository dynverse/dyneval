
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://api.travis-ci.org/dynverse/dyneval.svg)](https://travis-ci.org/dynverse/dyneval)
[![codecov](https://codecov.io/gh/dynverse/dyneval/branch/master/graph/badge.svg)](https://codecov.io/gh/dynverse/dyneval)
<img src="man/figures/logo.png" align="right" />

# Metrics to compare two trajectories

This R package implements several metrics for comparing two single-cell
trajectories.

These include:

  - **Specific metrics**, metrics which look at the similarity of a
    specific part of the trajectory, such as the topology or the
    cellular
ordering

| Name                       | Long name                          | Category      |
| :------------------------- | :--------------------------------- | :------------ |
| cor<sub>dist</sub>         | Geodesic distance correlation      | ordering      |
| MSE<sub>rf</sub>           | Random Forest MSE                  | neighbourhood |
| NMSE<sub>rf</sub>          | Random Forest Normalised MSE       | neighbourhood |
| R<sup>2</sup><sub>rf</sub> | Random Forest R²                   | neighbourhood |
| NMSE<sub>lm</sub>          | Linear regression Normalised MSE   | neighbourhood |
| MSE<sub>lm</sub>           | Linear regression MSE              | neighbourhood |
| R<sup>2</sup><sub>lm</sub> | Linear regression R²               | neighbourhood |
| edgeflip                   | Edge flip                          | topology      |
| HIM                        | Hamming-Ipsen-Mikhailov similarity | topology      |
| Isomorphic                 | isomorphic                         | topology      |
| F1<sub>branches</sub>      | Overlap between the branches       | clustering    |
| F1<sub>milestones</sub>    | Overlap between the milestones     | clustering    |

  - **Application metrics**, which assess the accuracy of some
    downstream analyses of trajectories

| Name                     | Long name                               |
| :----------------------- | :-------------------------------------- |
| cor<sub>features</sub>   | Feature importance correlation          |
| wcor<sub>features</sub>  | Feature importance weighted correlation |
| ks<sub>feature</sub>     | Feature importance enrichment ks        |
| wilcox<sub>feature</sub> | Feature importance enrichment wilcox    |

  - **Overall metrics**, which combine several specific and/or
    application metrics to analyse the overall similarity of two
    trajectories

| Name                      | Long name       |
| :------------------------ | :-------------- |
| mean<sub>harmonic</sub>   | Harmonic mean   |
| mean<sub>geometric</sub>  | Geometric mean  |
| mean<sub>arithmetic</sub> | Arithmetic mean |

## Latest changes

Check out `news(package = "dyneval")` or [NEWS.md](inst/NEWS.md) for a
full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Latest changes in dyneval 0.2.0 (30-10-2018)

  - DOCUMENTATION: Add NEWS.md.

  - BUG FIX: Fix long\_name of rf\_mse and rf\_nmse in metrics.tsv.

  - FEATURE: Allow expression\_source to be specified for feature
    importance metrics.

  - BUG FIX: Fix category and type of lm\_nmse in metrics.tsv.

  - MINOR CHANGE: update dyneval for dynfeature 0.2.0.

### Latest changes in dyneval 0.1.0 (10-03-2017)

  - INITIAL RELEASE: dyneval provides an evaluation pipeline and the
    required metrics for evaluating trajectory inference methods.
    Contains the following metrics:
      - Ordering:
          - Geodesic distance correlation: `correlation`
      - Predicting trajectory from expression:
          - Random Forest MSE: `rf_mse`
          - Random Forest Normalised MSE `rf_nmse`
          - Random Forest R²: `rf_rsq`
          - Linear Model MSE: `lm_mse`
          - Linear Model Normalised MSE `lm_nmse`
          - Linear Model R²: `lm_rsq`
      - Topology:
          - Edge flip: `edge_flip`
          - Hamming-Ipsen-Mikhailov similarity: `him`
          - Isomorphic: `isomorphic`
      - Feature importance:
          - Feature importance correlation: `featureimp_cor`
          - Feature importance weighted correlation: `featureimp_wcor`
          - Feature importance enrichment ks: `featureimp_ks`
          - Feature importance enrichment wilcox: `featureimp_wilcox`
      - Clustering:
          - Overlap between the branches: `F1_branches`
          - Overlap between the milestones: `F1_milestones`

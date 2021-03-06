# dyneval 0.9.9 (07-04-2019)

* MINOR CHANGE: Update categories of metrics

* MINOR CHANGE: Move metric functions from dyneval to dynutils.

* MINOR CHANGE: Bump version number for dynwrap 1.0.0

# dyneval 0.2.2 (21-11-2018)

* DOCUMENTATION: Update and expand documentation.

* BUG FIX: Update dyneval for dynwrap 0.3.1 (#50).

* MINOR CHANGE: Clean up imports and suggests.

# dyneval 0.2.1 (14-11-2018)

* BUG FIX: Fix for feature importance scores when only a limited number of features are present

# dyneval 0.2.0 (30-10-2018)

* DOCUMENTATION: Add NEWS.md.

* BUG FIX: Fix long_name of rf_mse and rf_nmse in metrics.tsv.

* FEATURE: Allow expression_source to be specified for feature importance metrics.

* BUG FIX: Fix category and type of lm_nmse in metrics.tsv.

* MINOR CHANGE: update dyneval for dynfeature 0.2.0.

# dyneval 0.1.0 (10-03-2017)

* INITIAL RELEASE: dyneval provides an evaluation pipeline and the required metrics for evaluating trajectory inference methods.
  Contains the following metrics:
   - Ordering:
     * Geodesic distance correlation: `correlation`
   - Predicting trajectory from expression:
     * Random Forest MSE: `rf_mse`
     * Random Forest Normalised MSE `rf_nmse`
     * Random Forest R²: `rf_rsq`
     * Linear Model MSE: `lm_mse`
     * Linear Model Normalised MSE `lm_nmse`
     * Linear Model R²: `lm_rsq`
   - Topology:
     * Edge flip: `edge_flip`
     * Hamming-Ipsen-Mikhailov similarity: `him`
     * Isomorphic: `isomorphic`
   - Feature importance:
     * Feature importance correlation: `featureimp_cor`
     * Feature importance weighted correlation: `featureimp_wcor`
     * Feature importance enrichment ks: `featureimp_ks`
     * Feature importance enrichment wilcox: `featureimp_wilcox`
   - Clustering:
     * Overlap between the branches: `F1_branches`
     * Overlap between the milestones: `F1_milestones`

metrics <- tribble(
  ~metric_id, ~long_name, ~category, ~perfect, ~worst,
  "correlation", "Geodesic distance correlation", "ordering", 1, 0,
  "rf_nmse", "Random Forest MSE", "neighbourhood", 1, 0,
  "rf_rsq", "Random Forest R²", "neighbourhood", 1, 0,
  "lm_nmse", "Linear regression MSE", "neighbourhood", 1, 0,
  "lm_rsq", "Linear regression R²", "neighbourhood", 1, 0,
  "edge_flip", "Edge flip", "topology", 1, 0,
  "featureimp_cor", "Feature importance correlation", "Feature importance", 1, 0
)

devtools::use_data(metrics, overwrite = TRUE)

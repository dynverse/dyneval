library(tibble)

metrics <- tribble(
  ~metric_id, ~name, ~long_name, ~category, ~perfect, ~worst,
  "correlation", "cor[dist]", "Geodesic distance correlation", "ordering", 1, 0,
  "rf_nmse", "NMSE[rf]", "Random Forest MSE", "neighbourhood", 1, 0,
  "rf_rsq", "R[rf]^2", "Random Forest R²", "neighbourhood", 1, 0,
  "lm_nmse", "NMSE[lm]", "Linear regression MSE", "neighbourhood", 1, 0,
  "lm_rsq", "R[lm]^2", "Linear regression R²", "neighbourhood", 1, 0,
  "edge_flip", "edgeflip", "Edge flip", "topology", 1, 0,
  "featureimp_cor", "cor[features]", "Feature importance correlation", "Feature importance", 1, 0
)

devtools::use_data(metrics, overwrite = TRUE)

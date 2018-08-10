library(tibble)

metrics <- tribble(
  ~metric_id, ~plotmath, ~latex, ~long_name, ~category, ~perfect, ~worst,
  "correlation", "cor[dist]", "\\mathit{cor}_{\\textrm{dist}}", "Geodesic distance correlation", "ordering", 1, 0,
  "rf_nmse", "NMSE[rf]", "\\mathit{NMSE}_{rf}", "Random Forest MSE", "neighbourhood", 1, 0,
  "rf_mse", "MSE[rf]", "\\mathit{MSE}_{rf}", "Random Forest Normalised MSE", "neighbourhood", 0, 0.3,
  "rf_rsq", "R[rf]^2", "R^{2}_{rf}", "Random Forest R²", "neighbourhood", 1, 0,
  "lm_nmse", "NMSE[lm]", "\\mathit{NMSE}_{lm}", "Linear regression Normalised MSE", "neighbourhood", 1, 0,
  "lm_mse", "MSE[lm]", "\\mathit{MSE}_{lm}", "Linear regression MSE", "neighbourhood", 0, 0.3,
  "lm_rsq", "R[lm]^2", "R^{2}_{lm}", "Linear regression R²", "neighbourhood", 1, 0,
  "edge_flip", "edgeflip", "\\textrm{edgeflip}", "Edge flip", "topology", 1, 0,
  "featureimp_cor", "cor[features]", "\\mathit{cor}_{\\textrm{features}}", "Feature importance correlation", "Feature importance", 1, 0,
  "F1_branches", "F1[branches]", "\\mathit{F1}_{\\textit{branches}}", "Overlap between the branches", "mapping", 1, 0,
  "F1_milestones", "F1[milestones]", "\\mathit{F1}_{\\textit{milestones}}", "Overlap between the milestones", "mapping", 1, 0,
  "harm_mean", "harmonic mean", "\\textrm{harmonic mean}", "Harmonic mean", "Average", 1, 0
)

devtools::use_data(metrics, overwrite = TRUE)

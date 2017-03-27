imp_trajectory_inference <- list(
  wrap_method(
    method_name = "princurve_princurve",
    method_type = "trajectory_inference",
    method_function = function(space) {
      fit <- princurve::principal.curve(space, trace = F, stretch = 0, plot.true = F)
      path <- fit$s[fit$tag,,drop=F]
      time <- setNames(fit$lambda, rownames(space))
      list(
        pseudotime = time, trajectory = path
      )
    },
    parameter_sets = list(list()),
    required_namespaces = c("princurve")
  ),
  wrap_method(
    method_name = "princurve_SCORPIUS",
    method_type = "trajectory_inference",
    method_function = function(space, num_clusters) {
      fit <- SCORPIUS::infer.trajectory(space, k = num_clusters)
      list(
        pseudotime = fit$time, trajectory = fit$path
      )
    },
    parameter_sets = list(list(num_clusters = 4)),
    required_namespaces = c("SCORPIUS")
  )
)

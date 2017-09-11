#' Description for mfa
#' @export
description_mfa <- function() create_description(
  name = "mfa",
  short_name = "mfa",
  package_required = c("mfa"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "b", lower = 1, upper = 10, default = 2),
    makeIntegerParam(id = "iter", lower = 20, upper = 5000, default = 2000),
    makeIntegerParam(id = "thin", lower = 1, upper = 20, default = 1),
    makeIntegerParam(id = "pc_initialise", lower = 1, upper = 5, default = 1),
    makeNumericParam(id = "prop_collapse", lower = 0, upper = 1, default = 0),
    makeLogicalParam(id = "scale_input",default = TRUE),
    makeLogicalParam(id = "zero_inflation",default = FALSE)
  ),
  properties = c(),
  run_fun = run_mfa,
  plot_fun = plot_mfa
)

run_mfa <- function(
  counts,
  b = 2,
  iter=2000,
  thin=1,
  zero_inflation=FALSE,
  pc_initialise=1,
  prop_collapse=0,
  scale_input=TRUE
) {
  requireNamespace("mfa")

  m <- mfa::mfa(counts, b=b, iter=iter, thin=thin, zero_inflation=zero_inflation, pc_initialise = pc_initialise, prop_collapse=prop_collapse, scale_input = scale_input)
  ms <- summary(m)

  milestone_network <- tibble(from="M0", to=paste0("M", seq_len(b)), length=1)
  progressions <- ms %>% mutate(
    from="M0",
    to=paste0("M", branch),
    percentage=(pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime)),
    cell_id=rownames(counts)
  )

  wrap_ti_prediction(
    ti_type = "multifurcating",
    id = "mfa",
    cell_ids = rownames(counts),
    milestone_ids = unique(c(milestone_network$from, milestone_network$to)),
    milestone_network = milestone_network,
    progressions = progressions %>% dplyr::select(cell_id, from, to, percentage),
    m_fit = m,
    ms = ms
  )
}

#' @import ggplot2
plot_mfa <- plot_default

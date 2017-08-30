wrap <- function(milestone_network, progressions, name="toy") {
  task <- dyneval::wrap_ti_prediction(
    "toy",
    name,
    unique(progressions$cell_id),
    unique(c(milestone_network$from, milestone_network$to)),
    milestone_network,
    progressions=progressions
  )
}

generate_toy_datasets <- function() {
  settings <- expand.grid(ti_type=c("linear", "bifurcating", "cycle"), replicate=1:5) %>% as.data.frame()

  purrrlyr::invoke_rows(function(ti_type, ...) {
    milestone_network = dyngen::generate_toy_milestone_network(ti_type)
    progressions = dyngen::random_progressions_tented(milestone_network)

    wrap(milestone_network, progressions)
  }, settings) %>% pull(.out)
}

tasks <- generate_toy_datasets()
saveRDS(tasks, paste0(tempdir(), "/dyneval_test_datasets.rds"))

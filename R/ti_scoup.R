#' Description for SCOUP
#' @export
description_scoup <- function() create_description(
  name = "SCOUP",
  short_name = "SCOUP",
  package_required = c(),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "nbranch", lower = 1L, upper = 20L, default = 3),
    makeIntegerParam(id = "m", lower = 2L, upper = 1000L, default = 50),
    makeIntegerParam(id = "M", lower = 2L, upper = 10000L, default = 50)
    ## TODO.... :(
  ),
  properties = c(),
  run_fun = run_scoup,
  plot_fun = plot_scoup
)

run_scoup <- function(
  counts,
  nbranch = 2,
  m = 20,
  M = 20,
  ndim=3
  ) {
  tmp_dir <- paste0(tempfile(), "/")
  dir.create(tmp_dir)
  scoup_dir <- "/home/wouters/thesis/tools/SCOUP"

  nbranch <- 2

  start_group <- cell_grouping %>% filter(cell_id == special_cells$start_cell_id) %>% pull(group_id)
  start_ix <- cell_grouping %>% filter(group_id==start_group) %>% pull(cell_id)

  vars <- apply(counts[start_ix,], 2, var)
  means <- apply(counts[start_ix,], 2, mean)
  distr.df <- data.frame(i = seq_along(vars) - 1, means, vars)

  write.table(t(counts), file = paste0(tmp_dir, "1"), sep = "\t", row.names = F, col.names = F)
  write.table(distr.df, file = paste0(tmp_dir, "2"), sep = "\t", row.names = F, col.names = F)

  system(glue::glue("{scoup_dir}/sp {tmp_dir}1 {tmp_dir}2 {tmp_dir}3 {tmp_dir}4 {ncol(counts)} {nrow(counts)} {ndim}"))
  system(glue::glue("{scoup_dir}/scoup -k {nbranch} {tmp_dir}1 {tmp_dir}2 {tmp_dir}3 {tmp_dir}4 {tmp_dir}5 {tmp_dir}6 {ncol(counts)} {nrow(counts)} -m {m} -M {M}"))

  model <- read.table(paste0(tmp_dir, "5"))
  colnames(model) <- c("time", paste0("M", seq_len(ncol(model)-1)+1))
  model <- model %>% mutate(cell_id = rownames(counts))

  maxtime <- max(model$time)
  progressions <- model %>%
    gather("to", "percentage", -cell_id, -time) %>%
    mutate(percentage = percentage * (time/maxtime)) %>%
    mutate(from="M1")

  milestone_network <- tibble(from="M1", to=paste0("M", seq_len(nbranch) + 1), length=1)

  wrap_ti_prediction(
    ti_type = "branching",
    id = "SCOUP",
    cell_ids = rownames(counts),
    milestone_ids = unique(c(milestone_network$from, milestone_network$to)),
    milestone_network = milestone_network %>% select(from, to, length),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    model = model
  )
}


plot_scoup <- function(prediction) {
  model %>% gather("endstate", "percentage", -cell_id, -time) %>%
    ggplot() + geom_point(aes(time, percentage, color=endstate)) + facet_grid(endstate~.)
}

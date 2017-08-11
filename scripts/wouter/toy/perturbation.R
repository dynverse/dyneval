perturb_gs <- function(task) {
  task
}


perturb_switch_n_cells <- function(task, n=length(task$cell_ids)) {
  the_chosen_ones <- sample(task$cell_ids, n)

  mapper <- set_names(task$cell_ids, task$cell_ids)
  mapper[match(the_chosen_ones, mapper)] <- rev(the_chosen_ones)

  task$milestone_percentages$cell_id <- mapper[task$milestone_percentages$cell_id]
  task$progression$cell_id <- mapper[task$progression$cell_id]

  recreate_task(task)
}

perturb_switch_two_cells <- function(task) perturb_switch_n_cells(task, 2)

perturb_switch_all_cells <- function(task) perturb_switch_n_cells(task, length(unique(task$progressions$cell_id)))

recreate_task <- function(task) {
  task <- wrap_ti_prediction(
    task$ti_type,
    task$id,
    task$cell_ids,
    task$milestone_ids,
    task$milestone_network,
    progression=task$progression
  )
}


rename_toy <- function(task, toy_id) {task$id<-toy_id;task}

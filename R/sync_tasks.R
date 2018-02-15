#' Save a tibble of tasks as separate files, locally and on a remote
#'
#' @param tasks A tibble of tasks
#' @param local_tasks_folder An empty directory on the local device
#' @param remote_tasks_folder An empty directory on a remote device
#'
#' @export
#'
#' @importFrom testthat expect_length
sync_tasks <- function(
  tasks,
  local_tasks_folder,
  remote_tasks_folder
) {
  requireNamespace("PRISM")
  config <- PRISM::override_qsub_config()
  remote <- config$remote

  PRISM:::mkdir_remote(local_tasks_folder, remote = "")
  PRISM:::mkdir_remote(remote_tasks_folder, remote = remote)

  testthat::expect_length(PRISM:::ls_remote(local_tasks_folder, remote = ""), 0)
  testthat::expect_length(PRISM:::ls_remote(remote_tasks_folder, remote = remote), 0)

  for (ti in seq_len(nrow(tasks))) {
    task_id <- tasks$id[[ti]]
    folder <- str_replace(task_id, "[^/]*$", "")

    PRISM:::mkdir_remote(paste0(local_tasks_folder, "/", folder), remote = "")

    task <- tasks %>% slice(ti)
    write_rds(task, paste0(local_tasks_folder, "/", task_id, ".rds"))
  }

  PRISM:::rsync_remote(
    remote_src = "",
    path_src = paste0(local_tasks_folder, "/"),
    remote_dest = remote,
    path_dest = paste0(remote_tasks_folder, "/")
  )

  tasks$id
}

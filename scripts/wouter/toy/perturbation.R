## Again gold standard
perturb_gs <- function(task) {
  task
}

# Switch cells
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



# Break cycle
# this should be easier, no?
perturb_break_cycles <- function(task) {
  # task <- generate_cycle()

  net <- task$milestone_network %>% mutate(linkid = seq_len(n()))
  visited <- c()
  remove_links <- c()

  walker <- function(net, curmilestone, curlinkid) {
    if(curmilestone %in% visited) {
      remove_links <<- c(remove_links, curlinkid)
    } else {
      visited <<- c(visited, curmilestone)
      net %>%
        filter(from == curmilestone) %>%
        {walk2(.$to, .$linkid, function(to, linkid) {walker(net, to, linkid)})}
    }
  }

  walker(net, net$from[[1]], NULL)

  new_milestone_n <- 1
  for(remove_link in remove_links) {
    from <- task$milestone_network[remove_link, ]$from
    to <- task$milestone_network[remove_link, ]$to
    length <- task$milestone_network[remove_link, ]$length

    newto <- paste0("NM", new_milestone_n)
    new_milestone_n <- new_milestone_n + 1

    task$progressions[(task$progressions$from == from) & (task$progressions$to == to),]$to = newto
    task$milestone_network <- task$milestone_network %>% add_row(from=from, to=newto, length=length)

    task$milestone_ids <- c(task$milestone_ids, newto)
  }

  task$milestone_network <- task$milestone_network[-remove_links, ]

  recreate_task(task)
}

# Join linear
perturb_join_linear <- function(task) {
  # task <- generate_linear()
  if(nrow(task$milestone_network) != 1) {stop("joining non-linear trajectories not supported")}

  length <- task$milestone_network$length

  task$milestone_network <- tibble::tribble(
    ~from, ~to, ~length,
    "M1", "M2", length/3,
    "M2", "M3", length/3,
    "M3", "M1", length/3
  )

  task$milestone_ids <- c("M1", "M2", "M3")

  task$progressions <- task$progressions %>%
    mutate(
      from = task$milestone_ids[as.numeric(percentage*3) + 1],
      to = task$milestone_ids[((as.numeric(percentage*3) +1) %% 3 + 1)],
      percentage = (percentage * 3) %% 1
    )

  recreate_task(task)
}


# Split linear to bifurcation
perturb_split_linear <- function(task) {
  # task <- generate_linear()
  if(nrow(task$milestone_network) != 1) {stop("joining non-linear trajectories not supported")}

  length <- task$milestone_network$length

  task$milestone_network <- tibble::tribble(
    ~from, ~to, ~length,
    "M1", "M2", length/2,
    "M2", "M3", length/2,
    "M2", "M4", length/2
  )

  task$milestone_ids <- c("M1", "M2", "M3", "M4")

  task$progressions <- task$progressions %>%
    mutate(
      from = ifelse(percentage > 0.5, "M2", "M1"),
      to = ifelse(percentage > 0.5, sample(c("M3", "M4"), n(), replace=TRUE), "M2"),
      percentage = (percentage * 2) %% 1
    )
  recreate_task(task)
}



# Warping the times
# very quick and dirty way to wrap, but it works :p
perturb_warp <- function(task) {
  # task <- generate_linear()

  task$progressions$percentage <- task$progressions$percentage^(2^runif(1, -2, 2))

  recreate_task(task)
}


### Some helper functions-------------------

# Recreate task, forcing a reculaculation of geodesic distances
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

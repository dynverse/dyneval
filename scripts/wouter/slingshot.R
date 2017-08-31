#library(tidyverse)
#library(magrittr)

#.datasets_location = "../dyngen/results/4/"
#dataset <- dyngen::load_dataset(dyngen::overviewer("datasets")$id[[30]])

#' @importFrom magrittr %<>%
run_slingshot <- function(
  counts,
  ndim = 3,
  nclus = 10,
  dimred_name = "pca"
) {

  dimred_func <- match.fun(paste0("dimred_", dimred_name))

  # normalization
  # from the vignette of slingshot
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

  expression <- FQnorm(t(counts))

  # dimensionality reduction
  space <- dimred_func(counts, ndim=ncomp)

  # clustering
  labels <- kmeans(space, centers = nclus)$cluster

  #plot(space, col = rainbow(nclus)[labels], pch=16, asp = 1)

  sds <- slingshot::slingshot(space, labels)

  pt <- slingshot::pseudotime(sds)
  #pt %>% {.[is.na(.)] = 0;.} %>% pheatmap::pheatmap(scale="none", cluster_cols=FALSE)

  #plot(sds)

  # goal is to divide the different curves in bundles, defined by the combination of lineages
  # we do this recursively and build up the network between the lineage combinations
  sds@lineages
  special_clusters <- sds@connectivity %>% apply(1, sum) %>% keep(~. != 2) %>% names

  bundles <- map(sds@lineages, function(x) {
    prevstart <- 1
    bundles <- list()
    for(i in seq_along(x)[-1]) {
      if(x[[i]] %in% special_clusters) {
        bundles <- c(bundles, list(x[prevstart:i]))
        prevstart <- i +1
      }
    }
    bundles
  })

  combinedbundles <- list()
  bundleid <- 1
  for(subbundles_id in seq_along(bundles)) {
    subbundles <- bundles[[subbundles_id]]
    prevbundle <- 0
    for (bundle in subbundles) {
      already_found <- map_lgl(combinedbundles, ~setequal(.$states, bundle))
      if(any(already_found)) {
        prevbundle <- combinedbundles[[which(already_found)[[1]]]]$id

        combinedbundles[[which(already_found)[[1]]]]$lineages %<>% c(names(bundles)[[subbundles_id]])
        print(prevbundle)
      } else {
        combinedbundles <- c(combinedbundles, list(list(states = bundle, id = bundleid, prevbundle = prevbundle, lineages = names(bundles)[[subbundles_id]])))

        prevbundle <- bundleid

        bundleid <- bundleid + 1
      }
    }
  }
  combinedbundles <- combinedbundles %>% list_as_tibble()

  # to which set of lineages does each cell belong
  combinations_cell <- apply(pt, 1, function(x) names(sds@lineages)[which(!is.na(x))])

  # map the known possible combinations of trajectories (each corresponding to a "bundle" of trajectories) to the combinations of the cells
  tos <- combinations_cell %>% map_int(function(combination) {
    as.integer(combinedbundles$id[map_lgl(combinedbundles$lineages, setequal, combination)][[1]])
  })

  # now create the milestone net and progressions
  milestone_network <- combinedbundles %>%
    rename(from=prevbundle, to=id) %>%
    select(from, to) %>%
    mutate(from=paste0("M", from), to=paste0("M", to))

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  progressions <- tibble(
    time = pt %>% apply(1, min, na.rm=TRUE),
    to = paste0("M", tos),
    cell_id = rownames(pt)
  ) %>%
    group_by(to) %>%
    mutate(percentage = (time - min(time)) / (max(time) - min(time))) %>%
    left_join(milestone_network, by="to") %>%
    ungroup()

  # add length using progressions
  milestone_network <- left_join(
    milestone_network,
    progressions %>% group_by(from, to) %>% summarise(length=max(time)-min(time)),
    by=c("from", "to")
  )

  task <-  wrap_ti_prediction(
    ti_type = "slingshot",
    id = "slingshot",
    cell_ids = colnames(expression),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network %>% select(from, to, length),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    dimred_samples = space,
    dimred_clust = labels
  )
}

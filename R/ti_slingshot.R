#' Description for slingshot
#' @export
description_slingshot <- function() create_description(
  name = "slingshot",
  short_name = "slngsht",
  package_loaded = c("slingshot"),
  package_required = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", lower = 2L, upper = 20L, default = 3L),
    makeIntegerParam(id = "nclus", lower = 2L, upper = 40L, default = 5L),
    makeNumericParam(id = "shrink", lower = 0, upper = 1, default=1),
    makeLogicalParam(id = "reweight", default=TRUE),
    makeLogicalParam(id = "drop.multi", default=TRUE),
    makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "loess", "periodic.lowess")),
    makeDiscreteParam(id = "shrink.method", default = "cosine", values = c("cosine", "tricube", "density")),
    makeDiscreteParam(id = "dimred_name", values = names(list_dimred_methods()), default="pca")

  ),
  properties = c(),
  run_fun = run_slingshot,
  plot_fun = plot_slingshot
)

run_slingshot <- function(
  counts,
  start_cell_id = NULL,
  end_cell_ids = NULL,
  ndim = 3,
  nclus = 5,
  dimred_name = "pca",
  shrink=1,
  reweight=TRUE,
  drop.multi=TRUE,
  thresh=0.001,
  maxit=15,
  stretch=2,
  smoother = "smooth.spline",
  shrink.method = "cosine"
) {
  requireNamespace("slingshot")

  # normalization & preprocessing --------------------
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
  space <- dimred(counts, method = dimred_name, ndim = ndim)

  # clustering
  labels <- kmeans(space, centers = nclus)$cluster

  # actual slingshot algorithm ----------------
  if(!is.null(start_cell_id)) {
    start.clus <- labels[[start_cell_id]]
  } else {
    start.clus <- NULL
  }
  if(!is.null(end_cell_ids)) {
    end.clus <- unique(labels[end_cell_ids])
  } else {
    end.clus <- NULL
  }


  sds <- slingshot::slingshot(
    space,
    labels,
    start.clus=start.clus,
    end.clus=end.clus,
    shrink=shrink,
    reweight=reweight,
    drop.multi=drop.multi,
    thresh=thresh,
    maxit=maxit,
    stretch=stretch,
    smoother = smoother,
    shrink.method = shrink.method
  )
  pt <- slingshot::pseudotime(sds)

  #pt %>% {.[is.na(.)] = 0;.} %>% pheatmap::pheatmap(scale="none", cluster_cols=FALSE)
  #plot(sds)

  # postprocessing ---------------------------

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

  # to which set of lineages does each cell belong, this does not work because sometimes strange combinations are chosen here outside of the known bundle
  combinations_cell <- apply(pt, 1, function(x) names(sds@lineages)[which(!is.na(x))])
  # we therefore work at state (clustering) level
  states_cell <- sds@clusterLabels

  # map the known possible combinations of trajectories (each corresponding to a "bundle" of trajectories) to the combinations of the cells
  tos <- states_cell %>% map_int(function(state) {
    #as.integer(combinedbundles$id[map_lgl(combinedbundles$lineages, setequal, combination)][[1]])
    as.integer(combinedbundles$id[map_lgl(combinedbundles$states, ~state %in% .)][[1]])
  })

  # now create the milestone net and progressions
  milestone_network <- combinedbundles %>%
    rename(from=prevbundle, to=id) %>%
    select(from, to) %>%
    mutate(from=paste0("M", from), to=paste0("M", to), directed=TRUE)

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  progressions <- tibble(
    time = pt %>% apply(1, min, na.rm=TRUE),
    to = paste0("M", tos),
    cell_id = rownames(pt)
  ) %>%
    group_by(to) %>%
    mutate(percentage = dynutils::scale_minmax(time)) %>%
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
    milestone_network = milestone_network %>% select(from, to, length, directed),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    dimred_samples = space,
    dimred_clust = labels,
    sds = sds
  )
}

#' @importFrom graphics pairs
plot_slingshot <- function(ti_predictions) {
  graphics::pairs(ti_predictions$sds, horInd=2, verInd=1)
}


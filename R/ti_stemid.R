#' @import ParamHelpers
#' @export
description_stemid <- function() {
  list(
    name = "StemID",
    short_name = "StemID",
    package = c("tidyverse"),
    par_set = makeParamSet(
      makeIntegerParam(id = "clustnr", lower = 20L, default = 30L, upper = 100L),
      makeIntegerParam(id = "bootnr", lower = 20L, default = 50, upper = 100L),
      makeDiscreteParam(id = "metric", default = "pearson", values = c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
      makeDiscreteParam(id = "num_cluster_method", default = "sat", values = c("sat", "gap", "manual")),
      makeDiscreteParam(id = "SE.method", default = "Tibs2001SEmax", values = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")),
      makeNumericParam(id = "SE.factor", default = .25, lower = 0, upper = 1),
      makeIntegerParam(id = "B.gap", lower = 20L, default = 50, upper = 100L),
      makeIntegerParam(id = "cln", lower = 20L, default = 30, upper = 100L),
      makeDiscreteParam(id = "FUNcluster", default = "kmedoids", values = c("kmedoids", "kmeans", "hclust")),
      makeDiscreteParam(id = "dimred_method", default = "tsne", values = c("tsne", "sammon", "tsne_initcmd")),
      makeIntegerParam(id = "outminc", lower = 0, default = 5, upper = 100L),
      makeIntegerParam(id = "outlg", lower = 0, default = 2, upper = 100L),
      makeNumericParam(id = "probthr", lower = -10, default = -3, upper = -1, trafo = function(x) 10^x),
      makeNumericParam(id = "thr_lower", lower = -100, default = -40, upper = -1),
      makeNumericParam(id = "thr_upper", lower = -100, default = -40, upper = -1),
      makeNumericParam(id = "outdistquant", lower = 0, default = .95, upper = 1),
      makeLogicalParam(id = "nmode", default = F),
      makeNumericParam(id = "pdishuf", lower = 2, default = log10(2000), upper = 3.5, trafo = function(x) ceiling(10)^x),
      makeNumericParam(id = "pthr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
      makeNumericParam(id = "pethr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x)
    ),
    properties = c(),
    run_fun = run_stemid,
    plot_fun = plot_stemid,
    forbidden = quote(thr_lower > thr_lower)
  )
}

#' @export
run_stemid <- function(counts,
                       clustnr = 30,
                       bootnr = 50,
                       metric = "pearson",
                       num_cluster_method = "sat",
                       SE.method = "Tibs2001SEmax",
                       SE.factor = .25,
                       B.gap = 50,
                       cln = 0,
                       FUNcluster = "kmedoids",
                       dimred_method = "tsne", # tsne, sammon, tsne_initcmd
                       outminc = 5,
                       outlg = 2,
                       probthr = 1e-3,
                       thr_lower = -40,
                       thr_upper = -1,
                       outdistquant = .95,
                       nmode = F,
                       pdishuf = 2000,
                       pthr = .01,
                       pethr = .01
                       ) {
  ## load class definition and functions
  code_path <- paste0(path.package("dyneval"), "/extra_code/StemID")
  source(paste0(code_path, "/RaceID2_StemID_class.R"))

  # initialize SCseq object with transcript counts
  sc <- SCseq(data.frame(t(counts), check.names = F, stringsAsFactors = F))

  # filtering of expression data
  sc <- filterdata(sc, mintotal = 1, minexpr = 0, minnumber = 0, maxexpr = Inf, downsample = TRUE, dsn = 1)

  # k-medoids clustering
  do_gap <- num_cluster_method == "gap"
  do_sat <- num_cluster_method == "sat"
  sc <- clustexp(sc, clustnr = clustnr, bootnr = bootnr, metric = metric, do.gap = do_gap,
                 sat = do_sat, SE.method = SE.method, SE.factor = SE.factor,
                 B.gap = B.gap, cln = cln, FUNcluster = FUNcluster)

  # compute t-SNE map
  sammonmap <- dimred_method == "sammon"
  initial_cmd <- dimred_method == "tsne_initcmd"
  sc <- comptsne(sc, sammonmap = sammonmap, initial_cmd = initial_cmd)

  # detect outliers and redefine clusters
  thr <- 2^(thr_lower:thr_upper)
  sc <- findoutliers(sc, outminc = outminc, outlg = outlg, probthr = probthr, thr = thr, outdistquant = outdistquant)

  # initialization
  ltr <- Ltree(sc)

  # computation of the entropy
  ltr <- compentropy(ltr)

  # computation of the projections for all cells
  ltr <- projcells(ltr, cthr = 0, nmode = nmode)

  # computation of the projections for all cells after randomization
  ltr <- projback(ltr, pdishuf = pdishuf, nmode = nmode)

  # assembly of the lineage tree
  ltr <- lineagetree(ltr, pthr = pthr, nmode = nmode)

  # get names of the states and the samples
  state_names <- paste0("state_", ltr@ldata$m)
  ids <- rownames(counts)

  # can't find any cluster distances, so calculating manually
  medoids <- compmedoids(ltr@sc@fdata, ltr@sc@cpart)
  cent <- ltr@sc@fdata[,medoids]
  dc <- as.data.frame(1 - cor(cent))
  dimnames(dc) <- list(state_names, state_names)
  trl <- spantree(dc[ltr@ldata$m,ltr@ldata$m])

  gr_df <- data_frame(from = state_names[seq_along(trl$kid)+1], to = state_names[trl$kid], weight = dc[cbind(from, to)])
  gr <- igraph::graph_from_data_frame(gr, directed = F, vertices = state_names)
  state_network <- igraph::distances(gr) %>% reshape2::melt(varnames = c("from", "to"), value.name = "length") %>% filter(from != to)

  # calculating the cell-to-cluster distances manually
  trproj_res <- ltr@trproj$res %>% as_data_frame() %>% rownames_to_column("id") %>% dplyr::select(id, closest = o, furthest = l)
  state_percentages <- bind_rows(lapply(seq_len(nrow(trproj_res)), function(i) {
    id <- trproj_res$id[[i]]
    closest <- trproj_res$closest[[i]]
    furthest <- trproj_res$furthest[[i]]

    if (!is.na(furthest)) {
      dist_closest <- 1 - cor(cent[,closest], ltr@sc@fdata[,id])
      dist_furthest <- 1 - cor(cent[,furthest], ltr@sc@fdata[,id])

      data_frame(id = id, state = state_names[c(closest, furthest)], percentage = 1 - c(dist_closest, dist_furthest) / (dist_closest + dist_furthest))
    } else {
      data_frame(id = id, state = state_names[[closest]], percentage = 1)
    }
  }))

  wrap_ti_prediction(
    ti_type = "linear",
    name = "StemID",
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    dimred_samples = ltr@ltcoord,
    ltr = ltr
  )
}

#' @export
plot_stemid <- function(ti_predictions) {
  # recreating plotmapprojections(ti_predictions$ltr)
  object <- ti_predictions$ltr
  cent <- object@sc@fdata[,compmedoids(object@sc@fdata,object@sc@cpart)]
  dc <- as.data.frame(1 - cor(cent))
  names(dc) <- sort(unique(object@sc@cpart))
  rownames(dc) <- sort(unique(object@sc@cpart))
  trl <- spantree(dc[object@ldata$m,object@ldata$m])

  u <- object@ltcoord[,1]
  v <- object@ltcoord[,2]
  cnl <- object@ldata$cnl
  lp <- object@ldata$lp

  lines <- data.frame(from = cnl[seq_along(trl$kid)+1,], to = cnl[trl$kid,])

  ggplot() +
    geom_point(aes(u, v), col = "gray", size = 2) +
    geom_text(aes(u, v, label = lp, colour = factor(lp))) +
    geom_point(aes(V1, V2), cnl, shape = 1, size = 2) +
    geom_segment(aes(x = from.V1, xend = to.V1, y = from.V2, yend = to.V2), lines) +
    geom_text(aes(V1, V2, label = seq_len(nrow(cnl))), cnl, size = 10) +
    labs(x = "Dim 1", y = "Dim 2") +
    scale_colour_manual(values = object@sc@fcol) +
    theme(legend.position = "none")
}



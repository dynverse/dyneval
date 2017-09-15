#' Description for StemID
#' @export
description_stemid <- function() create_description(
  name = "StemID",
  short_name = "StemID",
  package_loaded = c(),
  package_required = c("StemID", "igraph", "reshape2"),
  par_set = makeParamSet(
    makeIntegerParam(id = "clustnr", lower = 20L, default = 30L, upper = 100L),
    makeIntegerParam(id = "bootnr", lower = 20L, default = 50L, upper = 100L),
    makeDiscreteParam(id = "metric", default = "pearson", values = c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
    makeDiscreteParam(id = "num_cluster_method", default = "sat", values = c("sat", "gap", "manual")),
    makeDiscreteParam(id = "SE.method", default = "Tibs2001SEmax", values = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")),
    makeNumericParam(id = "SE.factor", default = .25, lower = 0, upper = 1),
    makeIntegerParam(id = "B.gap", lower = 20L, default = 50L, upper = 100L),
    makeIntegerParam(id = "cln", lower = 20L, default = 30L, upper = 100L),
    makeDiscreteParam(id = "FUNcluster", default = "kmedoids", values = c("kmedoids", "kmeans", "hclust")),
    makeDiscreteParam(id = "dimred_method", default = "tsne", values = c("tsne", "sammon", "tsne_initcmd")),
    makeIntegerParam(id = "outminc", lower = 0L, default = 5L, upper = 100L),
    makeIntegerParam(id = "outlg", lower = 0L, default = 2L, upper = 100L),
    makeNumericParam(id = "probthr", lower = -10, default = -3, upper = -1, trafo = function(x) 10^x),
    makeNumericParam(id = "thr_lower", lower = -100, default = -40, upper = -1),
    makeNumericParam(id = "thr_upper", lower = -100, default = -40, upper = -1),
    makeNumericParam(id = "outdistquant", lower = 0, default = .95, upper = 1),
    makeLogicalParam(id = "nmode", default = FALSE),
    makeNumericParam(id = "pdishuf", lower = 2, default = log10(2000), upper = 3.5, trafo = function(x) ceiling(10)^x),
    makeNumericParam(id = "pthr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "pethr", lower = -4, default = -2, upper = 0, trafo = function(x) 10^x),
    forbidden = quote(thr_lower > thr_upper)
  ),
  properties = c(),
  run_fun = run_stemid,
  plot_fun = plot_stemid
)

run_stemid <- function(
  counts,
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
  nmode = FALSE,
  pdishuf = 2000,
  pthr = .01,
  pethr = .01
) {
  requireNamespace("StemID")
  requireNamespace("igraph")
  requireNamespace("reshape2")

  # initialize SCseq object with transcript counts
  sc <- StemID:::SCseq(data.frame(t(counts), check.names = FALSE, stringsAsFactors = FALSE))

  # filtering of expression data
  sc <- StemID:::filterdata(sc, mintotal = 1, minexpr = 0, minnumber = 0, maxexpr = Inf, downsample = TRUE, dsn = 1)

  # k-medoids clustering
  do_gap <- num_cluster_method == "gap"
  do_sat <- num_cluster_method == "sat"
  sc <- StemID:::clustexp(sc, clustnr = clustnr, bootnr = bootnr, metric = metric, do.gap = do_gap,
                 sat = do_sat, SE.method = SE.method, SE.factor = SE.factor,
                 B.gap = B.gap, cln = cln, FUNcluster = FUNcluster)

  # compute t-SNE map
  sammonmap <- dimred_method == "sammon"
  initial_cmd <- dimred_method == "tsne_initcmd"
  sc <- StemID:::comptsne(sc, sammonmap = sammonmap, initial_cmd = initial_cmd)

  # detect outliers and redefine clusters
  thr <- 2^(thr_lower:thr_upper)
  sc <- StemID:::findoutliers(sc, outminc = outminc, outlg = outlg, probthr = probthr, thr = thr, outdistquant = outdistquant)

  # initialization
  ltr <- StemID:::Ltree(sc)

  # computation of the entropy
  ltr <- StemID:::compentropy(ltr)

  # computation of the projections for all cells
  ltr <- StemID:::projcells(ltr, cthr = 0, nmode = nmode)

  # computation of the projections for all cells after randomization
  ltr <- StemID:::projback(ltr, pdishuf = pdishuf, nmode = nmode)

  # assembly of the lineage tree
  ltr <- StemID:::lineagetree(ltr, pthr = pthr, nmode = nmode)

  # compute a spanning tree
  ltr <- StemID:::compspantree(ltr)

  dc <- ltr@trl$dc
  milestone_ids <- rownames(dc)
  trl <- ltr@trl$trl
  cent <- ltr@trl$cent

  # converting to milestone network
  gr_df <- data_frame(from = seq_along(trl$kid)+1, to = trl$kid, weight = dc[cbind(from, to)])
  gr <- igraph::graph_from_data_frame(gr_df, directed = FALSE, vertices = milestone_ids)
  milestone_network <- igraph::distances(gr) %>%
    {.[upper.tri(., diag=TRUE)] = NA; .} %>% # no bidirectionality
    reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
    filter(!is.na(length)) %>%
    filter(from != to) %>%
    mutate(from = as.character(from), to = as.character(to), directed=FALSE)

  # calculating the cell-to-cluster distances manually
  trproj_res <- ltr@trproj$res %>% as_data_frame() %>% rownames_to_column("id") %>% select(id, closest = o, furthest = l)
  milestone_percentages <- bind_rows(lapply(seq_len(nrow(trproj_res)), function(i) {
    id <- trproj_res$id[[i]]
    closest <- trproj_res$closest[[i]]
    furthest <- trproj_res$furthest[[i]]

    if (!is.na(furthest)) {
      dist_closest <- 1 - cor(cent[,closest], ltr@sc@fdata[,id])
      dist_furthest <- 1 - cor(cent[,furthest], ltr@sc@fdata[,id])

      data_frame(cell_id = id, milestone_id = milestone_ids[c(closest, furthest)], percentage = 1 - c(dist_closest, dist_furthest) / (dist_closest + dist_furthest))
    } else {
      data_frame(cell_id = id, milestone_id = milestone_ids[[closest]], percentage = 1)
    }
  }))

  milestone_network$from <- paste0("Milestone", milestone_network$from)
  milestone_network$to <- paste0("Milestone", milestone_network$to)
  milestone_ids <- paste0("Milestone", milestone_ids)
  milestone_percentages$milestone_id <- paste0("Milestone", milestone_percentages$milestone_id)

  wrap_ti_prediction(
    ti_type = "linear",
    id = "StemID",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = ltr@ltcoord,
    objects = list(
      trl = ltr@trl,
      ltcoord = ltr@ltcoord,
      cnl = ltr@ldata$cnl,
      lp = ltr@ldata$lp,
      fcol = ltr@sc@fcol
    )
  )
}

plot_stemid <- function(ti_predictions) {
  trl <- ti_predictions$objects$trl$trl
  ltcoord <- ti_predictions$objects$ltcoord
  cnl <- ti_predictions$objects$cnl
  lp <- ti_predictions$objects$lp
  fcol <- ti_predictions$objects$fcol

  u <- ltcoord[,1]
  v <- ltcoord[,2]

  lines <- data.frame(from = cnl[seq_along(trl$kid)+1,], to = cnl[trl$kid,])

  ggplot() +
    geom_point(aes(u, v), col = "gray", size = 2) +
    geom_text(aes(u, v, label = lp, colour = factor(lp))) +
    geom_point(aes(V1, V2), cnl, shape = 1, size = 2) +
    geom_segment(aes(x = from.V1, xend = to.V1, y = from.V2, yend = to.V2), lines) +
    geom_text(aes(V1, V2, label = seq_len(nrow(cnl))), cnl, size = 10) +
    labs(x = "Dim 1", y = "Dim 2") +
    scale_colour_manual(values = fcol) +
    theme(legend.position = "none")
}



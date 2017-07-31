#' @import ParamHelpers
#' @export
description_stemid <- function() {
  list(
    name = "StemID",
    short_name = "StemID",
    package = c("R.matlab", "readr", "tidyverse"),
    par_set = makeParamSet(),
    properties = c(),
    run_fun = run_stemid,
    plot_fun = plot_stemid
  )
}

#' @export
run_stemid <- function(counts) {
  ## load class definition and functions
  code_path <- paste0(path.package("dyneval"), "/extra_code/StemID")
  source(paste0(code_path, "/RaceID2_StemID_class.R"))

  # initialize SCseq object with transcript counts
  sc <- SCseq(data.frame(t(counts), check.names = F, stringsAsFactors = F))
  # filtering of expression data
  sc <- filterdata(sc, mintotal=1, minexpr=0, minnumber=0, maxexpr=Inf, downsample=TRUE, dsn=1)
  # k-medoids clustering
  sc <- clustexp(sc,clustnr=30,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,FUNcluster="kmedoids")
  # compute t-SNE map
  sc <- comptsne(sc)
  # detect outliers and redefine clusters
  sc <- findoutliers(sc, outminc=0, outlg=2, probthr=1e-3, thr=2**-(1:40), outdistquant=.95)

  # initialization
  ltr <- Ltree(sc)
  # computation of the entropy
  ltr <- compentropy(ltr)
  # computation of the projections for all cells
  ltr <- projcells(ltr,cthr=2,nmode=FALSE)
  # computation of the projections for all cells after randomization
  ltr <- projback(ltr,pdishuf=2000,nmode=FALSE)
  # assembly of the lineage tree
  ltr <- lineagetree(ltr,pthr=0.01,nmode=FALSE)
  # determination of significant differentiation trajectories
  ltr <- comppvalue(ltr,pethr=0.01,nmode=FALSE)


  # assuming column o is the cluster of origin, l is the cluster of destination,
  # and h is weighted between 0 and 1 according to how far along it is in its transition from o to l,
  # and h<0 or h>1 means it is an outlier with respect to the transition
  # (h<0 means )
  state_percentages <- ltr@trproj$res %>%
    rownames_to_column("id") %>%
    na.omit %>%
    gather(type, st, o, l) %>%
    mutate(
      state = paste0("state_", st),
      h2 = ifelse(type == "o", 1-h, h),
      percentage = ifelse(h2 > 1, 1, ifelse(h2 < 0, 0, h2))
    ) %>%
    dplyr::select(id, state, percentage)

  state_names <- paste0("state_", ltr@ldata$m)

  ## diagnostic plots
  # histogram of ratio between cell-to-cell distances in the embedded and the input space
  plotdistanceratio(ltr)
  # t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
  plotmap(ltr)
  # visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
  plotmapprojections(ltr)
  # lineage tree showing the projections of all cells in t-SNE space
  plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=0)
  # lineage tree without showing the projections of all cells
  plottree(ltr,showCells=FALSE,nmode=FALSE,scthr=0)
  # heatmap of the enrichment p-values for all inter-cluster links
  plotlinkpv(ltr)
  # heatmap of the link score for all inter-cluster links
  plotlinkscore(ltr)
  # heatmap showing the fold enrichment (or depletion) for significantly enriched or depleted links
  projenrichment(ltr)

  ## extract projections onto all links for all cells in a given cluster i
  x <- getproj(ltr,i=1)
  # heatmap of all projections for cluster i
  pheatmap(x$pr)
  # heatmap of z-score for all projections for cluster i
  pheatmap(x$prz)

  ## extracting all cells on two branches sharing the same cluster and computing differentially expressed genes between these two branches
  x <- branchcells(ltr,list("1.3","1.2"))
  # z-scores for differentially expressed genes
  head(x$diffgenes$z)
  # plotting the cells on the two branches as additional clusters in the t-SNE map
  plottsne(x$scl)


  ## computing the StemID score
  x <- compscore(ltr,nn=1)
  #plotting the StemID score
  plotscore(ltr,1)




  # create final output
  ids <- rownames(counts)
  state_names <- tree$cluster_name
  state_network <- tree %>%
    filter(parent != 0) %>%
    select(to = cluster_name, from_ix = parent) %>%
    left_join(tree %>% select(from = cluster_name, from_ix = cluster), by = "from_ix") %>%
    mutate(length = 1) %>%
    select(from, to, length)
  state_percentages <- data_frame(id = projection_data$Time, state = paste0("Cluster", treemat_data$s[1,]), percentage = 1)

  # remove temporary output
  unlink(dataset_path, recursive = T)

  wrap_ti_prediction(
    ti_type = "linear",
    name = "StemID",
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages
    # ... add more output
  )
}

#' @export
plot_stemid <- function(ti_predictions) {
  #qplot(percent_rank(ti_predictions$state_percentages[,1]), ti_predictions$state_percentages[,1], colour = data$sample_info$group.name)
}


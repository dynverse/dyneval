list_dimred_methods <- function() {
  list(
    pca = dimred_pca,
    mds = dimred_mds,
    tsne = dimred_tsne,
    # dp = dimred_dp,
    ica = dimred_ica,
    lle = dimred_lle
  )
}

dimred <- function(x, method, ...) {
  methods <- list_dimred_methods()
  if (method %in% names(methods)) {
    meth <- methods[[method]]
    params <- list(x = x, ...)
    do.call(meth, params)
  } else {
    stop("Method ", sQuote(method), " not found.")
  }
}

process_dimred <- function(space, rn) {
  space <- as.matrix(space)
  dimnames(space) <- list(rn, paste0("Comp", 1:ncol(space)))
  space
}


dimred_pca <- function(x, ndim = 3) {
  space <- prcomp(t(x))$rotation[,seq_len(ndim)]
  process_dimred(space, rownames(x))
}

# dimred_simlr <- function(x, ndim = 3, nclusters = 4) {
#   requireNamespace("SIMLR")
#   requireNamespace("tsn")
#   result <- SIMLR::SIMLR(t(x), nclusters)
#   S <- result$S
#   space <- tsne::tsne(as.dist(max(S)-S), k = ndim)
#   process_dimred(space, rownames(x))
# }

dimred_mds <- function(x, ndim = 3) {
  requireNamespace("SCORPIUS")
  space <- SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(x), ndim = ndim)
  process_dimred(space, rownames(x))
}

# dimred_mds_sammon <- function(x, ndim = 3) {
#   requireNamespace("SCORPIUS")
#   requireNamespace("MASS")
#   dist <- SCORPIUS::correlation_distance(x)
#   space <- MASS::sammon(dist, k = ndim)$points
#   process_dimred(space, rownames(x))
# }
#
# dimred_mds_isomds <- function(x, ndim = 3) {
#   requireNamespace("SCORPIUS")
#   requireNamespace("MASS")
#   dist <- SCORPIUS::correlation_distance(x)
#   space <- MASS::isoMDS(dist, k = ndim)$points
#   process_dimred(space, rownames(x))
# }
#
# dimred_lmds <- function(x, ndim = 3) {
#   mds.out <- dambiutils::mds_withlandmarks(x %>% as.data.frame, SCORPIUS::correlation_distance, k = ndim, landmark.method = "naive", num.landmarks = min(1000, round(nrow(x)*0.1)), num.seed.landmarks = 10, pca.normalisation = F)
#   process_dimred(mds.out$S, rownames(x))
# }
#
# dimred_mds_smacof <- function(x, ndim = 3) {
#   requireNamespace("SCORPIUS")
#   requireNamespace("smacof")
#   dist <- SCORPIUS::correlation_distance(x)
#   space <- smacof::mds(as.dist(dist), type = "ratio", ndim = ndim)$conf
#   process_dimred(space, rownames(x))
# }

dimred_tsne <- function(x, ndim = 3) {
  requireNamespace("SCORPIUS")
  requireNamespace("Rtsne")
  requireNamespace("stats")
  space <- Rtsne::Rtsne(stats::as.dist(SCORPIUS::correlation_distance(x)), dims = ndim, is_distance = TRUE)$Y
  process_dimred(space, rownames(x))
}

# dimred_dp <- function(x, ndim = 3, neigen = 3) {
#   requireNamespace("SCORPIUS")
#   requireNamespace("diffusionMap")
#   requireNamespace("stats")
#   space <- diffusionMap::diffuse(stats::as.dist(SCORPIUS::correlation_distance(x)), neigen=neigen)
#   process_dimred(space$X[,seq_len(ndim)], rownames(x))
# }

dimred_ica <- function(x, ndim = 3) {
  requireNamespace("fastICA")
  space <- fastICA::fastICA(t(scale(t(x))), ndim)$S
  process_dimred(space, rownames(x))
}

dimred_lle <- function(x, ndim = 3) {
  requireNamespace("lle")
  poss_k <- lle::calc_k(t(scale(t(x))), ndim)
  k <- poss_k$k[which.min(poss_k$rho)]
  space <- lle::lle(t(scale(t(x))), ndim, k)$Y
  process_dimred(space, rownames(x))
}

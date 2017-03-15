clean_dimred_output <- function(space, rownames) {
  space <- as.matrix(space)
  dimnames(space) <- list(rownames, paste0("Comp", seq_len(ncol(space))))
  space
}

# TODO: implement LMDS
# lmds = function(x, ndim=3) {
#   mds.out <- dambiutils::mds_withlandmarks(x %>% as.data.frame, SCORPIUS::correlation.distance, k = ndim, landmark.method = "naive", num.landmarks = min(1000, round(nrow(x)*0.1)), num.seed.landmarks = 10, pca.normalisation = F)
#   process_dimred(mds.out$S)
# }

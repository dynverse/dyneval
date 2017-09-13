# lmds = function(x, ndim=3) {
#   mds.out <- dambiutils::mds_withlandmarks(x %>% as.data.frame, SCORPIUS::correlation_distance, k = ndim, landmark.method = "naive", num.landmarks = min(1000, round(nrow(x)*0.1)), num.seed.landmarks = 10, pca.normalisation = F)
#   process_dimred(mds.out$S)
# }

# simlr = function(x, ndim=3, nclusters=4) {
#   result = SIMLR::SIMLR(t(x), nclusters)
#   S = result$S
#   space = tsne::tsne(as.dist(max(S)-S), k = ndim)
#   process_dimred(space)
# }

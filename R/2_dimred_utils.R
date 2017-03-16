clean_dimred_output <- function(space, rownames) {
  space <- as.matrix(space)
  dimnames(space) <- list(rownames, paste0("Comp", seq_len(ncol(space))))
  space
}

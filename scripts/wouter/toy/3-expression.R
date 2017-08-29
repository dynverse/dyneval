library(tidyverse)
library(dyneval)

source("scripts/wouter/toy/generation.R")
source("scripts/wouter/toy/perturbation.R")

# generate a simple linear toy dataset
task <- generate_linear(ncells = 50)

timepoints <- task$progressions$percentage
cell_ids <- task$progressions$cell_id

ngenes <- 50
splinefuns <- map(seq_len(ngenes), function(gene) {
  x <- c(0, runif(sample(seq_len(1))), 1)
  y <- runif(length(x), -1, 1)

  interpol <- approxfun(x, y)

  function(z) {
    interpol(z)
  }
})
# coefs <- runif(ngenes, -1, 1)
# intercepts <- runif(ngenes, -2, 2)

# expression <- t(coefs %*% t(timepoints)) %>% apply(1, function(x) x + intercepts) %>% t
expression <- map(splinefuns, ~.(timepoints)) %>% invoke(rbind, .) %>% t
expression <- expression + rnorm(length(expression), mean = 0, sd=0.01)
dimnames(expression) <- list(cell_ids, seq_len(ncol(expression)))
expression[order(timepoints),] %>% t %>% pheatmap::pheatmap(cluster_cols=F)
plot(timepoints, expression[, 1])

counts <- (2^expression)*10
counts[order(timepoints),] %>% t %>% pheatmap::pheatmap(cluster_cols=F)
task$counts <- counts
tasks <- dyneval:::list_as_tibble(list(task))

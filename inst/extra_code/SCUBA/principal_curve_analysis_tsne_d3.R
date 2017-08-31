#
# R file to do principal curve analysis on tsne output
#

library(princurve)


args <- commandArgs()

#ARGS
baseName <- args[6]


infile <- paste0(baseName, "/intermediate_files/tsne_d3.csv")
pcvout <- paste0(baseName, "/intermediate_files/tsne_d3_pcv.csv")
lambdaout <- paste0(baseName, "/intermediate_files/tsne_d3_lambda.csv")

x <- read.csv(file=infile, header=TRUE, sep=",")

y <- data.matrix(x)

# code for not initializing the starting point lowess version
fitpc <- principal.curve(y, plot = TRUE, smoother = "lowess", maxit = 200)

write.table(file=pcvout, sep=",", fitpc$s)
write.table(file=lambdaout, sep=",", fitpc$lambda)


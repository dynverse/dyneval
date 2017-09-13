#BiocInstaller::biocLite("cellTree")


library(cellTree)
library(SCORPIUS)
library(tidyverse)
library(igraph)
library(magrittr)

expression = ginhoux$expression


.datasets_location = "/home/wouters/thesis/projects/dyngen/results/"
datasetidoi = dataset_id
dataset = dyngen::load_dataset(datasetidoi)
goldstandard = dataset$gs
expression = dataset$counts
cellinfo = goldstandard$cellinfo %>% slice(match(rownames(expression), cell))
milestone_percentages = goldstandard$milestone_percentages %>% slice(match(rownames(expression), id))


E = task$expression
#expression <- t(SCORPIUS::quant.scale(t(dataset$counts), 0))
lda.results = compute.lda(t(E) + min(E) + 1, k.topics=4, method="maptpx", log.scale=F, sd.filter=0)
#dists = get.cell.dists(lda.results)
mst.tree = compute.backbone.tree(lda.results, width.scale.factor=1.5)
#mst.tree = compute.backbone.tree(lda.results, only.mst = T)
#ct.plot.topics(mst.tree)
#ct.plot.grouping(mst.tree)

layout = cellTree:::.compute.tree.layout(pred_outputs$celltree$mst.tree, 5)
layout = cellTree:::.compute.tree.layout(mst.tree, 5)
plot(layout, vertex.size=2, vertex.shape="circle", vertex.color=cellinfo$piecestateid, vertex.label="")
plot(layout, vertex.size=2, vertex.shape="circle", vertex.color=milestone_percentages[, 1] > 0, vertex.label="")
plot(layout, vertex.size=2, vertex.shape="circle", vertex.color=viridis::viridis(100)[pred_outputs$celltree$state_percentages[, 1]*100], vertex.label="")


backbone = induced_subgraph(mst.tree, get.vertex.attribute(mst.tree, "is.backbone"))
tomerge = names(V(backbone))[degree(backbone) == 2]
backbone = backbone %>% as_data_frame
for(node in tomerge) {
  subgraph = backbone %>% filter((from == node) | (to == node))
  includeds = unlist(subgraph$included)
  newnodes = subgraph %>% {c(.$from, .$to)} %>% keep(~.!=node)

  backbone = backbone %>% filter((from != node) & (to != node)) %>% bind_rows(list(from=newnodes[[1]], to=newnodes[[2]], included=list(c(includeds, node))))
}
backbone = backbone %>% tidyr::unnest()

backbonenodes = names(V(mst.tree))[get.vertex.attribute(mst.tree, "is.backbone")]
sidenodes = names(V(mst.tree))[!get.vertex.attribute(mst.tree, "is.backbone")]
centralnodes = unique(c(backbone$from, backbone$to))
names(centralnodes) = seq_along(centralnodes)
sidenodes2backbone = mst.tree %>% as_data_frame %>% filter(to %in% sidenodes) %$% set_names(from, to)

percentages = tibble()
for(node in names(V(mst.tree))) {
  if (node %in% names(sidenodes2backbone)) {
    realnode = as.character(sidenodes2backbone[[node]])
  } else {
    realnode = node
  }

  if(realnode %in% centralnodes) {
    percentages = percentages %>% bind_rows(tibble(milestone=as.character(which(centralnodes == realnode)), cell=node, percentage=1))
  } else {
    centralnodesoi = backbone %>% filter(included == realnode) %>% {c(.$from, .$to)}
    distances = igraph::distances(mst.tree, realnode, centralnodesoi)
    percentages = percentages %>% bind_rows(tibble(milestone=as.character(match(centralnodesoi, centralnodes)), cell=node, percentage=1-distances[1, ]/sum(distances)))
  }
}
percentages = reshape2::acast(percentages %>% mutate(cell=as.numeric(cell)) %>% arrange(cell), cell~milestone, value.var="percentage", fill=0) %>%
  as.data.frame() %>% mutate(id=rownames(expression))
pheatmap::pheatmap(t(percentages), cluster_cols=F, cluster_rows=F, annotation_col = goldstandard$cellinfo %>% dplyr::select(piecestateid) %>% as.data.frame())

plot(mst.tree, vertex.color = c("red", "blue")[(percentages[4, ] > 0) + 1], vertex.shape="circle", layout=layout.mds(mst.tree), vertex.label="", vertex.size=3)


dists = SCORPIUS::correlation_distance(expression)
dists = SCORPIUS::euclidean_distance(t(scale(t(expression))))
space = SCORPIUS::reduce_dimensionality(dists)
SCORPIUS::draw.trajectory.plot(space, cellinfo$piecestateid)
SCORPIUS::draw.trajectory.plot(space, percentages[3, ])












node = "500"
distances = igraph::distances(mst.tree, node, backbone %>% filter(included == node) %>% {c(.$from, .$to)})



milestones = degree(backbone) != 2

vertex_attr_names(mst.tree)


degree(mst.tree) %>% hist



vertexoi = get.vertex.attribute(mst.tree, "is.backbone")
plot(mst.tree, vertex.color = c("red", "blue")[as.numeric(vertexoi)+1], vertex.shape="circle")


vertexoi = get.vertex.attribute(mst.tree, "is.backbone")




dists = euclidean.distance(expression)
space = SCORPIUS::reduce.dimensionality(dists)
SCORPIUS::draw.trajectory.plot(space, goldstandard$cellinfo$piecestateid)
SCORPIUS::draw.trajectory.plot(space, percentages[3, ])





library(PRISM)
qsub_lapply(1:5, function(x) Sys.sleep(5), qsub_config = override_qsub_config(num_cores = 2))

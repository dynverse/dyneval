library(dyneval)
library(tidyverse)

task <- dyneval:::read_ti_task_data("linear_ginhoux")

task$state_network
task$state_percentages

# ggplot(task$state_percentages) + geom_point(aes(MDP, CDP, colour = PreDC)) + scale_colour_distiller(palette = "RdBu")
# pheatmap::pheatmap(t(as.matrix(task$state_percentages[,-1])))

output <- dyneval:::trainLearner.ti.scorpius(
  task, .subset = NULL, num_dimensions = 3, num_clusters = 4)

output$state_network
output$state_percentages

# todo:
#  - compare task$state_network and task$state_percentages
#  - to output$state_network and output$state_percentages
gold <- as.matrix(task$state_percentages[,-1])
pred <- as.matrix(output$state_percentages[,-1])
grid <- expand.grid(i = seq_len(nrow(gold)), j = seq_len(nrow(gold)))
grid$gold.mse <- rowMeans(abs(gold[grid$i,] - gold[grid$j,]))
grid$pred.mse <- rowMeans(abs(pred[grid$i,] - pred[grid$j,]))
cor(grid$gold.mse, grid$pred.mse)

ggplot() +
  geom_point(aes(output$state_percentages$A, task$state_percentages$MDP, colour = "MDP")) +
  geom_point(aes(output$state_percentages$A, task$state_percentages$CDP, colour = "CDP")) +
  geom_point(aes(output$state_percentages$A, task$state_percentages$PreDC, colour = "PreDC"))
ggplot(grid) + geom_point(aes(gold.mse, pred.mse))


#########################
# starting from here on #
#########################
library(monocle)

bla<-readRDS("data/ti/linear/ginhoux.rds")
HSM <- t(as.matrix(bla$expression[,1:3500]))
HSMM<-newCellDataSet(as.matrix(HSM))
#HSMM_myo <- exprs(HSMM[1:100,])

HSMM_myo1 <- reduceDimension(HSMM, max_components=2)
HSMM_myo <- orderCells(HSMM_myo1, reverse=FALSE)
plot_cell_trajectory(HSMM_myo)#, color_by="Hours")


##igraph of the MST:
gr <- HSMM_myo@auxOrderingData$DDRTree$pr_graph_cell_proj_tree
root <- HSMM_myo@auxOrderingData$DDRTree$root_cell
library(igraph)
plot(gr, vertex.label = NA, layout = layout_as_tree(gr, root = root))
edf <- as_data_frame(gr)

## identifying specific cells based on their conectivity:
deg<-degree(gr, v = V(gr), mode = c("all"))
special<-which(deg!=2) ## vertices with degree 1 (endpoints) and degree 3 (branch points)
branch_points<-which(deg==3)

asp<-all_simple_paths(gr, from=roott, to = special, mode = c("total")) ##all simple paths in tree

le<-c()
#ordering lengths:
for (i in 1:length(asp)){
  le<-c(le,length(asp[[i]]))
}
ord<-order(le)

# inititation (root to first branching point)
firs<-asp[[ord[1]]][1]
las<-asp[[1]][length(asp[[ord[1]]])]
vectt<-data.frame(from=firs$name,to=las$name,length=1) # filling state_network

nodelist<-list(asp[[1]][c(which(asp[[1]]==firs):which(asp[[1]]==las))])
j<-2
# completing the tree structure
for(i in ord[-1]){
  branch_vect<-which(asp[[i]]%in%branch_points)
  if(asp[[i]][length(asp[[i]])]$name==asp[[i]][branch_vect]$name[length(asp[[i]][branch_vect]$name)]){
    firstt<-asp[[i]][branch_vect][length(asp[[i]][branch_vect]$name)-1]
  } else {firstt<-asp[[i]][branch_vect][length(asp[[i]][branch_vect]$name)]}
  vectt<-rbind(vectt,data.frame(from=firstt$name,to=V(gr)$name[asp[[i]][length(asp[[i]])]],length=1))

  las<-which(asp[[i]]==asp[[i]][length(asp[[i]])])
  nodelist[[j]]<-asp[[i]][c(which(asp[[i]]==firstt):las)] ## contains the nodes of this branch
  j<-j+1
}

# compute shortest paths between special cells:
for (i in 1:dim(vectt)[1]){
  vectt[i,3]<-distances(gr, v = as.character(vectt[i,1]), to = as.character(vectt[i,2]), mode = c("in"), algorithm = ("automatic"))
}

state_network<-vectt


# computing only branch probabilities (2 probs per cell)?
special_names<-unique(c(as.character(vectt[,1]),as.character(vectt[,2])))
state_percentages<-matrix(rep(0,length(special_names)*length(V(gr))),nrow=length(V(gr)))
colnames(state_percentages)<-special_names
rownames(state_percentages)<-V(gr)$name

for ( i in 1:length(nodelist)){
  for (j in 1:length(nodelist[[i]])){
    row2fil<-which(rownames(state_percentages)==nodelist[[i]][j]$name)
    col1<-which(colnames(state_percentages)==as.character(vectt[i,1]))
    col2<-which(colnames(state_percentages)==as.character(vectt[i,2]))
    dist_node<-c(distances(gr,v=nodelist[[i]][1],to=nodelist[[i]][j]),
           distances(gr,v=nodelist[[i]][j],to=nodelist[[i]][length(nodelist[[i]])]))
    state_percentages[row2fil,col1]<-1-(dist_node[1]/sum(dist_node))
    state_percentages[row2fil,col2]<-1-(dist_node[2]/sum(dist_node))
  }
}

head(state_percentages)

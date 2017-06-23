library(dyneval)
library(tidyverse)

task <- dyneval:::read_ti_task_data("linear_ginhoux")

task$state_network
task$state_percentages

library(destiny)
library(rgl)
library(ggplot2)
library(ade4)

### load data ###
setwd("~/Documents/VIB/Projects/dyneval")
bla<-readRDS("data/ti/linear/ginhoux.rds")
ginhoux_data<-as.matrix(bla$expression)
pca<-dudi.pca(ginhoux_data,scale=T,center=T,scannf = F,nf=50) # pca, keep 50 main PCs
dm<-DiffusionMap(pca$li) # diffusion map on 50 main PCs
plot(dm,1:2)
dp<-DPT(dm)
plot(dp)#,col_by="branch")
summary(as.factor(dp$branch))

# find terminal points:
term<-find_tips(dp)

## state_network ##
state_network<-data_frame(from=rownames(ginhoux_data)[term[1]],to=rownames(ginhoux_data)[term[2]], dist=dm@transitions[term[1],term[2]])
state_network[2,]<-data_frame(from=rownames(ginhoux_data)[term[1]],to=rownames(ginhoux_data)[term[3]],
                              dist=dm@transitions[term[1],term[3]])
state_network[3,]<-data_frame(from=rownames(ginhoux_data)[term[2]],to=rownames(ginhoux_data)[term[3]],
                              dist=dm@transitions[term[2],term[3]])

## state_percentage ##
state_percentage<-matrix(rep(0, nrow(ginhoux_data)*3), ncol=3)

for (i in 1:nrow(ginhoux_data)){
  state_percentage[i,]<-c(dm@transitions[i,term[1]],dm@transitions[i,term[2]],dm@transitions[i,term[3]])
  row_sum<-sum(state_percentage[i,])
  state_percentage[i,]<-state_percentage[i,]/row_sum
}
rownames(state_percentage)<-rownames(ginhoux_data)
colnames(state_percentage)<-c(rownames(ginhoux_data)[term[1]],rownames(ginhoux_data)[term[2]],
                              rownames(ginhoux_data)[term[3]])

## taking into account which branch DPT assigned the cell to :
br[br==1]<-dp_tips[1]

## nb of branches:
length(unique(dp$Branch[!is.na(unique(dp$Branch))]))


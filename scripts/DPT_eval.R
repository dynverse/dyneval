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

## state_network ##
branches<-unique(dp$branch)
#br<-branches[complete.cases(branches),]
dp_tips<-which(dp@tips==T)
br[br==1]<-dp_tips[1]

dim(branches)

# find terminal points:
term<-find_tips(dp)

length(unique(dp$Branch[!is.na(unique(dp$Branch))]))

#tips of branch 2:
trues<-which(dp@tips==T)
plot(dp,col_by="branch") # shows only 3 tips

ggplot(data=dp)+
  geom_point(mapping=aes(x=DC1,y=DC2))

ggplot(data=dp)+
  geom_smooth(mapping=aes(x=DC1,y=DC2,linetype=branch)) # doesn't find Branch !!


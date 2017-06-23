library(tidyverse)
if(!exists("foo", mode="function")) source("~/Documents/VIB/Projects/dyneval/scripts/waterfall.R")

setwd("~/Documents/VIB/Projects/dyneval")
bla<-readRDS("data/ti/linear/ginhoux.rds")
ginhoux_data<-as.matrix(bla$expression)

wf<-waterfall(ginhoux_data,k=10)

## state_network ##
state_network<-data_frame(from=names(which(wf$value==min(wf$value))),
                          to=names(which(wf$value==max(wf$value))),
                          dist=1)

## state_percentages ##
state_percentages<-matrix(rep(0,length(wf$value)*2),ncol=2)
state_percentages[,1]<-1-(as.numeric(wf$value))
state_percentages[,2]<-as.numeric(wf$value)
colnames(state_percentages)<-c(names(which(wf$value==min(wf$value))),
                               names(which(wf$value==max(wf$value))))
rownames(state_percentages)<-names(wf$value)

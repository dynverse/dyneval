library(tidyverse)
if(!exists("foo", mode="function")) source("~/Documents/VIB/Projects/dyneval/scripts/waterfall.R")

setwd("~/Documents/VIB/Projects/dyneval")
bla<-readRDS("data/ti/linear/ginhoux.rds")
ginhoux_data<-as.matrix(bla$expression)

wf<-waterfall(ginhoux_data,k=10)

## milestone_network ##
milestone_network<-data_frame(from=names(which(wf$value==min(wf$value))),
                          to=names(which(wf$value==max(wf$value))),
                          dist=1)

## milestone_percentages ##
milestone_percentages<-matrix(rep(0,length(wf$value)*2),ncol=2)
milestone_percentages[,1]<-1-(as.numeric(wf$value))
milestone_percentages[,2]<-as.numeric(wf$value)
colnames(milestone_percentages)<-c(names(which(wf$value==min(wf$value))),
                               names(which(wf$value==max(wf$value))))
rownames(milestone_percentages)<-names(wf$value)

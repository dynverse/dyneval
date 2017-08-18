pseudotime.foo <-function(gene_name, all=all, span=.75, unit=.025, all.col=all.col, ...){

  pdf(file=paste0(gene_name,"_pseudotime_exp.pdf"),width=5,height=4)
  pseudotime_gene_plot(gene_name,span=span, all=all, all.col=all.col, ...)
  dev.off()

  # HMM
  library(RHmm)
  exp <-as.numeric(all[grep(gene_name,rownames(all)),match(rownames(pseudotime.df),colnames(all))])
  onoff <-state.foo(exp,unit=unit)#; onoff[which(onoff==2/3)] <-1; onoff[which(onoff==4/3)] <-2

  pdf(file=paste0(gene_name,"_pseudotime_HMM.pdf"),width=8,height=.5)
  pheatmap(t(as.matrix(onoff)),cluster_row=F,cluster_col=F,color=c("#000000","#FFEE00"),legend=F,show_rownames=F,show_colnames=F,cellwidth=6.25,cellheight=6.25)
  dev.off()
}


pseudotime_gene_plot <-function(gene_name,span, all=all, all.col=all.col, ...){
  xpt <-as.numeric(pseudotime.df$pseudotime)
  ypt <-as.numeric(all[,match(rownames(pseudotime.df),colnames(all))][grep(gene_name,rownames(all)),])

  y.loess <-loess(ypt~xpt,span=span)
  ypt.predict <-predict(y.loess,data.frame(x=xpt),se=T)
  plot(pseudotime.df$pseudotime,ypt,bty="n",xlab="Pseudotime",ylab="Expression (TPM)",pch=20,cex=1,col=paste0(all.col[rownames(pseudotime.df)],"99"),main=gene_name,...)
  lines(xpt,ypt.predict$fit,lwd=1,col="#BC2B2B")

  y.polygon <- c((ypt.predict$fit+1.96*ypt.predict$se.fit), rev(c(ypt.predict$fit-1.96*ypt.predict$se.fit)))
  x.polygon <- c(xpt, rev(xpt))
  polygon(x.polygon, y.polygon, col="#00529533", border=NA)
}

sin_scale2.foo <-function(m,limit=200){
  scale_row.foo <-function(X){
    medX <-max(median(X),200)
    maxX <-max(max(X),2*200)
    Y <-rep(0,length(X))
    l.idx <-which(X>=medX)
    s.idx <-which(X<medX)

    Y[s.idx] <- 0.5*sin((X[s.idx]-medX)/medX*pi/2)+.5
    Y[l.idx] <- 0.5*sin((X[l.idx]-medX)/(max(X)-medX)*pi/2)+.5

    return(Y)
  }
  t(apply(m,1,scale_row.foo))
}

scale_row.foo <-function(X,limit=200){
  medX <-max(median(X),limit)
  maxX <-max(max(X),2*limit)
  Y <-rep(0,length(X))
  l.idx <-which(X>=medX)
  s.idx <-which(X<medX)

  Y[s.idx] <- 0.5*sin((X[s.idx]-medX)/medX*pi/2)+.5
  Y[l.idx] <- 0.5*sin((X[l.idx]-medX)/(max(X)-medX)*pi/2)+.5

  return(Y)
}

marker_plot.foo <-function(gene_name){
  pdf(file=paste0("marker_",gene_name,".pdf"),height=5.3,width=5)
  scala_exp <-as.numeric(all[grep(gene_name,rownames(all)),])
  size_ramp <-scale_row.foo(scala_exp)*5+.01
  plot(pca$x[,1:2],col=paste0(all.col,"90"),cex=size_ramp,pch=19)
  dev.off()
}

mst.of.classification <- function (x, k=6,color,seed=1,...) {
  #   x <- t(x)
  #   x <- t( t(x) - apply(x,2,mean) )
  # c <- cor(x, method="pearson"); d <- dist(c); hr <- hclust(d, method = "ward", members=NULL)
  pca <- prcomp(as.data.frame(t(x)), cor=T)
  y <- pca$x

  r <- prcomp(t(x))

  # kmeans
  set.seed(seed)
  r <- kmeans(y[,1:2],k)
  z <- r$centers
  z <- z[order(z[,1]),]
  rownames(z) <-paste0("t",1:nrow(z))
  m <- mst(dist(z))
  plot(y[,1:2], col=paste0(color), cex=3, pch=20, bty="n", ...)
  points(z[,1:2], col="#FF0000", pch=19,cex=1)
  #    text(y[,1:2],labels=rownames(y))
  w <- which(m!=0)
  i <- as.vector(row(m))[w]
  j <- as.vector(col(m))[w]
  segments( z[i,1], z[i,2], z[j,1], z[j,2], col="#FF0000" ,lwd=5)
}

intersect.foo <-function(tot1,tot2,ab){
  # intersect between tot1-tot2 segment and ab point
  x1 = tot1[1];y1=tot1[2];x2=tot2[1];y2=tot2[2];a=ab[1];b=ab[2]
  slope = (y2-y1)/(x2-x1)

  int.x=(slope*x1 + a/slope + b -y1)/(slope+1/slope)
  int.y=(-x1 + a + b*slope + y1/slope)/(slope+1/slope)
  return(c(int.x,int.y))
}

inside_check.foo <-function(tot1,tot2,ab){
  int <-intersect.foo(tot1,tot2,ab)
  all((tot1<=int& int<=tot2)|(tot2<=int& int<=tot1))
}

distance.foo <-function(tot1,tot2,ab){
  int <-intersect.foo(tot1,tot2,ab)
  return(dist(rbind(int,ab)))
}

`%notin%` <- function(x,y) !(x %in% y)
new.foo <-function(x,y){x[which(x %notin% y)]}

unit_vector.foo <-function(x,y){
  x <-as.numeric(x);y <-as.numeric(y)
  (y-x)/(dist(rbind(x,y)))
}

crossvec <- function(x,y){
  if(length(x)!=2 |length(y)!=2) stop('bad vectors')
  cv <-  x[1]*y[2]-x[2]*y[1]
  return(invisible(cv))
}

crossvec_direction <-function(tot1,tot2,ab){
  if ((crossvec(tot2-tot1,ab-tot1))>=0){
    return(1)
  } else{
    return(-1)
  }
}

pseudotimeprog.foo <- function(x, k=10, color, x.reverse=F, plot=F, ...){
  r <- prcomp(t(x))
  y <- r$x*matrix(r$sdev^2/sum(r$sdev^2),nrow=nrow(r$x),ncol=ncol(r$x),byrow=T)
  #y <-y[order(y[,1]),]
  #u <- r$rotation

  # kmeans
  r <- kmeans(y,k)
  z <- r$centers
  z <- z[order(z[,1]),]
  rownames(z) <-paste0("t",1:nrow(z))
  m <- mst(dist(z))

  t.names <-names(which(colSums(m!=0)==1))[1] # There are two ends, then use the left most one.
  for (i in 1:nrow(m)){
    t.names <-append(t.names,names(which(m[t.names[i],]==1))[which(names(which(m[t.names[i],]==1)) %notin% t.names)])
  }

  y2d <-y[,1:2]
  #y2d <-y2d[order(y2d[,1]),]
  z2d <-z[,1:2]
  z2d <-z2d[t.names,]

  time_start.i <-0
  updatethis.dist <-rep(Inf,nrow(y2d))
  updatethis.time <-rep(0,nrow(y2d))
  update.updown <-rep(0,nrow(y2d))
  pseudotime.flow <-c(0)

  for (i in 1:(nrow(z2d)-1)){
    # distance between this z2d.i and all y2d
    dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})

    # distance between this z2d.i-z2d.i+1 segment and "insider" y2d
    inside_this_segment <-which(apply(y2d,1,function(X){inside_check.foo(z2d[i,],z2d[i+1,],X)}))
    seg.dist.i <-rep(Inf,nrow(y2d))
    seg.dist.i[inside_this_segment] <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})[inside_this_segment]

    # intersect coordinate between this z2d.i-z2d.i+1 segment and all y2d
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))

    # this z2d.i-z2d.i+1 segment's unit vector
    seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])

    # UPDATE
    # 2. idx for the shortest distance at this round (either dot or seg)
    update.idx <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,which.min)
    # 3. update the pseudotime for y2ds with the short distance from the z2d.i
    updatethis.time[which(update.idx==1)] <-time_start.i
    # 4. update the pseudotime for y2ds with the short distance from the z2d.i-z2d.i+1 segment
    relative_cordinates <-t(apply(intersect.i[which(update.idx==2),],1,function(X){seg_unit_vector%*%(X-z2d[i,])}))
    updatethis.time[which(update.idx==2)] <-time_start.i + relative_cordinates
    # 1. update the shortest distance
    updatethis.dist <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,min)

    update.updown[which(update.idx==1)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*dot.dist.i)[which(update.idx==1)]
    update.updown[which(update.idx==2)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[which(update.idx==2)]

    # update time for the next round
    time_start.i <-time_start.i + dist(rbind(z2d[i,],z2d[i+1,]))
    pseudotime.flow <-append(pseudotime.flow,time_start.i)
  }

  # For the y2ds that are closest to the starting z2d
  i=1
  dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
  if (length(start.idx <-which(dot.dist.i <= updatethis.dist))>0) {
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))
    seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])
    relative_cordinates <-0 + t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])}))[start.idx]
    updatethis.time[start.idx] <-relative_cordinates
    seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})
    update.updown[start.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[start.idx]
  }
  # For the y2ds that are closest to the arriving z2d
  i=nrow(z2d)
  dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
  if (length(arrive.idx <-which(dot.dist.i <= updatethis.dist))>0) {
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i-1,],z2d[i,],X)}))
    seg_unit_vector <-unit_vector.foo(z2d[i-1,],z2d[i,])
    relative_cordinates <-time_start.i + as.numeric(t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])})))[arrive.idx]
    updatethis.time[arrive.idx] <-relative_cordinates
    seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i-1,],z2d[i,],X)})
    update.updown[arrive.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i-1,],z2d[i,],X)})*seg.dist.i)[arrive.idx]
  }

  pseudotime <-updatethis.time
  pseudotime.y <-update.updown
  pseudotime.flow <-pseudotime.flow

  if (x.reverse){
    pseudotime <- -pseudotime
    pseudotime.flow <- -pseudotime.flow
  }

  pseudotime_range <-max(pseudotime)-min(pseudotime)

  pseudotime.flow <-pseudotime.flow-min(pseudotime)
  pseudotime.flow <-pseudotime.flow/pseudotime_range

  pseudotime <-pseudotime-min(pseudotime)
  pseudotime <-pseudotime/pseudotime_range

  if (plot) {
    plot(pseudotime, pseudotime.y,col=paste0(color),cex=3,pch=20,bty="n", ...)
    points(pseudotime.flow,rep(0,length(pseudotime.flow)), col="#FF0000", pch=19,cex=1)
    segments(pseudotime.flow[1], 0, pseudotime.flow[length(pseudotime.flow)], 0, col="#FF0000" ,lwd=5)
  }

  return(data.frame(y[,1:2], pseudotime, pseudotime.y))
}

state.foo <-function(X,s_from=2,unit=.025){
  set.seed(1)
  s=s_from # number of states

  time <-pseudotime.df$pseudotime
  time.unit <-0:max(round((time-min(time))/unit))

  onoff <-rep(.5,length(time.unit))


  exp <-as.numeric(X)

  onoff_determine.foo <-function(exp,s){
    obs <-unlist(lapply(split(exp,round((time-min(time))/unit)),mean))
    ResFit <- HMMFit(obs, nStates=s)# Baum-Welch
    state <-viterbi(ResFit,obs)
    names(onoff) <-time.unit
    onoff[match(as.character(names(obs)),names(onoff))] <-state[[1]]

    onoff[which(onoff==.5)] <-onoff[which(onoff==.5)+1]
    onoff[which(onoff==.5)] <-onoff[which(onoff==.5)-1]
    onoff[which(onoff==.5)] <-onoff[which(onoff==.5)+2]
    onoff[which(onoff==.5)] <-onoff[which(onoff==.5)-2]

    state_order <-order(ResFit$HMM$distribution$mean)
    onoff <-sapply(onoff,function(X){state_order[X]})
  }

  onoff <-rep(.5,length(time.unit))
  onoff <-onoff_determine.foo(exp,2)

  endpoint=length(onoff)
  onethird=round(length(onoff)*1/3)
  twothird=round(length(onoff)*2/3)

  if (min(onoff)==max(onoff)| # If there is only one state exists
      !any(which(onoff==2) %in% (which(onoff==2)+1))| # No continuous state exists
      (cor(time,exp)> 0.2& length(which(onoff[1:twothird]==2))/2>length(which(onoff[(twothird+1):endpoint]==2)))| # If state estimation is way too off
      (cor(time,exp)< -.1& length(which(onoff[1:onethird]==2))<length(which(onoff[(onethird+1):endpoint]==2))/2) # If state estimation is way too off
  ){
    onoff <-rep(1,length(onoff))
    onoff <-rep(.5,length(time.unit))
    onoff <-onoff_determine.foo(exp,4)
  }

  if (min(onoff)==max(onoff)| # If there is only one state exists
      !any(which(onoff==2) %in% (which(onoff==2)+1))| # No continuous state exists
      (cor(time,exp)> 0.2& length(which(onoff[1:twothird]==2))/2>length(which(onoff[(twothird+1):endpoint]==2)))| # If state estimation is way too off
      (cor(time,exp)< -.1& length(which(onoff[1:onethird]==2))<length(which(onoff[(onethird+1):endpoint]==2))/2) # If state estimation is way too off
  ){
    onoff <-rep(1,length(onoff))
    onoff <-rep(.5,length(time.unit))
    onoff <-onoff_determine.foo(exp,3)
  }

  onoff <-onoff/max(onoff)*2
  onoff <-as.integer(ceiling(onoff))
  return(onoff)
}

diff_gene.foo <-function(exp.df=all,group.name1,group.name2){
  groupA.idx <-which(colnames(exp.df) %in% group.name1)
  groupB.idx <-which(colnames(exp.df) %in% group.name2)

  meanA <-rowMeans(exp.df[,groupA.idx])
  meanB <-rowMeans(exp.df[,groupB.idx])
  maxmean <-apply(cbind(meanA,meanB),1,max)

  # coefficient of variation
  covar1 <-apply(exp.df[,groupA.idx],1,sd)/meanA
  covar2 <-apply(exp.df[,groupB.idx],1,sd)/meanB

  test.wilcox <-function(X){
    return(wilcox.test(X[groupA.idx],X[groupB.idx])$p.value)}
  test.diff <-meanA/meanB

  pval <-apply(exp.df,1,test.wilcox)
  log2FC <-log2(test.diff)

  #sig.idx <-which(pval<.05& log2FC>0& maxmean>50& covar1<=1& meanA>=50)

  return(
    data.frame(name=rownames(exp.df),pval=pval,log2FC=log2FC,covar1=covar1,covar2=covar2,meanA=meanA,meanB=meanB)
  )
}

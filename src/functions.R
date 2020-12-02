genes.plot.pseudo0 <- function (object,genes.use, pseudo="pseudo.time",pch.use = 16, cex.use = 1.5,
                               spline.span = 0.75,name.x=NULL,col.use=NULL,do.logit=FALSE,
                               inset=c(0,0),ylim=NULL,do.label=T,lwd=1,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }

  cell.ids = set.ifnull(cell.ids, object@cell.names)
  xlab=paste(name.x,"Pseudotime")
  data.use = data.frame(t(FetchData(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.scaled = FALSE)))
  pseudo.order=order(data.use[pseudo,cell.ids])
  pseudo.data=unlist(data.use[pseudo,cell.ids][pseudo.order])
  object@scale.data=object@scale.data[,pseudo.order]
  ylim=set.ifnull(ylim,range(data.use[genes.use,cell.ids]))
  
  plot(0,0,xlim = range(pseudo.data),ylim = ylim,type = "l",xlab=xlab,ylab="log(FPKM)",
       cex = cex.use, pch = pch.use, font.lab=2,cex.lab=1.2,bty="l",...)
  col.use=set.ifnull(col.use,rainbow(length(genes.use)))
  
  for (i in 1:length(genes.use)){
    g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
    loess.fit = loess(g2 ~ pseudo.data, span = spline.span)
    lines(pseudo.data, loess.fit$fitted,type="l",col=col.use[i],lwd=lwd)
  }
  if(do.logit){
    gene.norm=apply(object@data[genes.use[i],cell.ids],1, FUN = function(X) (X - min(X))/diff(range(X)))
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm[pseudo.order]~pseudo.data,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn)]
    turn.mat=c(turn.mat,turn.point)
    abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
  }
  if(do.label){
    legend("topright",legend = genes.use, 
           lty=1, col=col.use, cex=.75,xpd=T,inset = inset)
  }
}

genePlot.pseudo0 <- function (object,pseudo="pseudo.time", gene, cell.ids = NULL, col.use = NULL, 
                             pch.use = 16, cex.use = 1.5, do.ident = FALSE, conf.int=FALSE,
                             do.spline = FALSE,spline.span = 0.75,name.x=NULL, do.logit=FALSE,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = use.scale)))
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  name.x = set.ifnull(name.x, pseudo)
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = as.numeric(data.use[gene, cell.ids])
  ident.use = as.factor(object@ident[cell.ids])
  pseudo.data=as.matrix(FetchData(object,pseudo))
  
  if(do.logit){
    gene.norm=apply(object@data[gene,colnames(object@data)],1, FUN = function(X) (X - min(X))/diff(range(X)))
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm~pseudo.data,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn),1]
  }
  
  if (length(col.use) > 1) {
    col.use = col.use[as.numeric(ident.use)]
  }
  else {
    col.use = set.ifnull(col.use, as.numeric(ident.use))
  }

  plot(g1, g2, xlab = name.x, ylab = gene, col = col.use, cex = cex.use, 
       main = "", pch = pch.use, font.lab=2,cex.lab=1.2, ...)
  
  if (do.logit){
    abline(v=turn.point,lwd=4,col="purple",lty="longdash")
  }
  
  if (do.spline) {
    loess.fit = loess(g2 ~ g1, span = spline.span)
    lines(g1[order(pseudo.data)], loess.fit$fitted[order(pseudo.data)],col = "black",lwd=3)
    if(conf.int){
      prid = predict(loess.fit,se=T)
      lines(g1[order(g1)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(g1)],col = "darkblue",lwd=1.2,lty=2)
      lines(g1[order(g1)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(g1)],col = "darkblue",lwd=1.2,lty=2)
    }
  }
}

genes.plot.pseudo <- function (object,genes.use, cell.ids = NULL, pseudo="pseudo.time",pch.use = 16, cex.use = 1.5, do.induct=TRUE,
                               do.scale=FALSE,use.imputed = FALSE, method = "logit", spline.span = 0.75,name.x=NULL,do.line=FALSE,col.use=NULL, bin = 25, nState = 2,
                               do.smooth=FALSE,knn=10,inset=c(0,0),ylim=NULL,conf.int=FALSE,do.label=F,do.ret.mat=F,lwd=1,do.ret.pt=F,use.loess=F,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  nn.count.add <- function(mat, k, dmat) {
    dmat[dmat == 0] <- Inf
    knn <- t(apply(dmat, 1, function(x) rank(x) <= k))
    mat <- mat + mat %*% t(knn)
    return(mat)
  }
  
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  ident.use = as.factor(object@ident[cell.ids])
  xlab=paste(name.x,"Pseudotime")
  ret.mat=c()
  turn.mat=c()
  turn.pt=list(on=c(),off=c())
  if(do.induct & method == "hmm") {
    library(RHmm)
    object.pseudo = pseudoBin(object, bin.size = bin)}
  
  if (do.smooth){
    pt.ncc=fetch.data(object,pseudo)
    dmat = as.matrix(dist(pt.ncc))
    count.mat = object@scale.data
    object@scale.data = nn.count.add(count.mat, knn, dmat)
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }
  
  if(do.scale){
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }else{
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = FALSE)))
  }
  
  pseudo.order=order(data.use[pseudo,cell.ids])
  pseudo.data=unlist(data.use[pseudo,cell.ids][pseudo.order])
  object@scale.data=object@scale.data[,pseudo.order]
  ylim=set.ifnull(ylim,range(data.use[genes.use,cell.ids]))
  
  plot(0,0,xlim = range(pseudo.data),ylim = ylim,type = "l",xlab=xlab,ylab="log(FPKM)",
       cex = cex.use, pch = pch.use, font.lab=2,cex.lab=1.2,bty="l",...)
  col.use=set.ifnull(col.use,rainbow(length(genes.use)))
  if(do.line){
    for (i in 1:length(genes.use)){
      g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
      points(pseudo.data, g2,col=col.use[i],pch=20)
      lines(pseudo.data, g2,col=col.use[i],lwd=lwd)
    }
  }else{
    for (i in 1:length(genes.use)){
      g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
      loess.fit = loess(g2 ~ pseudo.data, span = spline.span)
      lines(pseudo.data, loess.fit$fitted,type="l",col=col.use[i],lwd=lwd)
      ret.mat=rbind(ret.mat,loess.fit$fitted)
      
      if(conf.int){
        prid = predict(loess.fit,se=T)
        lines(pseudo.data[order(pseudo.data)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(pseudo.data)],col=col.use[i],lty=2)
        lines(pseudo.data[order(pseudo.data)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(pseudo.data)],col=col.use[i],lty=2)
      }
      
      if(do.induct & method == "logit"){
        if(do.smooth){
          gene.norm=apply(matrix(object@scale.data[genes.use[i],cell.ids]),2, FUN = function(X) (X - min(X))/diff(range(X)))
        }
        
        if(use.loess){
          logit.use <- loess.fit$fitted
          gene.norm=(logit.use - min(logit.use))/diff(range(logit.use))}
        
        if(!use.loess & !do.smooth){
          gene.norm=apply(object@data[genes.use[i],cell.ids],1, FUN = function(X) (X - min(X))/diff(range(X)))
        }
        
        gene.norm[which(gene.norm>0.5)]=1
        gene.norm[which(gene.norm<0.5)]=0
        model <- glm(gene.norm[pseudo.order]~pseudo.data,family=binomial(link='logit'))
        model.fit=abs(fitted(model)-0.5)
        turn=min(model.fit)
        turn.point=pseudo.data[which(model.fit==turn)]
        turn.mat=c(turn.mat,turn.point)
        abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
      }
      if(method == "hmm"){
        if(use.loess){hmm.use <- loess.fit$fitted} else {hmm.use <- unlist(object@data[genes.use[i],])}
        if(max(hmm.use) <= 2){
          turn=1
          turn.point=pseudo.data[turn]
          names(turn.point)=genes.use[i]
          abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
          turn.pt$off=c(turn.pt$off,turn.point)
        } else {
          hmm.fit=HMMFit(hmm.use,nStates = nState)
          hmm.vit=viterbi(hmm.fit,hmm.use)
          if(length(unique(hmm.vit$states)) == 1){
            turn=1
            turn.point=pseudo.data[turn]
            names(turn.point)=genes.use[i]
            abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
            if(loess.fit$fitted[turn] >=1){
              turn.pt$on=c(turn.pt$on,turn.point)
            } else{
              turn.pt$off=c(turn.pt$off,turn.point)
            }
          } else {
            counter = 1
            turn = NULL
            for (j in 1:(length(hmm.use)-1)){
              if(hmm.vit$states[j] - hmm.vit$states[j+1] != 0){
                turn.point=pseudo.data[j]
                names(turn.point)=genes.use[i]
                abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
                if(hmm.use[j]-hmm.use[j+1]>0){
                  turn.pt$off = c(turn.pt$off, turn.point)
                } else {turn.pt$on=c(turn.pt$on, turn.point)}
              }
            }
          }
        }
      }
    }
  }
  if(do.label){
    legend("topright",legend = genes.use, 
           lty=1, col=col.use, cex=.75,xpd=T,inset = inset)
  }
  if(do.ret.mat){ 
    rownames(ret.mat)=genes.use
    colnames(ret.mat)=cell.ids[pseudo.order]
    return(ret.mat)}
  if(do.ret.pt){
    if(method=="hmm") {return(turn.pt)} else {
      names(turn.mat)=genes.use
      return(turn.mat)
    }
  }
}

genePlot.pseudo <- function (object,gene,pseudo="pseudo.time", cell.ids = NULL, col.use = NULL, 
                             pch.use = 16, cex.use = 1.5, use.imputed = FALSE, do.ident = FALSE, 
                             do.spline = FALSE,do.line=FALSE,spline.span = 0.75,name.x=NULL, do.logit=FALSE,
                             do.smooth=FALSE,knn=10,use.scale=FALSE,do.return=FALSE,conf.int=FALSE,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  nn.count.add <- function(mat, k, dmat) {
    dmat[dmat == 0] <- Inf
    knn <- t(apply(dmat, 1, function(x) rank(x) <= k))
    mat <- mat + mat %*% t(knn)
    return(mat)
  }
  
  if (do.smooth){
    pt.ncc=fetch.data(object,pseudo)
    dmat = as.matrix(dist(pt.ncc))
    count.mat = object@scale.data
    object@scale.data = nn.count.add(count.mat, knn, dmat)
    data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }else{
    data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = use.scale)))
  }
  
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  name.x = set.ifnull(name.x, pseudo)
  pseudo.order=order(data.use[pseudo,cell.ids])
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = as.numeric(data.use[gene, cell.ids])
  ident.use = as.factor(object@ident[cell.ids])
  pseudo.data=as.matrix(fetch.data(object,pseudo))
  
  if(do.logit){
    if(do.smooth){
      gene.norm=apply(matrix(object@scale.data[gene,colnames(object@scale.data)]),2, FUN = function(X) (X - min(X))/diff(range(X)))
    }
    
    if(do.spline){
      loess.fit = loess(g2 ~ g1, span = spline.span)
      gene.norm = (loess.fit$fitted - min(loess.fit$fitted))/diff(range(loess.fit$fitted))
    }
    
    if(!do.smooth & !do.spline){
      gene.norm=apply(object@data[gene,colnames(object@data)],1, FUN = function(X) (X - min(X))/diff(range(X)))
    }
    
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm[pseudo.order]~g1,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn),1]
  }
  
  if (length(col.use) > 1) {
    col.use = col.use[as.numeric(ident.use)]
  }
  else {
    col.use = set.ifnull(col.use, as.numeric(ident.use))
  }
  
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = unlist(object@data[gene, cell.ids])
  gene.cor = round(cor(g1, g2), 2)
  plot(g1, g2, xlab = name.x, ylab = gene, col = col.use, cex = cex.use, 
       main = "", pch = pch.use, font.lab=2,cex.lab=1.2, ...)
  
  if (do.logit){
    abline(v=turn.point,lwd=4,col="purple",lty="longdash")
  }
  
  if (do.spline) {
    loess.fit = loess(g2 ~ g1, span = spline.span)
    lines(g1[order(pseudo.data)], loess.fit$fitted[order(pseudo.data)],col = "black",lwd=4)
    if(conf.int){
      prid = predict(loess.fit,se=T)
      lines(g1[order(g1)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(g1)],col = "black",lwd=1.2,lty=2)
      lines(g1[order(g1)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(g1)],col = "black",lwd=1.2,lty=2)
    }
  }
  
  if(do.line){
    lines(g1[order(pseudo.data)], g2[order(pseudo.data)],col = "purple",lwd=3)
  }
  
  if (do.ident) {
    return(identify(g1, g2, labels = cell.ids))
  }
  
  if (do.return){
    if (do.spline){
      return(data.frame(g1,loess.fit$fitted))}else{
        return(data.frame(g1,g2))
      }
  }
}


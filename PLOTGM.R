# GM Production Code
# Filename: PlotGM.R
# Author: CA

# Libraries
library(shape)

# Plot graph function
plotgraph <- function(result, G=NULL, altnames, ynames,Ycols="white", Ytextcols=0, Xcols="white", Xtextcols=0) {
  
  # Code for drawing chain graphs
  zigzag <- function(n) {
    z <- (rep(c(-1,1), n%/%2+1) * rep(1:(n%/%2+1), each=2))[1:n]
    -z
  }
  
  vars <- colnames(result)[!grepl(":",colnames(result))]
  ycoords <- data.frame(x=zigzag(nrow(result))/2, y=(nrow(result):1-(nrow(result)+1)/2)*max(1,length(vars)/nrow(result)))
  xcoords <- data.frame(x=min(ycoords$x)-3, y=(length(vars):1-((length(vars)+1)/2))*max(1,nrow(result)/length(vars)))
  
  par(mar=rep(0,4))
  plot(NULL, xlim=range(c(xcoords$x, ycoords$x))+1/c(-2,2),  ylim=range(c(xcoords$y, ycoords$y))+1/c(-2,2), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  if(!is.null(G)) {
    for (i in 2:nrow(result))
      for (j in 1:(i-1))
        if (G[i,j]==1)
          lines(ycoords$x[c(i,j)], ycoords$y[c(i,j)], lwd=2)
  }

  
  # Get rid of effects implied by higher-order interactions
  gresult <- result
  for (i in 1:(ncol(gresult)-1))
    for (j in 1:nrow(gresult))
      if (gresult[j,i]) {       
        cvars1 <- strsplit(colnames(gresult)[i],":")[[1]]
        for (k in (i+1):ncol(gresult))
          if (gresult[j,k]) {                    
            cvars2 <- strsplit(colnames(gresult)[k],":")[[1]]
            if (all(cvars1 %in% cvars2))
              gresult[j,i] <- FALSE
          }
      }
  
  
  connections <- list()
  for (i in 1:ncol(result))
    for (j in 1:nrow(result))
      if (gresult[j,i])             
        connections <- c(connections, list(list(from=which(vars %in% strsplit(colnames(result)[i],":")[[1]]), to=j, label=colnames(result)[i])))
  
  if(length(connections)==0)  {
    
    points(ycoords$x, ycoords$y, cex=10, pch=21, bg=Ycols)
    text(ycoords$x, ycoords$y, ynames, cex=0.75, col=Ytextcols)
    
    points(xcoords$x, xcoords$y, cex=10, pch=21, bg=Xcols)
    text(xcoords$x, xcoords$y, altnames, cex=0.75, col=Xtextcols)
    
    if (!missing(altnames))
      vars <- altnames[vars]
    text(xcoords$x, xcoords$y, vars, cex=0.75, col=Xtextcols)
  }else{
    delta <- 0.15
    
    fromslotcount <- integer(length(vars))
    fromslotused <- integer(length(vars))
    fromslotcoords <- list()
    for (i in 1:length(vars)) {
      fromslotcount[i] <- sum(sapply(connections,function(a) i %in% a$from))
      if (fromslotcount[i]==1) {
        fromslotcoords[[i]] <- xcoords$y[i]
      } else {
        fromslotcoords[[i]] <- xcoords$y[i] + seq(-delta, delta, len=fromslotcount[i])
      }
    }
    
    for (k in order(ycoords$y[sapply(connections, function(a) min(a$to))])) {
      connections[[k]]$fromy <- c()
      for (l in connections[[k]]$from) {
        fromslotused[l] <- fromslotused[l]+1
        connections[[k]]$fromy <- c(connections[[k]]$fromy, fromslotcoords[[l]][fromslotused[l]])
      }
    }
    
    
    toslotcount <- integer(nrow(result))
    toslotused <- integer(nrow(result))
    toslotcoords <- list()
    for (i in 1:nrow(result)) {
      toslotcount[i] <- sum(sapply(connections,function(a) i %in% a$to))
      if (toslotcount[i]==1) {
        toslotcoords[[i]] <- ycoords$y[i] 
      } else {
        toslotcoords[[i]] <- ycoords$y[i] + seq(-delta, delta, len=toslotcount[i])
      }
    }
    
    for (k in order(xcoords$y[sapply(connections, function(a) min(a$from))])) {
      connections[[k]]$toy <- c()
      for (l in connections[[k]]$to) {
        toslotused[l] <- toslotused[l]+1
        connections[[k]]$toy <- c(connections[[k]]$toy, toslotcoords[[l]][toslotused[l]])
      }
    }
    
    lspace <- 0.25
    rspace <- 0.5
    for (i in 1:length(connections)) {
      for (j in 1:length(connections[[i]]$fromy))
        lines(c(xcoords$x[connections[[i]]$from[j]]+lspace, ycoords$x[connections[[i]]$to]-2*rspace), c(connections[[i]]$fromy[j], connections[[i]]$toy), lwd=2)
      Arrows(ycoords$x[connections[[i]]$to]-2*rspace, connections[[i]]$toy, ycoords$x[connections[[i]]$to]-rspace, connections[[i]]$toy, lwd=2)
      if (length(connections[[i]]$fromy)>1)
        points(ycoords$x[connections[[i]]$to]-2*rspace, connections[[i]]$toy, pch=16)
    }
    
    points(ycoords$x, ycoords$y, cex=10, pch=21, bg=Ycols)
    text(ycoords$x, ycoords$y, ynames, cex=0.75, col=Ytextcols)
    
    points(xcoords$x, xcoords$y, cex=10, pch=21, bg=Xcols)
    text(xcoords$x, xcoords$y, altnames, cex=0.75, col=Xtextcols)
    
    if (!missing(altnames))
      vars <- altnames[vars]
    text(xcoords$x, xcoords$y, vars, cex=0.75, col=Xtextcols)
  }
  
}
# Functionality to plot grpahs from given model output
# Filename: PlotGraph.R
# Author: CA

# Source files
source("PlotGM.R")
source("ChordalGraph.R")

# Generic plot graph function, using PlotGM
PlotGraph <- function(x,select.model=1,mult.models=FALSE)  {
  
  # Obtain rankings of models
  rnd <- rnorm(length(x$actives[1,,]))
  code <- numeric(length(x$beta))
  for (i in 1:(length(x$beta)))
    code[i] <- sum(rnd*x$actives[i,,])+x$Graphs[i]
  t <- sort(table(code), decreasing=TRUE)
  
  altnames <- x$fixed
  Xcols <- rep("pink",length(altnames))
  
  if(ncol(x$Y)==1)  {
    Ycols <- c(rgb(237/255,28/255,36/255))
    
    if(mult.models==TRUE) {
      par(mfrow=c(2,2))
      for(j in 1:4) {
        id <- which.min(abs(code-as.numeric(names(t)[j])))
        result <- (x$actives[id,,])
        if(is.null(ncol(result))) {
          result <- matrix(result,nrow=1)
        }
        colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
        plotgraph(result, G=NULL, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
        text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[j]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
      }
    }else{
      par(mfrow=c(1,1))
      id <- which.min(abs(code-as.numeric(names(t)[select.model])))
      result <- (x$actives[id,,])
      if(is.null(ncol(result))) {
        result <- matrix(result,nrow=1)
      }
      colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
      plotgraph(result, G=NULL, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
      text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[select.model]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
    }
 }
  
  if(ncol(x$Y)==2)  {
    Ycols <- c(rgb(237/255,28/255,36/255),
               rgb(0/255,0/162,232/255))
    
    if(mult.models==TRUE) {
      par(mfrow=c(2,2))
      for(j in 1:4) {
        id <- which.min(abs(code-as.numeric(names(t)[j])))
        result <- (x$actives[id,,])
        if(is.null(ncol(result))) {
          result <- matrix(result,nrow=1)
        }
        colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
        G <- graphs2[[x$Graphs[id]]]
        plotgraph(result, G=G, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
        text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[j]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
      }
    }else{
      id <- which.min(abs(code-as.numeric(names(t)[select.model])))
      result <- (x$actives[id,,])
      if(is.null(ncol(result))) {
        result <- matrix(result,nrow=1)
      }
      
      colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
      G <- graphs2[[x$Graphs[id]]]
      
      plotgraph(result, G, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
      text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[select.model]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
    }
  }
  
  if(ncol(x$Y)==3)  {
    Ycols <- c(rgb(237/255,28/255,36/255),
               rgb(0/255,0/162,232/255),
               rgb(34/255,177/255, 76/255))

    if(mult.models==TRUE) {
      par(mfrow=c(2,2))
      for(j in 1:4) {
        id <- which.min(abs(code-as.numeric(names(t)[j])))
        result <- (x$actives[id,,])
        if(is.null(ncol(result))) {
          result <- matrix(result,nrow=1)
        }
        colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
        G <- graphs3[[x$Graphs[id]]]
        plotgraph(result, G=G, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
        text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[j]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
      }
    }else{
      par(mfrow=c(1,1))
      id <- which.min(abs(code-as.numeric(names(t)[select.model])))
      result <- (x$actives[id,,])
      if(is.null(ncol(result))) {
        result <- matrix(result,nrow=1)
      }
      
      colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
      G <- graphs3[[x$Graphs[id]]]
      
      plotgraph(result, G, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
      text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[select.model]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
    }
  }
  
  if(ncol(x$Y)>3) {
    Ycols <- rep("red",ncol(x$Y))

    if(mult.models==TRUE) {
      par(mfrow=c(2,2))
      for(j in 1:4) {
        id <- which.min(abs(code-as.numeric(names(t)[j])))
        result <- (x$actives[id,,])
        if(is.null(ncol(result))) {
          result <- matrix(result,nrow=1)
        }
        colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
        plotgraph(result, G=NULL, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
        text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[j]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
      }
    }else{
      par(mfrow=c(1,1))
      id <- which.min(abs(code-as.numeric(names(t)[select.model])))
      result <- (x$actives[id,,])
      if(is.null(ncol(result))) {
        result <- matrix(result,nrow=1)
      }
      
      colnames(result) <- sapply(x$cmodel$blocks,"[[","label")
      
      plotgraph(result, G=NULL, altnames, x$resp, Ycols, Xcols, Ytextcols="white", Xtextcols = "black")
      text(par()$usr[2], par()$usr[3], paste0("Posterior probability of model: ", round(1000*mean(abs(code-as.numeric(names(t)[select.model]))<1e-8))/10,"%"), col="black", adj=c(1,-0.1))
    }

  }
}

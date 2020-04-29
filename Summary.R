##Summary table function for CG Model

# Libraries
library(coda)

Summarise<-function(x){
  
  ## Format the arguments    
  n <- length(x$beta)  # Number of samples
  n.data <- nrow(x$Y)  # Number of observations
  n.resp <- length(x$beta[[1]]) # Number of responses
  p <- length(x$beta[[1]][[1]]) # number of terms
  q <- length(x$rLambda[[1]])   # number of random effects
  n.refs <- c()
  for(j in 1:q)
    n.refs[j] <- nrow(x$reffs[[1]][[j]])
  

  ##BETA##

  ##Calculate summary quantities of interest
  
  # Rescale beta to matrix
  betamat <- matrix(unlist(x$beta),ncol=n.resp*p,byrow=T)
  #Posterior medians
  #Beta
  beta.median <- apply(betamat,2,median)  
  
  #95% credible intervals
  #Beta
  beta.interval<-apply(betamat,2,quantile,probs=c(0.025,0.975))  
  
  #Effective sample size
  #Beta
  beta.ess<-apply(mcmc(betamat),2,effectiveSize)
  
  
  #Geweke diagnostic
  #Beta
  geweke.beta<-apply(mcmc(betamat),2,geweke.diag)
  

  ## MODEL ERROR
  
  # Obtain variance estimates from list 
  model.varmat <- matrix(1/unlist(lapply((x$Lambda.eps),diag)),ncol=n.resp,byrow=T)
  
  # Posterior medians
  model.varmedian <- apply(model.varmat,2,median)
  
  #95% credible intervals
  model.varinterval<-apply(model.varmat,2,quantile,probs=c(0.025,0.975))  
  
  #Effective sample size
  model.varess<-apply(mcmc(model.varmat),2,effectiveSize)
  
  ## RANDOM EFFECT ERROR
  ref.varmat <- vector(length=q, mode='list')
  for(j in 1:q){
    ref.varmat[[j]] <- matrix(0,nrow=n,ncol=n.resp)
    for(k in 1:n){
      ref.varmat[[j]][k,] <- (1/(diag((x$rLambda[[k]][[j]]))))
    }
  }
    
    
  
  ref.varmedian <- vector(length=q, mode='list')
  ref.varinterval <- vector(length=q, mode='list')
  ref.varess <- vector(length=q, mode='list')
  
  for(j in 1:q){
    # Posterior medians
    ref.varmedian[[j]] <- apply(ref.varmat[[j]],2,median)
    
    #95% credible intervals
    ref.varinterval[[j]]<-apply(ref.varmat[[j]],2,quantile,probs=c(0.05,0.95))  
    
    #Effective sample size
    ref.varess[[j]] <-apply(mcmc(ref.varmat[[j]]),2,effectiveSize)
  }
  
  

  ##Beta Summary
  median<-c(beta.median)
  
  #add sigma2.interval as a new column in beta.interval
  interval<-cbind(beta.interval)
  
  # Add in ESS for each coefficient
  ess <- c(beta.ess)
  
  # Add in Geweke statistic for each coeff
  #access z only from geweke diagnostic for beta and save as number
  z.score<-sapply(geweke.beta, "[[", 1)
  z.beta<-as.numeric(z.score)
  
 
  # Create summary for fixed effects
  summary<-data.frame(rbind(median, interval, ess,z.beta))
  
  #rename column names for posterior summary matrix
  #beta
  varnames <- list()
  for(i in 1:n.resp)
    varnames[[i]] <- paste(x$resp[i],"-",colnames(x$cmodel$X))
      
  varnames <- unlist(varnames)
  for(i in 1:(p*n.resp)){
    colnames(summary)[i]<-varnames[i]
  }
  rownames(summary) <- c("post.median","lower-95% CI","upper-95% CI","ESS","Geweke-R")

  ## Create summary for variance parameters
  model.varsummary <- data.frame(rbind(model.varmedian,model.varinterval,model.varess))
  colnames(model.varsummary) <- model.1$resp
  rownames(model.varsummary) <- c("post.median","lower-95% CI","upper-95% CI","ESS")
  
  ## Create summary for random effects variance parameters
  ref.varsummary <- vector(length=q,mode='list')
  
  for(j in 1:q) {
    ref.varsummary[[j]] <- suppressWarnings(data.frame(rbind(ref.varmedian[[j]],(ref.varinterval[[j]]),c(ref.varess[[j]]))))
    colnames(ref.varsummary[[j]]) <- model.1$resp
    rownames(ref.varsummary[[j]]) <- c("post.median","lower-95% CI","upper-95% CI","ESS")
   }

  
  
  #Print the number of samples and summary info
  
  cat("\n")
  cat("Number of iterations:", n )
  cat("\n")
  cat("Model Formula:")
  print(x$model)
  cat("\n")
  cat("Number of observations:",n.data)
  cat("\n")
  cat("Number of response variables:",n.resp)
  cat("\n")
  cat("\n")
  
  ## Print the posterior summary for the random effects variance
  for(j in 1:q) {
    cat("Posterior Summary of", x$rand[j])
    cat("\n")
    print(t(ref.varsummary[[j]]))
    cat("Number of levels for",x$rand[j],"-",n.refs[j])
    cat("\n")
    cat("\n")
  }
  
  ## Print the posterior summary for the model error
  cat("Posterior Summary of Residuals")
  cat("\n")
  print(t(model.varsummary))
  
  cat("\n")
  ##Print the following posterior summary Beta
  cat("Posterior Summary of Fixed Effects:")
  cat("\n")
  print(t(summary))
  cat("\n")
}

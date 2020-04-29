# Full sampler code with data sorting and simple functionality
# Filename - SimplifiedSampler.R
# Author - CA

# Libraries
library(MASS)
library(mvtnorm)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(foreach)
library(doParallel)
library(MatrixModels)
library(MCMCpack)
registerDoParallel(cores=detectCores())

# Source files
sourceCpp("GWishGen.cpp")
source("MixModelFunctions.R")
source("ChordalGraph.R")
source("DoubleMHMultReffs.R")

##       ##
# Sampler #
##       ##

Simplified_sampler <- function(data,model,random,iterations=500,interactions=1,model.selection=TRUE,Graphs=TRUE,
                        tau,a0.tau=NULL,b0.tau=NULL,d0.reffs=NULL,V0.reffs=NULL,d0.eps=NULL,V0.eps=NULL,G=NULL,
                        nesting=TRUE,nest.terms=NULL,nest.reff=NULL,
                        param.ex=TRUE,phi.alpha=NULL,phi.beta=NULL,
                        parallel.blocks=detectCores()) {
  
  # Split formula into response, fixed and random effects
  split.model  <- strsplit(Reduce(paste,deparse(model)),"~" )[[1]]
  feffects     <- trimws(strsplit(split.model[2],"\\+")[[1]])
  response     <- trimws(strsplit(split.model[1],"\\+")[[1]])
  split.random <- strsplit(Reduce(paste,deparse(random)),"~")[[1]]
  reffects     <- trimws(strsplit(split.random[2],"\\+")[[1]])
  
  # Creation of input for data
  cmodel <- prepare.data(data[,feffects], max.interaction=interactions)
  
  Y <- as.matrix(data[, response],drop=FALSE)
  
  rdata <- data[,reffects]
  
  rlookup <- list()
  rrlookup <- list()
  rcounts <- list()
  if(is.null(ncol(rdata))==TRUE)  {
    vals <- unique(rdata)
    rlookup[[1]] <- sapply(vals, function(a) which(rdata==a),simplify=FALSE)
    names(rlookup[[1]]) <- vals
    rcounts[[1]] <- sapply(rlookup[[1]],length)
    rrlookup[[1]] <- integer(nrow(Y))
    for(k in 1:length(vals))  {
      rrlookup[[1]][rlookup[[1]][[k]]] <- k
    }
  }else{
    for (j in 1:ncol(rdata)) {
      vals <- unique(rdata[,j])
      rlookup[[j]] <- sapply(vals, function(a) which(rdata[,j]==a))
      names(rlookup[[j]]) <- vals
      rcounts[[j]] <- sapply(rlookup[[j]], length)
      rrlookup[[j]] <- integer(nrow(Y))
      for (k in 1:length(vals)) {
        rrlookup[[j]][rlookup[[j]][[k]]] <- k
      }
    }
  }
  
  
  
  main <- non.interaction.terms(cmodel)
  active <- list()
  for (i in 1:ncol(Y))
    active[[i]] <- main
  
  X <- lapply(active, getX, cmodel=cmodel)
  
  # Generate priors if not specified
  if(missing(tau))
    tau <- 1 
  if(length(tau)==1)
    tau <- rep(tau,length(X))
  if(is.null(a0.tau)==TRUE)
    a0.tau <- 1e-4
  if(is.null(b0.tau)==TRUE)
    b0.tau <- 1e-4
  if(is.null(d0.eps)==TRUE)
    d0.eps <- ncol(Y)
  if(is.null(V0.eps)==TRUE)
    V0.eps <- diag(ncol(Y))
  if(is.null(d0.reffs)==TRUE)
    d0.reffs <- rep(ncol(Y),length(rlookup))
  if(is.null(V0.reffs)==TRUE) {
    V0.reffs <- list()
    for(j in 1:length(rlookup))
      V0.reffs[[j]] <- diag(ncol(Y))
  }
  if(param.ex==TRUE){
    if(is.null(phi.alpha)==TRUE)
      phi.alpha <- 10
    if(is.null(phi.beta)==TRUE)
      phi.beta <- 10
  }
  if(is.null(G)==TRUE)  
    G <- diag(0,ncol(Y))
  
  # Generate initial estimates for parameters
  beta <- list()
  for (j in 1:ncol(Y))
    beta[[j]] <- numeric(ncol(X[[j]]))
  Lambda.eps <- diag(ncol(Y))
  Yhat <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
  Ywork <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
  
  
  reffs <- lapply(rlookup, function(a) { result <- matrix(0, nrow=length(a), ncol=ncol(Y)); rownames(result) <- names(a); result})
  
  rLambda <- rep(list(diag(ncol(Y))), length(reffs))
  
  if(nesting==T)  {
    if(is.null(nest.terms)==TRUE) {
      nest.terms <- list()
      nest.reff <- vector()
      for(j in 1:length(reffs)) {
        counts <- as.data.frame(data %>% group_by(get(reffects[j]))  %>% summarise_all(n_distinct) %>% summarise_if(max, .predicate=is.numeric))
        vars <- names(counts)[counts==1]
        if(length(vars)==0)  {
          nest.terms[[j]] <- NULL
        }else{
          nest.terms[[j]] <- intersect(vars,feffects)
          nest.reff <- c(nest.reff,j)
        }
      }
    }
    if(length(nest.terms)==0)
      nesting <- FALSE
  }
  
  # Storage objects
  betas <- vector(length=iterations, mode='list')
  Lambda.epss <- vector(length=iterations,mode='list')
  reffss <- vector(length=iterations, mode='list')
  rLambdas <- vector(length=iterations, mode='list')
  actives <- array(FALSE, dim=c(iterations, ncol(Y), length(cmodel$blocks)))
  Gs <- integer(iterations)
  acceptance.rate <- vector(length=iterations,mode='list')

  # Initialise progress bar
  pb <- txtProgressBar(0,iterations)
  
  
  for(ix in 1:iterations) {
    setTxtProgressBar(pb,ix) # Set state in progress bar
    
    # Compute current working Y for beta calculation by removing all random effects
    Ywork <- Y
    for(j in 1:length(reffs)) {
      Ywork <- Ywork - reffs[[j]][rrlookup[[j]],]
    }
    
    # Update all betas at once with new parameters from previous iteration
    
    
    # Sample beta update
    # Generate tau prior on beta covariance
    for(j in 1:length(beta)) {
      if(ix!=1) # Don't update tau on the first iteration
        tau[j] <- 1/(rgamma(1,a0.tau + (length(beta[[j]])/2), b0.tau + sum(beta[[j]]^2)/2))
      
      # Update beta
      result <- compute.posterior.beta(Ywork, X, X[[j]], beta, j, tau, Lambda.eps)
      beta[[j]] <- simulate.posterior.beta(result)
      Yhat[,j] <- drop(Ywork[,j]-X[[j]]%*%beta[[j]])
      
      
      # Joint update on beta with adding/removing candidate terms
      if(model.selection==TRUE) {
        
        # Try to add a covariate
        cands <- candidate.add.terms(active[[j]], cmodel$blocks)
        if (length(cands)>0) {
          current.lmodelev <- compute.lmodelev(result,Y,Yhat, X[[j]], j, tau, Lambda.eps)
          k <- sample(length(cands),1)
          new.active <- c(active[[j]], cands[k])
          new.X <- getX(cmodel, new.active)
          new.result <- compute.posterior.beta(Ywork, X, new.X, beta, j, tau, Lambda.eps)
          new.lmodelev <- compute.lmodelev(new.result,Y, Yhat, new.X, j, tau, Lambda.eps)
          backs <- length(candidate.drop.terms(new.active, cmodel$blocks))
          logacc <- new.lmodelev - current.lmodelev + log(length(cands)) - backs
          if (log(runif(1))<=logacc) {
            active[[j]] <- new.active
            result <- new.result
            X[[j]] <- new.X
            beta[[j]] <- simulate.posterior.beta(new.result)
          }
        }
        
        # Try to remove covariate
        cands <- candidate.drop.terms(active[[j]], cmodel$blocks)
        if (length(cands)>0) {
          current.lmodelev <- compute.lmodelev(result,Y,Yhat, X[[j]], j, tau, Lambda.eps)
          k <- sample(length(cands),1)
          new.active <- setdiff(active[[j]], cands[k])
          new.X <- getX(cmodel, new.active)
          new.result <- compute.posterior.beta(Ywork, X, new.X, beta, j, tau, Lambda.eps)
          new.lmodelev <- compute.lmodelev(new.result,Y,Yhat, new.X, j, tau, Lambda.eps)
          backs <- length(candidate.add.terms(new.active, cmodel$blocks))
          logacc <- new.lmodelev - current.lmodelev + log(length(cands)) - backs
          if (log(runif(1))<=logacc) {
            active[[j]] <- new.active
            result <- new.result
            X[[j]] <- new.X
            beta[[j]] <- simulate.posterior.beta(new.result)
          }
        }
      }
    }

    # Compute random effects coefficients
    for (k in 1:ncol(Y))
      Ywork[,k] <- drop(Y[,k]-X[[k]]%*%beta[[k]]) # Obtain working Y by removing betas
    
    for (j in 1:length(reffs)) {
      Ywork2 <- Ywork
      for (l in 1:length(reffs))
        if (j!=l)
          Ywork2 <- Ywork2 - reffs[[l]][rrlookup[[l]],]  # Remove all other random effects for working Y 
        
        
        # For parralelisation, split counts and lookup tables
        rlevel <- nrow(reffs[[j]])
        rcount.list <- list()
        rlookup.list <- list()
        for(k in 1:parallel.blocks){
          if(k==parallel.blocks)  {
          rcount.list[[k]] <- rcounts[[j]][(((k-1)*(rlevel%/%parallel.blocks))+1):(k*(rlevel%/%parallel.blocks) + rlevel%%parallel.blocks)]
          rlookup.list[[k]] <- rlookup[[j]][(((k-1)*(rlevel%/%parallel.blocks))+1):(k*(rlevel%/%parallel.blocks)+ rlevel%%parallel.blocks)]
        }else{
          rcount.list[[k]] <- rcounts[[j]][(((k-1)*(rlevel%/%parallel.blocks))+1):(k*(rlevel%/%parallel.blocks))]
          rlookup.list[[k]] <- rlookup[[j]][(((k-1)*(rlevel%/%parallel.blocks))+1):(k*(rlevel%/%parallel.blocks))]
          }
        }
         reffs[[j]] <- (foreach(k=1:parallel.blocks,.combine = "rbind",.export = "getReffs") %dopar%
                          getReffs(rLambda[[j]],Lambda.eps,rcount.list[[k]],rlookup.list[[k]],Ywork2)
        )
    }
    
    ## Nesting step ##
    if(nesting==T)  {
      for(j in 1:length(nest.reff)) {
        if(length(nest.terms[[j]])!=0)  {
          # Generate nested design matrix
          Nest.res <- getNestedX(nest.terms[[j]],active,rlookup[[nest.reff[j]]],cmodel)
          X.alpha <- Nest.res$Xnest
          Xalpha <- blockdiag(X.alpha)
          cols <- Nest.res$columns
          beta.nest <- list()
          betaexp <- expand.beta(beta,active,cmodel)
          
          # Form alpha
          for(k in 1:ncol(Y)) {
            beta.nest[[k]] <- betaexp[[k]][cols[[k]]] 
          }
          
          alpha <- Xalpha%*%unlist(beta.nest) + c(reffs[[nest.reff[j]]])
          alphamat <- matrix(alpha,ncol=ncol(Y))
          
          betanest <- list()
          for(k in 1:ncol(Y)) {
           result <- compute.posterior.beta(alphamat, X.alpha, X.alpha[[k]], beta.nest, k, tau, rLambda[[nest.reff[j]]])
           betanest[[k]] <- simulate.posterior.beta(result)
          }

          
          # Update beta
          for(k in 1:ncol(Y)) {
            betaexp[[k]][cols[[k]]] <- betanest[[k]]
          }
          
          for(k in 1:ncol(Y)) {
            beta[[k]] <- betaexp[[k]][betaexp[[k]]!=0]
          }
          
          # Updata corresponding random effect
          reffs[[nest.reff[j]]] <- matrix(as.vector(alpha) - Xalpha%*%unlist(betanest),ncol=ncol(Y))
        }
      }
    }
    
    if(Graphs==TRUE)  {
      # Precision estimates using GM step
      if(ncol(Y)==1)  {
        for (j in 1:length(reffs)) 
          Ywork <- Ywork - reffs[[j]][rrlookup[[j]],] 
        
        Lambda.eps <- as.matrix(rwish(v=nrow(Y)+d0.eps,S=solve(crossprod(Ywork)+V0.eps)))
        
        for(j in 1:length(reffs)) {
          rLambda[[j]] <- as.matrix(rwish(v=nrow(reffs[[j]])+d0.reffs[j],S=solve(crossprod(reffs[[j]])+V0.reffs[[j]])))
        }
      }
      
      if(ncol(Y)==2)  {
        for(j in 1:length(reffs))
          Ywork <- Ywork - reffs[[j]][rrlookup[[j]],] 
        
        evs <- numeric(length(graphs2))
        for (k in 1:length(graphs2)) {
          evs[k] <- graphlgev2(Ywork, graphs2[[k]], d0.eps, V0.eps)
          for (j in 1:length(reffs)) 
            evs[k] <- evs[k]+graphlgev2(reffs[[j]], graphs2[[k]], d0.reffs[j], V0.reffs[[j]])
        }
        evs <- evs - max(evs)
        Gs[[ix]] <- sample(length(graphs2), 1, prob=exp(evs))
        
        G <- graphs2[[Gs[[ix]]]]
        Lambda.eps <- rGWish_sampler(1, nrow(Y)+d0.eps,  solve(crossprod(Ywork)+V0.eps),G,100)[,,1]
        for (j in 1:length(reffs)) {
          rLambda[[j]] <- rGWish_sampler(1, d0.reffs[j]+nrow(reffs[[j]]), solve(crossprod(reffs[[j]])+V0.reffs[[j]]),G,100)[,,1]
        }
      }
      
      if(ncol(Y)==3)  {
        for (j in 1:length(reffs)) 
          Ywork <- Ywork - reffs[[j]][rrlookup[[j]],]   
        
        evs <- numeric(length(graphs3))
        for (k in 1:length(graphs3)) {
          evs[k] <- graphlgev3(Ywork, graphs3[[k]], d0.eps, V0.eps)
          for (j in 1:length(reffs)) 
            evs[k] <- evs[k]+graphlgev3(reffs[[j]], graphs3[[k]], d0.reffs[j], V0.reffs[[j]])
        }
        evs <- evs - max(evs)
        Gs[[ix]] <- sample(length(graphs3), 1, prob=exp(evs))
        
        G <- graphs3[[Gs[[ix]]]]
        
        Lambda.eps <- rGWish_sampler(1, nrow(Y)+d0.eps,  solve(crossprod(Ywork)+V0.eps),G,100)[,,1]
        for (j in 1:length(reffs)) {
          rLambda[[j]] <- rGWish_sampler(1, d0.reffs[j]+nrow(reffs[[j]]), solve(crossprod(reffs[[j]])+V0.reffs[[j]]),G,100)[,,1]
        }
      }
      if(ncol(Y)>3){
        Ywork <- Ywork - reffs[[j]][rrlookup[[j]],]   
        
        prec.res <- double.MH(Ywork,reffs,iterations=1,b.eps=d0.eps,D.eps=V0.eps,b.reffs=d0.reffs,D.reffs=V0.reffs,G=G)
        
        Lambda.eps <- prec.res$Precision.eps[[1]]
        for(j in 1:length(reffs)) {
          rLambda[[j]] <- prec.res$Precision.reffs[[j]][[1]]
        }
        G <- prec.res$G[[1]]
      }
    }else{
      # Precision estimates using standard Wishart prior
      for (j in 1:length(reffs)) 
        Ywork <- Ywork - reffs[[j]][rrlookup[[j]],]   
      
      Lambda.eps <- as.matrix(rwish(v=nrow(Y)+d0.eps,S=solve(crossprod(Ywork)+V0.eps)))
      for (j in 1:length(reffs)) {
        rLambda[[j]] <-  as.matrix(rwish(v=nrow(reffs[[j]])+d0.reffs[j],S=solve(crossprod(reffs[[j]])+V0.reffs[[j]])))
      }
    }
    
    
    ## Parameter expansion step ##
    if(param.ex==T) {
      for(j in 1:length(reffs)) {
        # Draw phi from the inverse gamma distribution
        phi <- rgamma(1,phi.alpha,phi.beta)
        
        # Form effect.star and precision.star
        effect.star <- reffs[[j]]*phi
        prec.star <- rLambda[[j]]/(phi^2)
        
        # Compute acceptance rate
        original.rate <- dgamma(phi,phi.alpha,phi.beta,log=TRUE)+sum(dmvnorm(Ywork,rep(0,ncol(Y)),solve(Lambda.eps),log=TRUE))+
          sum(dmvnorm(reffs[[j]],rep(0,ncol(Y)),solve(rLambda[[j]]),log=TRUE))+
          log(dwish(rLambda[[j]],d0.reffs[j],V0.reffs[[j]]))
        
        Ystar <- NULL
        Ystar <- Ywork + reffs[[j]][rrlookup[[j]],]
        Ystar <- Ystar - effect.star[rrlookup[[j]],]
        
        star.rate <- dgamma(1/phi,phi.alpha,phi.beta,log=TRUE)+sum(dmvnorm(Ystar,rep(0,ncol(Y)),solve(Lambda.eps),log=TRUE))+ 
          sum(dmvnorm(effect.star,rep(0,ncol(Y)),solve(prec.star),log=TRUE))+
          log(dwish(prec.star,d0.reffs[j],V0.reffs[[j]]))
        
        factors <- ((length(reffs[[j]])+(ncol(Y)*(ncol(Y)+1))))*log(phi)
        
        acc.rate <- star.rate - original.rate + factors
        
        if(log(runif(1)) < acc.rate) {
          reffs[[j]] <- effect.star
          rLambda[[j]] <- prec.star
          acceptance.rate[[j]][ix] <- 1
        }else{
          acceptance.rate[[j]][ix] <- 0
        }
      }
    }
    
    reffss[[ix]] <- reffs
    rLambdas[[ix]] <- rLambda
    
    betas[[ix]] <- expand.beta(beta, active, cmodel)
    Lambda.epss[[ix]] <- Lambda.eps
    
    for (j in 1:length(active))
      actives[ix,j,active[[j]]] <- TRUE
  }
  return(list(beta=betas,reffs=reffss,Lambda.eps=Lambda.epss,rLambda=rLambdas,Graphs=Gs,actives=actives,acceptance=acceptance.rate,cmodel=cmodel,resp=response,fixed=feffects,rand=reffects,model=model,Y=Y))
}
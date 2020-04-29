# Functionality used within multivariate mixed model sampler
# Filename : MixModelFunctions.R
# Author: CA

##                                          ##
# Functions required in SimplifiedSampler.R  #
##                                          ##


# Functionality to prepare data matrix
prepare.data <- function(data, columns=colnames(data), max.interaction=1) {
  require(MatrixModels)
  formula <- as.formula(paste0("~",paste0(rep(".",max.interaction), collapse="*"))) # Initialise (max interaction sets interaction levels)
  term.labels <- attr(terms(formula, data=data[,columns]),"term.labels")   # Obtain term names (sets all possible combinations of interactions too)
  term.components <- lapply(term.labels,function(a) strsplit(a,":")[[1]])  # Put all possible terms in list object
  detailed.formula <- as.formula(paste0("~", paste0(term.labels, collapse="+"))) # Create full model from list object and adding terms
  X <- model.Matrix(detailed.formula, data=data, sparse=TRUE) # Generate model matrix, but from MatrixModels package exploiting sparsity
  X.blocks <- list()   # Generate block object
  for (i in unique(X@assign[-1])) {  # For i in (list the number of columns, dropping the intercept)
    # Here, we assign the information regarding each column in the data matrix in terms of column, the term label and in the case of interactions, which "parents" they arise from (their interacting term)
    info <- list()
    info$columns <-which(X@assign==i) # Assign column number
    info$label <- term.labels[i]      # Assign name corresponding to column number
    # Parents are one interaction level lower (parent level 0 at the non-interaction case)
    info$parents <- which((sapply(term.components, length)==length(term.components[[i]])-1) & sapply(term.components, function(a) all(a%in%term.components[[i]])))
    X.blocks[[i]] <- info
  }
  for (i in unique(X@assign[-1])) 
    X.blocks[[i]]$children <- which(sapply(X.blocks, function(a) i%in%a$parents)) # Here, we assign the "children" for the non interaction terms, indicating which blocks the term is present in (which interaction terms it's contained within)
  list(X=X, blocks=X.blocks)  # Return list object containing design matrix and the matrix structure in blocks
}


# Returns which terms in data matrix contain no interactions
non.interaction.terms <- function(cmodel)
  which(sapply(cmodel$blocks, function(a) length(a$parents))==0)


# Returns the design matrix x for a specified model with active selecting the respective blocks
getX <- function(cmodel, active) {
  col.active <- c(1, do.call(c, lapply(cmodel$blocks[active],"[[","columns")))
  cmodel$X[,col.active,drop=FALSE]
}

#  Functionality to obtain which terms can be considered for removal in the model selection step
candidate.drop.terms <- function(active.terms, X.blocks) {
  select <- sapply(X.blocks[active.terms], function(a) length(intersect(active.terms, a$children))==0)
  if (length(select)==0)
    return(integer(0))    
  active.terms[select]
}

# Functionality to obtain which terms can be considered for addition in the model selection step 
candidate.add.terms <- function(active.terms, X.blocks) {
  stage1 <- setdiff(which(!sapply(X.blocks, is.null)), active.terms) # Obtain the terms which can be added at the model selection step
  if (length(stage1)==0)  # If we have all terms added, return integer(0)
    return(stage1)
  stage1[sapply(X.blocks[stage1], function(a) all(a$parents%in%active.terms))]  # Select which terms to include, making sure all parents are included if an interaction is added
}

# Expand out beta to include terms which are now not active
expand.beta <- function(beta, active.terms, cmodel) {
  ebeta <- list()
  for (j in 1:length(beta)) {
    ebeta[[j]] <- numeric(ncol(cmodel$X))
    col.active <- c(1, do.call(c, lapply(cmodel$blocks[active.terms[[j]]],"[[","columns")))        
    ebeta[[j]][col.active] <- beta[[j]]
  }
  ebeta
}

# Compute the posterior for beta for ONE Y coefficient (Issue: doesn not quite take the cov structure into consideration)
# Output returns the mean and precision, and Cholesky decomposition of precision for generating a sample from the posterior of beta
compute.posterior.beta <- function(Y, X, Xj, beta, j, tau, Lambda.eps) {
  pprec <- 1/tau[j] * diag(ncol(Xj)) + Lambda.eps[j,j] * crossprod(Xj) # 1/tau[j]*I + omega[j,j]*Xj^TXj
  pprec.chol <- chol(pprec)   # Cholesky decomposition of precision
  pcov <- solve(pprec)        # Covariance estimate  
  z <- Lambda.eps[j,j]*Y[,j]  # 1/omega[j,j](Y-reffs)
  for (k in 1:length(beta))
    if (k!=j)   # Computation for the inclusion of the other Y values
      z <- z + Lambda.eps[j,k] * (Y[,k]-X[[k]]%*%beta[[k]]) # z + omega[j,k](prec vals between Yis)*(residuals[k])
  pmean <- as.numeric(backsolve(pprec.chol,forwardsolve(pprec.chol, t(Xj)%*%z, transpose=TRUE, upper.tri=TRUE)))
  list(pmean=pmean, pprec=pprec, pprec.chol=pprec.chol)
}

# Simulate from the posterior of beta using the output from the compute.posterior.beta functionality
simulate.posterior.beta <- function(posterior.result) {
  drop(posterior.result$pmean + backsolve(posterior.result$pprec.chol, rnorm(ncol(posterior.result$pprec.chol))))  # Simulates from MVnormal, but using backsolve and Chol. decompostions
}

# Calculate model evidence
compute.lmodelev <- function(posterior.result,Y, Yhat, Xj, j, tau, Lambda.eps) {
  Yhat.tmp <- Yhat # Set temporary Y.hat value
  Yhat.tmp[,j] <- as.numeric(Xj%*%posterior.result$pmean)
  # Calculate density value by summing over prior on beta and likelihood (removing xbeta - ugamma)
  sum(dnorm(posterior.result$pmean, 0, tau[j], log=TRUE)) + sum(dmvnorm(Y-Yhat.tmp, rep(0, ncol(Y)), solve(Lambda.eps), log=TRUE)) 
}

# Block diagonal matrix function
blockdiag <- function(M) {
  m.out <- list()
  for (i in 1:length(M)) {
    l <- list()
    for (j in 1:length(M))
      if (i==j)
        l[[j]] <- M[[j]]
      else
        l[[j]] <- matrix(0, nrow=nrow(M[[i]]), ncol=ncol(M[[j]]))
      m.out[[i]] <- do.call("cbind", l)
  }
  do.call("rbind", m.out)
}


# Obtain nested X for nested coefficients
getNestedX <- function(vars,active,lookup,cmodel)  {
  reff.lookup <- sapply(lookup,"[",1)
  X.nest <- list()
  cols <- list()
  for(j in 1:length(active))  {
    labels <- NULL
    cols[[j]] <- 1
    if(length(active[[j]])==0) {
      X.nest[[j]] <- as.matrix(cmodel$X[reff.lookup,cols[[j]]])
    }else{
      for(k in 1:length(active[[j]])) {
        labels[k] <- cmodel$blocks[[active[[j]][k]]]$label
      }
      active.terms <- intersect(vars,labels)
      if(length(active.terms>0)){
        for(k in 1:length(active.terms)) {
          for(l in 1:length(cmodel$blocks)) {
            if(active.terms[k]==cmodel$blocks[[l]]$label)
              select.terms <- cmodel$blocks[[l]]$columns
          }
          cols[[j]] <- c(cols[[j]],select.terms)
        }
        cols[[j]] <- cols[[j]][order(cols[[j]])]
        X.nest[[j]] <- as.matrix(cmodel$X[reff.lookup,cols[[j]]])
      }else{
        X.nest[[j]] <- as.matrix(cmodel$X[reff.lookup,cols[[j]]])
      }
    }
  }
  return(list(Xnest=X.nest,columns=cols))
}

# Obtain random effects coefficients 
getReffs <- function(RefLambda,EpsLambda,refcounts,reflookup,workY){
  reffect <- matrix(0,nrow=length(reflookup),ncol=ncol(workY))
  for(j in 1:length(reflookup))  {
    pprec <- RefLambda + EpsLambda * refcounts[j] # Compute pprec for each unique level of random effect
    pprec.chol <- chol(pprec)
    pmean <- as.numeric(backsolve(pprec.chol,forwardsolve(pprec.chol, EpsLambda%*%colSums(workY[reflookup[[j]],,drop=FALSE]), transpose=TRUE, upper.tri=TRUE))) # Compute pmean for current level of random effect
    reffect[j,] <- pmean + drop(backsolve(pprec.chol, rnorm(ncol(pprec.chol)))) # Obtain random effect coefficient
  }
  return(reffect)
}
# Functions to compute the density of a chordal graph (3x3)
# Filename - ChordalGraph.R
# Author - CA


# 2 Graph initialisation
Zero <- matrix(0,nrow=2,ncol=2)
One <- matrix(1,nrow=2,ncol=2) - diag(2)
graphs2 <- list(Zero,One)

# 3 Graph example initialisation
Zero <- matrix(0, nrow=3, ncol=3)
One <- rbind(c(0,1,0),
             c(1,0,0),
             c(0,0,0))
Two <- rbind(c(0,1,1),
             c(1,0,0),
             c(1,0,0))
Three <- matrix(1, nrow=3, ncol=3) - diag(3)
p2 <- c(1,3,2)
p3 <- c(2,3,1)
graphs3 <- list(Zero, One, One[p2,p2], One[p3,p3], Two, Two[p2,p2], Two[p3,p3], Three)

# Functionality to compute the log multivariate gamma function
lmvgamma <- function(a, p) {
  if (p==1)
    return(lgamma(a)) # Returns log of gamma function
  (p-1)/2*log(pi)+lmvgamma(a-0.5,p-1)+lgamma(a)
}

# Density of inverse Wishart distribution, 1st part
IW <- function(b, D) {
  p <- nrow(D)
  result <- log(2) * (b*p)/2 - b/2 *  determinant(D)$modulus +  lmvgamma(b/2, p)
  attributes(result) <- NULL
  result
}

# Full density of inverse wishart 
ldensW <- function(S, b, D) {
  (b-nrow(S)-1)/2*determinant(S)$modulus-sum(diag(S%*%D))/2-IW(b,D)
}

# Functionality to compute the density value for a given graph G, which is 3X3 (to allow chordality results)
IG3 <- function(b, D, G) {
  if (dim(G)[1]!=3)
    stop("Sorry, only 3x3")
  # Independence
  if (all(G==0))
    return(IW(b,D[1,1,drop=FALSE]) + IW(b,D[2,2,drop=FALSE]) + IW(b,D[3,3,drop=FALSE]))
  # Fully connected graph
  if (sum(G)==6)
    return(IW(b,D))
  # One link
  if (sum(G)==2) {
    pair <- min(which(G>0))
    pair <- c((pair-1)%/%3+1,(pair-1)%%3+1)
    return(IW(b,D[pair,pair]) + IW(b,D[-pair,-pair,drop=FALSE]))
  }
  # Remaining case is a chordal graph
  missing.pair <- min(which(G+diag(3)==0))
  missing.pair <- c((missing.pair-1)%/%3+1,(missing.pair-1)%%3+1)
  separator <- setdiff(1:3,missing.pair)
  clique.1 <- c(missing.pair[1], separator)
  clique.2 <- c(missing.pair[2], separator)
  return(IW(b,D[clique.1,clique.1]) + IW(b,D[clique.2,clique.2]) - IW(b,D[separator, separator, drop=FALSE]))
}

# Compute the log evidence for graph selection, combining the log likelihood of data and the G-wishart density
graphlgev3 <- function(Y, G, b, D) {
  S <- crossprod(Y)
  Omega <- rGWish_sampler(1, b+nrow(Y), solve(D+S), G, 100)[,,1]
  result <- sum(dmvnorm(Y, rep(0,ncol(Y)), solve(Omega), log=TRUE)) + (b-2)/2 * determinant(Omega)$modulus - 1/2 * sum(diag((D)%*%Omega)) - IG3(b, D, G) - (b+nrow(Y)-2)/2 * determinant(Omega)$modulus + 1/2 * sum(diag((D+S)%*%Omega)) + IG3(b+nrow(Y), D+S, G)
  attributes(result) <- NULL
  result
}

# 2 Graph case
IG2 <- function(b, D, G) {
  if(dim(G)[1]!=2)
    stop("Sorry! Only 2x2")
  
  # Independent graph
  if(all(G==0))
    return(IW(b,D[1,1,drop=FALSE]) + IW(b,D[2,2,drop=FALSE]))
  
  # Fully connected graph
  if(sum(G)==2)
    return(IW(b,D))
    
}

# Compute the log evidence for graph selection, combining the log likelihood of data and the G-wishart density
graphlgev2 <- function(Y, G, b, D) {
  S <- crossprod(Y)
  Omega <- rGWish_sampler(1, b+nrow(Y), solve(D+S), G, 100)[,,1]
  result <- sum(dmvnorm(Y, rep(0,ncol(Y)), solve(Omega), log=TRUE)) + (b-2)/2 * determinant(Omega)$modulus - 1/2 * sum(diag((D)%*%Omega)) - IG2(b, D, G) - (b+nrow(Y)-2)/2 * determinant(Omega)$modulus + 1/2 * sum(diag((D+S)%*%Omega)) + IG2(b+nrow(Y), D+S, G)
  attributes(result) <- NULL
  result
}
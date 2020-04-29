# Multivariate parameter expansion example
# Filename: ParameterExMultivar.R
# Author: CA
# Date: 28/01/18

# Libraries
library(MASS)
library(Matrix)
library(MCMCpack)
library(mvtnorm)

# Block diagonal function
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

# Generate simulated data with 1 constant fixed effect and a nested fixed effect with 1 random effect
true.Sigma <- cbind(c(10,0,0),c(0,10,0),c(0,0,10)) # Set covariance structure (Assuming 3 responses here)

n <- 100 # # of observations


X <- matrix(rep(1,n))
X.list <- list(X,X,X)
X.tilde <- blockdiag(X.list)

tau.beta.prior <- 0.0001
true.beta <- matrix(rnorm(3,sd=sqrt(1/tau.beta.prior)))

Y <- as.matrix(c(X.tilde%*%true.beta)+ matrix(rnorm(n*3),ncol=3)%*%chol(true.Sigma))


# Generate random effects
true.Sigma.gamma <- cbind(c(0.1,0,0),c(0,0.1,0),c(0,0,0.1))
true.gamma <- mvrnorm(0.9*n,mu=c(0,0,0),true.Sigma.gamma)

U1 <- diag(1,0.8*n)
U2 <- matrix(0,nrow=0.2*n,ncol=0.1*n)
for(i in 1:ncol(U2))  {
  U2[((2*i)-1):(2*i),i] <- 1
}

Ulist <- list(U1,U2)
U <- blockdiag(Ulist)
U.list <- list(U,U,U)
U.tilde <- blockdiag(U.list)



Y <- matrix(c(Y) + U.tilde%*%c(true.gamma),ncol=3)

##              ##
# Hybrid Sampler #
##              ##

# Initialise values
samples <- 1000
beta <- rep(0,nrow(true.beta)*ncol(Y))
gamma <- matrix(rep(0,nrow(true.gamma)*ncol(Y)),ncol=ncol(Y))
Omega.gamma <- diag(ncol(Y))
Omega.error <- diag(ncol(Y))

prior.prec.beta <- diag(tau.beta.prior,3)
d0.gamma <- 3
d0.error <- 3
V0.gamma <- diag(1,ncol(Y))
V0.error <- diag(1,ncol(Y))   


alpha <- matrix(rep(0,nrow(gamma)*ncol(Y)),ncol=ncol(Y))
X.alpha <- matrix(rep(1,nrow(gamma)),ncol=1)
Xalphalist <- list(X.alpha,X.alpha,X.alpha)
Xalpha <- blockdiag(Xalphalist)

beta.res <- matrix(0,nrow=samples,ncol=length(true.beta))
gamma.res <- matrix(0,nrow=samples,ncol=length(true.gamma))
Omega.gamma.res <- vector("list",samples)
Omega.error.res <- vector("list",samples)
alpha.res <- matrix(0,nrow=samples,ncol=length(true.gamma))
acceptance.rate <- rep(0,samples)

df.gamma <- nrow(gamma) + d0.gamma
df.error <- nrow(Y) + d0.error

phi.alpha <-10
phi.beta <- 10

# Sampler
for(i in 1:samples) {
  
  # Sample beta update
  sigma.beta <- solve(t(X.tilde)%*%kronecker(Omega.error,Diagonal(nrow(X)))%*%X.tilde + prior.prec.beta)
  mean.beta <- sigma.beta%*%(t(X.tilde)%*%kronecker(Omega.error,Diagonal(nrow(X)))%*%(as.vector(Y)-matrix(U.tilde%*%c(gamma))))
  beta <- mvrnorm(1,mean.beta,sigma.beta)
  beta <- matrix(beta,ncol=3)
  
  
  # Sample gamma update
  sigma.gamma <- solve(t(U.tilde)%*%kronecker(Omega.error,Diagonal(nrow(X)))%*%U.tilde + kronecker(Omega.gamma,Diagonal(nrow(gamma))))
  mean.gamma <- sigma.gamma%*%t(U.tilde)%*%kronecker(Omega.error,Diagonal(nrow(X)))%*%(as.vector(Y) - matrix(X.tilde%*%c(beta)))
  gamma <- mvrnorm(1,mean.gamma,sigma.gamma)
  gamma <- matrix(gamma,ncol=3)
  
  # Form alpha update
  #beta.nest <- beta
  #alpha <- matrix(Xalpha%*%c(beta.nest) + c(gamma),ncol=3)
  
  # Sample beta update conditional on alpha
  #sigma.beta.alpha <- solve(t(Xalpha)%*%(Xalpha%*%(diag(tau.beta.prior,3))%*%t(Xalpha)+kronecker(Omega.gamma,Diagonal(nrow(gamma))))%*%Xalpha+ diag(tau.beta.prior,3))
  #mean.beta.alpha <- sigma.beta.alpha%*%t(Xalpha)%*%(Xalpha%*%(diag(tau.beta.prior,3))%*%t(Xalpha)+kronecker(Omega.gamma,Diagonal(nrow(gamma))))%*%(c(alpha))
  #beta.nest <- mvrnorm(1,mean.beta.alpha,sigma.beta.alpha)
  #beta.nest <- matrix(beta.nest,ncol=3)
  
  # Update beta with beta nest coefficients
  #beta <- beta.nest
  
  
  # Obtain gamma coefficents
  #gamma <- c(alpha) - Xalpha%*%c(beta.nest)
  #gamma <- matrix(gamma,ncol=3)
  
  # Sample gamma precision update
  Scale.gamma <- solve(crossprod(gamma) + (V0.gamma))
  Omega.gamma <- rwish(df.gamma,(Scale.gamma))
  
  # Proposal step for new gamma updates via parameter expansion

  # Draw phi from the inverse gamma distribution
  phi <- rgamma(1,phi.alpha,phi.beta)
  D.phi <- diag(phi,ncol(Y))
  
  # Form effect.star and precision.star
  effect.star <- gamma*phi
  #prec.star <- t(D.phi)%*%Omega.gamma%*%D.phi
  prec.star <- (phi^2)*Omega.gamma
  
  Ywork <- matrix(c(Y) - X.tilde%*%c(beta) - U.tilde%*%c(gamma),ncol=3)
  
  # Compute acceptance rate
  original.rate <- dgamma(phi,phi.alpha,phi.beta,log=T)+sum(dmvnorm(Ywork,rep(0,ncol(Y)),solve(Omega.error)),log=TRUE)+
    sum(dmvnorm(gamma,rep(0,ncol(Y)),solve(Omega.gamma)),log=TRUE)#+
    #log(dwish(Omega.gamma,d0.gamma,V0.gamma))
  
  Ystar <- NULL
  Ystar <- c(Ywork) + U.tilde%*%c(gamma)
  Ystar <- matrix(Ystar - U.tilde%*%c(effect.star),ncol=3)
  
  star.rate <- dgamma(1/phi,phi.alpha,phi.beta,log=T)+sum(dmvnorm(Ystar,rep(0,ncol(Y)),solve(Omega.error)),log=TRUE)+ 
    sum(dmvnorm(effect.star,rep(0,ncol(Y)),solve(prec.star)),log=TRUE)#+
    #log(dwish(prec.star,d0.gamma,V0.gamma))
  
  factor <- ((ncol(Y)*(ncol(Y)+1))/2)*log(phi)
  
  acc.rate <- star.rate - original.rate + factor
  
  #acc.rate <- -Inf
  
  if(log(runif(1)) < acc.rate) {
    gamma <- effect.star
    Omega.gamma <- prec.star
    acceptance.rate[i] <- 1
  }else{
    acceptance.rate[i] <- 0
  }
  
  # Sample error precision update
  Yerror <- matrix(c(Y) - matrix(X%*%beta) - matrix(U.tilde%*%c(gamma)), ncol=3) 
  Scale.error <- solve(crossprod(Yerror) + (V0.error))
  Omega.error <- rWishart(1,df.error,(Scale.error))[,,1]
  
  # Store values
  beta.res[i,] <- as.vector(beta)
  gamma.res[i,] <- as.vector(gamma)
  alpha.res[i,] <- as.vector(alpha)
  Omega.gamma.res[[i]] <- Omega.gamma
  Omega.error.res[[i]] <- Omega.error
  
}

precs2 <- vector(length=3,mode='list')
for(i in 1:3) {
  for(j in 1:samples) {
    precs2[[i]][j] <- (Omega.gamma.res[[j]][i,i])
  }
}


par(mfrow=c(3,1))
for(i in 1:3)
  plot(precs2[[i]],type='l',xlab="",ylab=expression(sigma[3]^2),xaxt='n')


reffs2<- vector(length=89, mode='list')
for(i in 1:89)  {
  for(j in 1:samples)  {
    reffs2[[i]][j] <- gamma.res[j,i]
  }
}

par(mfrow=c(3,1))
for(i in 1:3)
  plot(reffs2[[3]],type='l',xlab="",ylab=expression(gamma[3]),xaxt='n')

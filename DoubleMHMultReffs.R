# Graphical model selection by double M-H algorithm for multiple precision estimates
# Filename - DoubleMHMultReffs.R
# Author - CA

# Source files
sourceCpp("GWishGen.cpp")

# Function for computing J(b,D,a)
J.b <- function(b,D,A)  {
  J.b1l <- 0.5*(log(2) + log(pi) - log(D[2,2])) + ((b-1)/2)*log(A[1,1])
  I.bDel <- -(b/2)*log(2) -log(D[2,2]) - lgamma((b/2))
  I.bDel <- -(b+1)*log(D[2,2])+lgamma(b+1)  
  J.b2l <- -0.5*(D[1,1] - ((1/D[2,2])*(D[1,2]^2)))*A[1,1]
  J.bl <- J.b1l + I.bDel + J.b2l
  return(J.bl)
}

# GM selection code
double.MH <- function(dataset,reffs,iterations=10000, D.eps=diag(ncol(dataset)), b.eps = ncol(dataset), D.reffs=rep(list(diag(ncol(dataset))),length(reffs)),b.reffs=rep(ncol(dataset),length(reffs)),G = NULL) {
  
  # Storage items
  p.reffs <- rep(0,length(reffs))
  n.reffs <- rep(0,length(reffs))
  b.star.reffs <- rep(0,length(reffs))
  D.star.reffs <- vector(length=length(reffs),mode='list')
  Omega.reffs.vals <- list()
  S.reffs <- vector(length=length(reffs),mode='list')
  Omega.reffs <- vector(length=length(reffs),mode='list')
  Omega.prime.reffs <- vector(length=length(reffs),mode='list')
  A.prime.reffs <- vector(length=length(reffs),mode='list')
  A.reffs <- vector(length=length(reffs),mode='list')
  I.bdl.reffs <- rep(0,length(reffs))
  f.Gprimel.reffs <- rep(0,length(reffs))
  J.b1.reffs <- rep(0,length(reffs))
  f.Gl.reffs <- rep(0,length(reffs))
  I.bnl.reffs <- rep(0,length(reffs))
  J.Dn.reffs <- rep(0,length(reffs))
  H.el.reffs <- rep(0,length(reffs))
  H.eOl.reffs <- rep(0,length(reffs))
  
  # Initialise values
  p.eps <- ncol(dataset)
  n.eps <- nrow(dataset)
  S.eps <- crossprod(dataset)  
  for(j in 1:length(reffs)) {
    p.reffs[j] <- ncol(reffs[[j]])
    n.reffs[j] <- nrow(reffs[[j]])
    S.reffs[[j]] <- crossprod(reffs[[j]])
  }
  
  # G-Wishart sampling parameters
  b.star.eps <- b.eps + n.eps
  D.star.eps <- (D.eps + S.eps)
  for(j in 1:length(reffs)) {
    b.star.reffs[j] <- n.reffs[j] + b.reffs[j]
    D.star.reffs[[j]] <- S.reffs[[j]] + D.reffs[[j]]
  }
  
  # Initialise parameters for storing results
  Omega.eps.vals <- vector(length=iterations,mode='list')
  for(j in 1:length(reffs)) {
    Omega.reffs.vals[[j]] <- vector(length=iterations,mode='list')
  }
  G.vals     <- vector(length=iterations,mode='list')
  acc.rate   <- NULL
  
  
  
  # Initial precision estimate
  Omega.eps <- rGWish_sampler(1,b.star.eps,solve(D.star.eps),G,100)[,,1]
  for(j in 1:length(reffs)) {
    Omega.reffs[[j]] <- rGWish_sampler(1,b.star.reffs[j],solve(D.star.reffs[[j]]),G,100)[,,1]
  }
  
  # Begin MCMC
  for(i in 1:iterations)  {
    
    # 1(a). Propose new graph G.prime
    sel <- sample(p.eps,2)
    J <- min(sel)
    K <- max(sel)
    G.prime <- G
    if(G[J,K] == 1) {
      G.prime[J,K]  <- G.prime[K,J] <- 0
      
    # 1(b). Generate Omega.prime from new graph G.prime
    Omega.prime.eps <- rGWish_sampler(1,b.eps,solve(D.eps),G.prime,100)[,,1]
    for(j in 1:length(reffs)) {
       Omega.prime.reffs[[j]] <- rGWish_sampler(1, b.reffs[j],solve(D.reffs[[j]]),G.prime,100)[,,1]
     }
      
    # 2. Updating 
    # Create Omega.prime.0
    Omega.prime.0.eps <- Omega.prime.eps
    Omega.prime.0.reffs <- Omega.prime.reffs
    Omega.prime.0.eps[J,K] <-Omega.prime.0.eps[K,J]<- 0
    for(j in 1:length(reffs)) {
      Omega.prime.0.reffs[[j]][J,K] <- Omega.prime.0.reffs[[j]][K,J] <- 0
    }
    Omega.prime.0.eps[K,K] <- Omega.prime.eps[K,-K]%*%solve(Omega.prime.eps[-K,-K],Omega.prime.eps[-K,K])
    for(j in 1:length(reffs)){
      Omega.prime.0.reffs[[j]][K,K] <- Omega.prime.reffs[[j]][K,-K]%*%solve(Omega.prime.reffs[[j]][-K,-K],Omega.prime.reffs[[j]][-K,K])
    }
      
    # Create Omega.prime.1
    Omega.prime.1.eps <- Omega.prime.eps
    Omega.prime.1.reffs <- Omega.prime.reffs
    Omega.prime.1.eps[c(J,K),c(J,K)] <- Omega.prime.eps[c(J,K),-c(J,K)]%*%solve(Omega.prime.eps[-c(J,K),-c(J,K)])%*%Omega.prime.eps[-c(J,K),c(J,K)]
    for(j in 1:length(reffs)) {
      Omega.prime.1.reffs[[j]][c(J,K),c(J,K)] <- Omega.prime.reffs[[j]][c(J,K),-c(J,K)]%*%solve(Omega.prime.reffs[[j]][-c(J,K),-c(J,K)])%*%Omega.prime.reffs[[j]][-c(J,K),c(J,K)]
    }
    
    # Calculate f(Omega.prime|G.prime)
    I.bdl.eps <- -(b.eps+1)*log(D.eps[K,K]) + lgamma(b.eps+1)
    for(j in 1:length(reffs)) {
      I.bdl.reffs[j] <- -(b.reffs[j]+1)*log(D.reffs[[j]][K,K]) + lgamma(b.reffs[j]+1)
    }
      
    f.Gprimel.eps <- I.bdl.eps + (((b.eps-2)/2)*determinant(Omega.prime.0.eps[-K,-K],drop=FALSE)$modulus) - (0.5*sum(diag(D.eps%*%Omega.prime.0.eps)))
    for(j in 1:length(reffs)) {
      f.Gprimel.reffs[j] <- I.bdl.reffs[j] + (((b.reffs[j]-2)/2)*determinant(Omega.prime.0.reffs[[j]][-K,-K],drop=F)$modulus) - (0.5*sum(diag(D.reffs[[j]]%*%Omega.prime.0.reffs[[j]])))
    }
    
    f.Gprimel <- f.Gprimel.eps + sum(f.Gprimel.reffs)
      
    # Calculate f(Omega.prime|G)
    A.prime.eps <- Omega.prime.eps[c(J,K),c(J,K)] - (Omega.prime.eps[c(J,K),-c(J,K),drop=F]%*%solve(Omega.prime.eps[-c(J,K),-c(J,K),drop=F])%*%Omega.prime.eps[-c(J,K),c(J,K),drop=F])
    for(j in 1:length(reffs)) {
      A.prime.reffs[[j]] <- Omega.prime.reffs[[j]][c(J,K),c(J,K)] - (Omega.prime.reffs[[j]][c(J,K),-c(J,K),drop=F]%*%solve(Omega.prime.reffs[[j]][-c(J,K),-c(J,K),drop=F])%*%Omega.prime.reffs[[j]][-c(J,K),c(J,K),drop=F])
    }
      
    J.b1.eps <- J.b(b.eps,D.eps[c(J,K),c(J,K)],A.prime.eps)
    for(j in 1:length(reffs)) {
      J.b1.reffs[j] <- J.b(b.reffs[j],D.reffs[[j]][c(J,K),c(J,K)],A.prime.reffs[[j]])
    }

    f.Gl.eps <- J.b1.eps +(((b.eps-2)/2)*determinant(Omega.prime.1.eps[-c(J,K),-c(J,K),drop=FALSE])$modulus) + (-0.5*sum(diag(D.eps%*%Omega.prime.1.eps)))
    for(j in 1:length(reffs)) {
      f.Gl.reffs[j] <- J.b1.reffs[j] +(((b.reffs[j]-2)/2)*determinant(Omega.prime.1.reffs[[j]][-c(J,K),-c(J,K),drop=F])$modulus) + (-0.5*sum(diag(D.reffs[[j]]%*%Omega.prime.1.reffs[[j]])))
    }
      
    f.Gl <- f.Gl.eps + sum(f.Gl.reffs)
      
    # Create Omega.0
    Omega.0.eps <- Omega.eps
    Omega.0.reffs <- Omega.reffs
      
    Omega.0.eps[J,K] <- Omega.0.eps[K,J] <- 0
    for(j in 1:length(reffs)) {
      Omega.0.reffs[[j]][J,K] <- Omega.0.reffs[[j]][K,J] <- 0
    }
    
    Omega.0.eps[K,K] <- Omega.eps[K,-K]%*%solve(Omega.eps[-K,-K],Omega.eps[-K,K])
    for(j in 1:length(reffs)) {
      Omega.0.reffs[[j]][K,K] <- Omega.reffs[[j]][K,-K]%*%solve(Omega.reffs[[j]][-K,-K],Omega.reffs[[j]][-K,K]) 
    }
      
    # Create Omega.1
    Omega.1.eps <- Omega.eps
    Omega.1.reffs <- Omega.reffs
    
    Omega.1.eps[c(J,K),c(J,K)] <- Omega.eps[c(J,K),-c(J,K)]%*%solve(Omega.eps[-c(J,K),-c(J,K)])%*%Omega.eps[-c(J,K),c(J,K)]
    for(j in 1:length(reffs)) {
      Omega.1.reffs[[j]][c(J,K),c(J,K)] <- Omega.reffs[[j]][c(J,K),-c(J,K)]%*%solve(Omega.reffs[[j]][-c(J,K),-c(J,K)])%*%Omega.reffs[[j]][-c(J,K),c(J,K)]
    }

    # Calculate H(e,Omega)
    A.eps    <- Omega.eps[c(J,K),c(J,K)] - (Omega.eps[c(J,K),-c(J,K),drop=F]%*%solve(Omega.eps[c(-J,-K),c(-J,-K),drop=F])%*%Omega.eps[-c(J,K),c(J,K),drop=F])
    for(j in 1:length(reffs)) {
      A.reffs[[j]] <- Omega.reffs[[j]][c(J,K),c(J,K)] - (Omega.reffs[[j]][c(J,K),-c(J,K),drop=F]%*%solve(Omega.reffs[[j]][c(-J,-K),c(-J,-K),drop=F])%*%Omega.reffs[[j]][-c(J,K),c(J,K),drop=F])
    }
      
    I.bnl.eps <- -(b.star.eps+1)*log(D.star.eps[K,K]) + lgamma(b.star.eps+1)
    for(j in 1:length(reffs)) {
      I.bnl.reffs[j] <- (b.star.reffs[j]+1)*log(D.star.reffs[[j]][K,K]) + lgamma(b.star.reffs[j]+1)
    }
      
    J.Dn.eps <- J.b(b.star.eps,D.star.eps[c(J,K),c(J,K)],A.eps)
    for(j in 1:length(reffs)) {
      J.Dn.reffs[j] <- J.b(b.star.reffs[j],D.star.reffs[[j]][c(J,K),c(J,K)],A.reffs[[j]])
    }

    H.el.eps <- (I.bnl.eps - J.Dn.eps) + ((b.star.eps-2)/2)*(determinant(Omega.0.eps[-K,-K,drop=FALSE])$modulus - determinant(Omega.1.eps[-c(J,K),-c(J,K),drop=FALSE])$modulus)
    for(j in 1:length(reffs)) {
      H.el.reffs[j] <- (I.bnl.reffs[j] - J.Dn.reffs[j]) + ((b.star.reffs[j]-2)/2)*(determinant(Omega.0.reffs[[j]][-K,-K,drop=F])$modulus - determinant(Omega.1.reffs[[j]][-c(J,K),-c(J,K),drop=F])$modulus)
    }
      
    H.eOl.eps <- H.el.eps - (0.5*sum(diag(D.star.eps%*%(Omega.0.eps-Omega.1.eps))))
    for(j in 1:length(reffs)) {
      H.eOl.reffs[j] <- H.el.reffs[j] - (0.5*sum(diag(D.star.reffs[[j]]%*%(Omega.0.reffs[[j]]-Omega.1.reffs[[j]]))))
    }

      
    H.eOl <- H.eOl.eps + sum(H.eOl.reffs)
      
      
    alpha.G <- (f.Gl - f.Gprimel) + H.eOl
      
    if(alpha.G > log(runif(1))) {
      
      # Update G
      G <- G.prime
      acc.rate[i] <- 1
      }else{
        acc.rate[i] <- 0
      }
      
    }else {
      
      # Else we add an edge 1(a).
      G.prime[J,K] <- G.prime[K,J] <- 1
      
      # 1(b). Generate Omega.prime from new graph G.prime
      Omega.prime.eps <- rGWish_sampler(1,b.eps,solve(D.eps),G.prime,100)[,,1]
      for(j in 1:length(reffs)) {
        Omega.prime.reffs[[j]] <- rGWish_sampler(1, b.reffs[j],solve(D.reffs[[j]]),G.prime,100)[,,1]
      }
      
      # 2. Updating 
      # Create Omega.prime.0
      Omega.prime.0.eps <- Omega.prime.eps
      Omega.prime.0.reffs <- Omega.prime.reffs
      Omega.prime.0.eps[J,K] <-Omega.prime.0.eps[K,J]<- 0
      for(j in 1:length(reffs)) {
        Omega.prime.0.reffs[[j]][J,K] <- Omega.prime.0.reffs[[j]][K,J] <- 0
      }
      Omega.prime.0.eps[K,K] <- Omega.prime.eps[K,-K]%*%solve(Omega.prime.eps[-K,-K],Omega.prime.eps[-K,K])
      for(j in 1:length(reffs)){
        Omega.prime.0.reffs[[j]][K,K] <- Omega.prime.reffs[[j]][K,-K]%*%solve(Omega.prime.reffs[[j]][-K,-K],Omega.prime.reffs[[j]][-K,K])
      }
      
      # Create Omega.prime.1
      Omega.prime.1.eps <- Omega.prime.eps
      Omega.prime.1.reffs <- Omega.prime.reffs
      Omega.prime.1.eps[c(J,K),c(J,K)] <- Omega.prime.eps[c(J,K),-c(J,K)]%*%solve(Omega.prime.eps[-c(J,K),-c(J,K)])%*%Omega.prime.eps[-c(J,K),c(J,K)]
      for(j in 1:length(reffs)) {
        Omega.prime.1.reffs[[j]][c(J,K),c(J,K)] <- Omega.prime.reffs[[j]][c(J,K),-c(J,K)]%*%solve(Omega.prime.reffs[[j]][-c(J,K),-c(J,K)])%*%Omega.prime.reffs[[j]][-c(J,K),c(J,K)]
      }
      
      # Calculate f(Omega.prime|G.prime)
      I.bdl.eps <- -(b.eps+1)*log(D.eps[K,K]) + lgamma(b.eps+1)
      for(j in 1:length(reffs)) {
        I.bdl.reffs[j] <- -(b.reffs[j]+1)*log(D.reffs[[j]][K,K]) + lgamma(b.reffs[j]+1)
      }
      
      f.Gprimel.eps <- I.bdl.eps + (((b.eps-2)/2)*determinant(Omega.prime.0.eps[-K,-K],drop=FALSE)$modulus) - (0.5*sum(diag(D.eps%*%Omega.prime.0.eps)))
      for(j in 1:length(reffs)) {
        f.Gprimel.reffs[j] <- I.bdl.reffs[j] + (((b.reffs[j]-2)/2)*determinant(Omega.prime.0.reffs[[j]][-K,-K],drop=F)$modulus) - (0.5*sum(diag(D.reffs[[j]]%*%Omega.prime.0.reffs[[j]])))
      }
      
      f.Gprimel <- f.Gprimel.eps + sum(f.Gprimel.reffs)
      
      # Calculate f(Omega.prime|G)
      A.prime.eps <- Omega.prime.eps[c(J,K),c(J,K)] - (Omega.prime.eps[c(J,K),-c(J,K),drop=F]%*%solve(Omega.prime.eps[-c(J,K),-c(J,K),drop=F])%*%Omega.prime.eps[-c(J,K),c(J,K),drop=F])
      for(j in 1:length(reffs)) {
        A.prime.reffs[[j]] <- Omega.prime.reffs[[j]][c(J,K),c(J,K)] - (Omega.prime.reffs[[j]][c(J,K),-c(J,K),drop=F]%*%solve(Omega.prime.reffs[[j]][-c(J,K),-c(J,K),drop=F])%*%Omega.prime.reffs[[j]][-c(J,K),c(J,K),drop=F])
      }
      
      J.b1.eps <- J.b(b.eps,D.eps[c(J,K),c(J,K)],A.prime.eps)
      for(j in 1:length(reffs)) {
        J.b1.reffs[j] <- J.b(b.reffs[j],D.reffs[[j]][c(J,K),c(J,K)],A.prime.reffs[[j]])
      }
      
      f.Gl.eps <- J.b1.eps +(((b.eps-2)/2)*determinant(Omega.prime.1.eps[-c(J,K),-c(J,K),drop=FALSE])$modulus) + (-0.5*sum(diag(D.eps%*%Omega.prime.1.eps)))
      for(j in 1:length(reffs)) {
        f.Gl.reffs[j] <- J.b1.reffs[j] +(((b.reffs[j]-2)/2)*determinant(Omega.prime.1.reffs[[j]][-c(J,K),-c(J,K),drop=F])$modulus + (-0.5*sum(diag(D.reffs[[j]]%*%Omega.prime.1.reffs[[j]]))))
      }
      
      f.Gl <- f.Gl.eps + sum(f.Gl.reffs)
      
      # Create Omega.0
      Omega.0.eps <- Omega.eps
      Omega.0.reffs <- Omega.reffs
      
      Omega.0.eps[J,K] <- Omega.0.eps[K,J] <- 0
      for(j in 1:length(reffs)) {
        Omega.0.reffs[[j]][J,K] <- Omega.0.reffs[[j]][K,J] <- 0
      }
      
      Omega.0.eps[K,K] <- Omega.eps[K,-K]%*%solve(Omega.eps[-K,-K],Omega.eps[-K,K])
      for(j in 1:length(reffs)) {
        Omega.0.reffs[[j]][K,K] <- Omega.reffs[[j]][K,-K]%*%solve(Omega.reffs[[j]][-K,-K],Omega.reffs[[j]][-K,K]) 
      }
      
      # Create Omega.1
      Omega.1.eps <- Omega.eps
      Omega.1.reffs <- Omega.reffs
      
      Omega.1.eps[c(J,K),c(J,K)] <- Omega.eps[c(J,K),-c(J,K)]%*%solve(Omega.eps[-c(J,K),-c(J,K)])%*%Omega.eps[-c(J,K),c(J,K)]
      for(j in 1:length(reffs)) {
        Omega.1.reffs[[j]][c(J,K),c(J,K)] <- Omega.reffs[[j]][c(J,K),-c(J,K)]%*%solve(Omega.reffs[[j]][-c(J,K),-c(J,K)])%*%Omega.reffs[[j]][-c(J,K),c(J,K)]
      }
      
      # Calculate H(e,Omega)
      A.eps    <- Omega.eps[c(J,K),c(J,K)] - (Omega.eps[c(J,K),-c(J,K),drop=F]%*%solve(Omega.eps[c(-J,-K),c(-J,-K),drop=F])%*%Omega.eps[-c(J,K),c(J,K),drop=F])
      for(j in 1:length(reffs)) {
        A.reffs[[j]] <- Omega.reffs[[j]][c(J,K),c(J,K)] - (Omega.reffs[[j]][c(J,K),-c(J,K),drop=F]%*%solve(Omega.reffs[[j]][c(-J,-K),c(-J,-K),drop=F])%*%Omega.reffs[[j]][-c(J,K),c(J,K),drop=F])
      }
      
      I.bnl.eps <- -(b.star.eps+1)*log(D.star.eps[K,K]) + lgamma(b.star.eps+1)
      for(j in 1:length(reffs)) {
        I.bnl.reffs[j] <- (b.star.reffs[j]+1)*log(D.star.reffs[[j]][K,K]) + lgamma(b.star.reffs[j]+1)
      }
      
      J.Dn.eps <- J.b(b.star.eps,D.star.eps[c(J,K),c(J,K)],A.eps)
      for(j in 1:length(reffs)) {
        J.Dn.reffs[j] <- J.b(b.star.reffs[j],D.star.reffs[[j]][c(J,K),c(J,K)],A.reffs[[j]])
      }
      
      H.el.eps <- (I.bnl.eps - J.Dn.eps) + ((b.star.eps-2)/2)*(determinant(Omega.0.eps[-K,-K,drop=FALSE])$modulus - determinant(Omega.1.eps[-c(J,K),-c(J,K),drop=FALSE])$modulus)
      for(j in 1:length(reffs)) {
        H.el.reffs[j] <- (I.bnl.reffs[j] - J.Dn.reffs[j]) + ((b.star.reffs[j]-2)/2)*(determinant(Omega.0.reffs[[j]][-K,-K,drop=F])$modulus - determinant(Omega.1.reffs[[j]][-c(J,K),-c(J,K),drop=F])$modulus)
      }
      
      H.eOl.eps <- H.el.eps - (0.5*sum(diag(D.star.eps%*%(Omega.0.eps-Omega.1.eps))))
      for(j in 1:length(reffs)) {
        H.eOl.reffs[j] <- H.el.reffs[j] - (0.5*sum(diag(D.star.reffs[[j]]%*%(Omega.0.reffs[[j]]-Omega.1.reffs[[j]]))))
      }
      
      
      H.eOl <- H.eOl.eps + sum(H.eOl.reffs)
      
      alpha.G <- (f.Gl - f.Gprimel) + H.eOl
      
      if(alpha.G > log(runif(1))) { 
        # Update G
        G <- G.prime
        acc.rate[i] <- 1
      }else{
        acc.rate[i] <- 0
      }
    }
    
    # Sample Omega from G-Wshart sampler with new G
    Omega.eps <- rGWish_sampler(1,b.star.eps,solve(D.star.eps),G,100)[,,1]
    for(j in 1:length(reffs)) {
      Omega.reffs[[j]] <- rGWish_sampler(1,b.star.reffs[j],solve(D.star.reffs[[j]]),G,100)[,,1]
    }
    
    # Store results
    Omega.eps.vals[[i]] <- Omega.eps
    for(j in 1:length(reffs)) {
      Omega.reffs.vals[[j]][[i]] <- Omega.reffs[[j]]
    }
    
    G.vals[[i]]    <- G
  }
  
  # Obtain acceptance rate
  acceptance.rate <- sum(acc.rate)/length(acc.rate)
  
  # Return output
  return(list(Precision.eps= Omega.eps.vals,Precision.reffs= Omega.reffs.vals, Graph = G.vals,acceptancerate = acceptance.rate))
}
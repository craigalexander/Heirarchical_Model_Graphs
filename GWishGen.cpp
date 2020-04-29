// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

cube rWish_sampler(int n, int delta, mat D) {
  //GetRNGstate();
  // Initialisation of parameters
  int p = D.n_rows;
  mat Tm = chol(D,"lower");
  cube res(p,p,n,fill::zeros);

//   Rcout << Tm;

  // Start simulation loop
  for(int i = 0; i < n; i++)  {

    mat K(p,p, fill::zeros);

    // Using the Bartlett decomposition, can sample directly a Wishart distribution
    // Sample the diagonal element gamma_i by taking the square root of a chi-squared distributed variate with parameter (delta - i +1)
    for(int j=0; j < p; j++)  {
      K(j,j) = sqrt(R::rchisq(delta - (j+1) +1));  //REMEMBER 0 INDEXING!!
    }

      // Sample the upper off-diagonal elements from a standard normal distribution
      for(int k=1; k < p; k++)  {
        for(int l=0; l < k; l++)  {
          K(k,l) = R::rnorm(0,1);
        }
      }
      mat C = Tm*K;
      res.slice(i) = C*C.t();//C.t() * C;

  }
  return(res);
//PutRNGstate();
}

//[[Rcpp::export]]

cube rGWish_sampler(int n, int delta, mat D, mat G, int conv) {
 //GetRNGstate();
  // Initialisation of parameters
  int p = G.n_rows;
  mat Tm = chol(D, "upper");
  cube res(p, p, n, fill::zeros);

  // Begin simulation loop
  for(int i = 0; i < n; i++)  {
    mat Sigma = inv(rWish_sampler(1,delta,D).slice(0));

    // Step 1: Set W = Sigma
    mat W = Sigma;

    // Step 2: For k=1,...,K (K = the total number of nodes)
    // j is used  here for convergence
    for(int j=0; j < conv; j++) {
      mat W_past = W;     // Set W.past to be W for convergence check
      for(int k=0; k < p; k++)  {
        // a. Form W_Nk and Sigma_Nk,k to obtain beta_k_star = W_Nk^-1 %*% Sigma_Nk,k
        // N_k denotes the set of neighbours of node k as defined in G
         uvec Nk = find(G.col(k) > 0);
         uvec veck(1);
         veck(0)=k;
         vec beta_k_star = inv(W.submat(Nk,Nk)) * Sigma.submat(Nk,veck);

         // b. Form beta_k by copying the elemnts of beta_k_star to the appropriate locations and
         // placing zeroes in those not connected to k in G
         vec beta_k(p,fill::zeros);
         beta_k(Nk) = beta_k_star;

         // c. Replace W[k,-k] and W[-k,k] with W[-j,-j]beta_k
	 vec replacement = W*beta_k;
	 replacement(k)=W(k,k);
	 W.row(k)=replacement.t();
	 W.col(k)=replacement;
      }

      if (max(max(abs(W-W_past)))<sqrt(DBL_EPSILON))
          break;
    }
    res.slice(i) = inv(W);
  }
  return res;
  //PutRNGstate();
}

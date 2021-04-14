// #define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


/*
  Pr(not detected) computed using Monte Carlo integration, by 
  averaged over possible paths. 
*/


// [[Rcpp::export]]
double get_qstar_cam_nomove(int nsim, arma::mat x, int K,
			    double lam0, double kappa,
			    double buffer, double maxDistSq,
			    arma::ivec trapOper,
			    int nthreads) {

  double xmin=x.col(0).min()-buffer;
  double xmax=x.col(0).max()+buffer;
  double ymin=x.col(1).min()-buffer;
  double ymax=x.col(1).max()+buffer;
  int J=x.n_rows;
  double qstar=0.0;
  omp_set_num_threads(nthreads);
  #pragma omp parallel for reduction(+:qstar)
  for(int i=0; i<nsim; i++) {
    double s1 = R::runif(xmin, xmax);
    double s2 = R::runif(ymin, ymax);
    double distSq = 0.0;
    double lam = 0.0;
    for(int j=0; j<J; j++) {
      distSq = pow(s1-x(j,0), 2) + pow(s2-x(j,1), 2);
      if(distSq>maxDistSq)
	continue; 
      lam += lam0*exp(-distSq/(2*kappa*kappa))*trapOper(j);
    }
    // Pr(y=0|s). Same as dpois(0, lam).
    qstar += exp(-lam);
  }
  // double qstar=mean(qstar_ind);
  return qstar/nsim;
}




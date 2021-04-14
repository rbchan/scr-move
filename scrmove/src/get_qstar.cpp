#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


/*
  Pr(not detected) computed using Monte Carlo integration, by 
  averaged over possible paths. 
*/


// [[Rcpp::export]]
double get_qstar_cam(int nsim, arma::mat x, int K,
		     double sigma, double rho, double lam0, double kappa,
		     double buffer, double maxDistSq, arma::imat notOper,
		     int nthreads) {

  double xmin=x.col(0).min()-buffer;
  double xmax=x.col(0).max()+buffer;
  double xdiff=xmax-xmin;
  double ymin=x.col(1).min()-buffer;
  double ymax=x.col(1).max()+buffer;
  double ydiff=ymax-ymin;
  int J=x.n_rows;
  // arma::cube u = arma::zeros<arma::cube>(nsim,K,2);
  // arma::vec qstar_ind = arma::zeros<arma::vec>(nsim);
  double qstar=0.0;
  double tau = std::sqrt(sigma*sigma-sigma*sigma*rho*rho);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for reduction(+:qstar)
  for(int i=0; i<nsim; i++) {
    // double s1 = R::runif(xmin, xmax);
    // double s2 = R::runif(ymin, ymax);
    double s1 = arma::randu<double>()*xdiff+xmin;
    double s2 = arma::randu<double>()*ydiff+ymin;
    double distSq = 0.0;
    double lam = 0.0;
    // u(i,0,0) = R::rnorm(s1, sigma);
    // u(i,0,1) = R::rnorm(s2, sigma);
    double u1 = arma::randn<double>()*sigma+s1;
    double u2 = arma::randn<double>()*sigma+s2;
    for(int j=0; j<J; j++) {
      if(notOper(j,0))
	continue;
      // distSq = pow(u(i,0,0)-x(j,0), 2) + pow(u(i,0,1)-x(j,1), 2);
      distSq = std::pow(u1-x(j,0), 2) + std::pow(u2-x(j,1), 2);
      if(distSq>maxDistSq)
	continue; 
      lam += lam0*std::exp(-distSq/(2*kappa*kappa));
    }
    for(int k=1; k<K; k++) { 
      // u(i,k,0) = R::rnorm(s1+rho*(u(i,k-1,0)-s1), tau);
      // u(i,k,1) = R::rnorm(s2+rho*(u(i,k-1,1)-s2), tau);
      u1 = arma::randn<double>()*tau+(s1+rho*(u1-s1));
      u2 = arma::randn<double>()*tau+(s2+rho*(u2-s2));
      for(int j=0; j<J; j++) {
	if(notOper(j,k))
	  continue;
	// distSq = pow(u(i,k,0)-x(j,0),2) + pow(u(i,k,1)-x(j,1),2);
	distSq = std::pow(u1-x(j,0),2) + std::pow(u2-x(j,1),2);
	if(distSq>maxDistSq)
	  continue; 
	lam += lam0*std::exp(-distSq/(2*kappa*kappa));
      }
    }
    // Pr(y=0|s,u). Same as dpois(0, lam).
    qstar += std::exp(-lam);
  }

  double qstar_bar = qstar/nsim;
  return qstar_bar;
}




// Pr(y=0|s)=\int_S p(y=0|u)p(u|s)du
// Condition on known value of activity center (s)
// There might be a closed-form solution here, as it's a convolution of
// two Gaussian kernels, but I'm be lazy and using MC integration

// [[Rcpp::export]]
double get_ld_ydet_aug(int nsim, double s1, double s2, arma::mat x, int K,
		       double sigma, double rho, double lam0, double kappa,
		       double maxDistSq, arma::imat notOper, int nthreads) {

  int J=x.n_rows;
  double qstar=0.0;
  double tau = std::sqrt(sigma*sigma-sigma*sigma*rho*rho);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for reduction(+:qstar)
  for(int i=0; i<nsim; i++) {
    double distSq = 0.0;
    double lam = 0.0;
    double u1 = arma::randn<double>()*sigma+s1;
    double u2 = arma::randn<double>()*sigma+s2;
    for(int j=0; j<J; j++) {
      if(notOper(j,0))
	continue;
      distSq = std::pow(u1-x(j,0), 2) + std::pow(u2-x(j,1), 2);
      if(distSq>maxDistSq)
	continue; 
      lam += lam0*std::exp(-distSq/(2*kappa*kappa));
    }
    for(int k=1; k<K; k++) { 
      u1 = arma::randn<double>()*tau+(s1+rho*(u1-s1));
      u2 = arma::randn<double>()*tau+(s2+rho*(u2-s2));
      for(int j=0; j<J; j++) {
	if(notOper(j,k))
	  continue;
	distSq = std::pow(u1-x(j,0),2) + std::pow(u2-x(j,1),2);
	if(distSq>maxDistSq)
	  continue; 
	lam += lam0*std::exp(-distSq/(2*kappa*kappa));
      }
    }
    // Pr(y=0|s,u). Same as dpois(0, lam).
    qstar += std::exp(-lam);
  }
  
  double ld_ydet = std::log(qstar/nsim);
  // shouldn't be necessary, but was getting NaN on rare occasions
  // that seemed to only occur (again, rarely) when qstar=0.0
  if(std::isnan(ld_ydet)) {
    Rprintf("NaN replaced with %f in get_ld_ydet_aug \n", 0.0);
    ld_ydet=0.0;
  }
  return ld_ydet;
}





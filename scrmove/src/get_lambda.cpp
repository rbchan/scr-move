#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
arma::cube get_lambda(double lam0, arma::cube distSq, double kappa,
		      arma::icube farOut, arma::imat notOper,
		      int nthreads) {
  int M=distSq.n_rows;
  int J=distSq.n_cols;
  int K=distSq.n_slices;
  arma::cube lambda = arma::zeros<arma::cube>(M, J, K);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(3) // shared(lambda,lam0,farOut,notOper,distSq,kappa)
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++) {
      for(int k=0; k<K; k++) {
	if(farOut(i,j,k) || notOper(j,k))
	  continue;
	lambda(i,j,k) = lam0*std::exp(-distSq(i,j,k)/(2*kappa*kappa));
      }
    }
  }

  return lambda;
}



// [[Rcpp::export]]
arma::cube get_lambda_nou(double lam0, arma::mat distSq, double sigma,
			  arma::imat farOut, arma::imat notOper,
			  int nthreads) {
  int M=distSq.n_rows;
  int J=distSq.n_cols;
  int K=notOper.n_cols;
  arma::cube lambda = arma::zeros<arma::cube>(M, J, K);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(3) // shared(lambda,lam0,farOut,notOper,distSq,kappa)
  for(int i=0; i<M; i++) {
    for(int j=0; j<J; j++) {
      if(farOut(i,j))
	continue;
      double tmp = lam0*std::exp(-distSq(i,j)/(2*sigma*sigma));
      for(int k=0; k<K; k++) {
	if(notOper(j,k))
	  continue;
	lambda(i,j,k) = tmp;
      }
    }
  }

  return lambda;
}



// [[Rcpp::export]]
arma::mat get_lambda_noK(double lam0, arma::mat distSq, double kappa,
			 arma::mat farOut, arma::ivec trapOper,
			 int nthreads) {
  int n=distSq.n_rows;
  int J=distSq.n_cols;
  arma::mat lambda = arma::zeros<arma::mat>(n, J);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(2) // shared(lambda,lam0,distSq,kappa,trapOper)
  for(int i=0; i<n; i++) {
    for(int j=0; j<J; j++) {
      lambda(i,j) = lam0*std::exp(-distSq(i,j)/(2*kappa*kappa))*trapOper(j);
    }
  }

  return lambda;
}


// [[Rcpp::export]]
arma::vec get_ld_y_n(arma::icube y, arma::cube lambda, int nthreads) {
  int n=y.n_rows;
  int J=y.n_cols;
  int K=y.n_slices;
  arma::vec ld_y = arma::zeros<arma::vec>(n);
  arma::vec Lambda = arma::zeros<arma::vec>(n);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(3) // shared(y,lambda) //reduction(+:ld_y)
  for(int i=0; i<n; i++) {
    bool goodlam = true;
    for(int j=0; j<J; j++) {
      if(!goodlam)
	continue;
      for(int k=0; k<K; k++) {
	if(lambda(i,j,k)>0.0) { // Avoid R API when result of dpois isn't finite
	  // may need to switch to stats::dpois
	  ld_y(i) += R::dpois(y(i,j,k), lambda(i,j,k), 1);
	} else if(y(i,j,k)>0) {
	  ld_y(i) = -std::numeric_limits<double>::infinity();
	  goodlam = false;
	  break;
	}
      }
    }
  }

  return ld_y;
}


// [[Rcpp::export]]
arma::mat get_ld_y_nK(arma::icube y, arma::cube lambda, int nthreads) {
  // Don't try to set ld_y to zero when dist>trim
  // Would have to allow for -Inf when y>0
  int n=y.n_rows;
  int J=y.n_cols;
  int K=y.n_slices;
  arma::mat ld_y = arma::zeros<arma::mat>(n,K);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(3) // shared(y,lambda) // reduction(+:ld_y)
  for(int i=0; i<n; i++) {
    for(int k=0; k<K; k++) {
      for(int j=0; j<J; j++) {
	if(lambda(i,j,k)>0.0) { // Avoid R API when result of dpois isn't finite
	  // may need to switch to stats::dpois
	  ld_y(i,k) += R::dpois(y(i,j,k), lambda(i,j,k), 1);
	} else if(y(i,j,k)>0) {
	  ld_y(i,k) = -std::numeric_limits<double>::infinity();
	  break;
	}
      }
    }
  }

  return ld_y;
}




// [[Rcpp::export]]
arma::vec get_ld_y(arma::icube y, arma::cube lambda,
		   arma::ivec z, int M) {
  int n=y.n_rows;
  int J=y.n_cols;
  int K=y.n_slices;
  arma::vec ld_y = arma::zeros<arma::vec>(M);
  arma::vec Lambda = arma::zeros<arma::vec>(M);
  for(int i=0; i<n; i++) {
    bool goodlam=true;
    for(int j=0; j<J; j++) {
      if(!goodlam)
	continue;
      for(int k=0; k<K; k++) {
	if(lambda(i,j,k)>0.0) { // Avoid R API when result of dpois isn't finite
	  // may need to switch to stats::dpois
	  ld_y(i) += R::dpois(y(i,j,k), lambda(i,j,k), 1);
	} else if(y(i,j,k)>0) {
	  ld_y(i) = -std::numeric_limits<double>::infinity();
	  goodlam = false;
	  break;
	}
      }
    }
  }
  for(int i=n; i<M; i++) {
    if(z(i)==1) {
      for(int j=0; j<J; j++) {
	for(int k=0; k<K; k++) {
	  Lambda(i) += lambda(i,j,k);
	}
      }
      ld_y(i) = -Lambda(i); // same as log(dpois(0,Lambda*K))
    } else {
      ld_y(i) = 0;
    }
  }
  return ld_y;
}



// [[Rcpp::export]]
arma::vec get_ld_y_noK(arma::imat y, arma::mat lambda, int nthreads) {
  int n=y.n_rows;
  int J=y.n_cols;
  arma::vec ld_y = arma::zeros<arma::vec>(n);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(2) // shared(y,lambda) // reduction(+:ld_y)
  for(int i=0; i<n; i++) {
    for(int j=0; j<J; j++) {
      if(lambda(i,j)>0.0) { // Avoid R API when result of dpois isn't finite
	ld_y(i) += R::dpois(y(i,j), lambda(i,j), 1);
      } else if(y(i,j)>0) {
	ld_y(i) = -std::numeric_limits<double>::infinity();
	break;
      }
    }
  }

  return ld_y;
}





// [[Rcpp::export]]
arma::cube get_distSq_ux(arma::cube u, arma::mat x,
			 int nthreads) {
  int n=u.n_rows;
  int J=x.n_rows;
  int K=u.n_cols;
  arma::cube distSq = arma::zeros<arma::cube>(n, J, K);

  omp_set_num_threads(nthreads);
  #pragma omp parallel for //collapse(3) // shared(u,x)
  for(int i=0; i<n; i++) {
    for(int j=0; j<J; j++) {
      for(int k=0; k<K; k++) {
	distSq(i,j,k) = std::pow(u(i,k,0)-x(j,0),2) + std::pow(u(i,k,1)-x(j,1),2);
      }
    }
  }

  return distSq;
}

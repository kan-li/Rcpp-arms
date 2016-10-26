# include <RcppArmadillo.h>
// [[Rcpp::depends( RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat inner(arma::vec x,
             arma::vec y){
  arma::mat ip=x.t()*y;
  return(ip);
}

double sq(double x){
  return(x*x);
}

// [[Rcpp::export]]
double sumsq_parallel(vec x, int ncores)
{
  double sum = 0.0;
  omp_set_num_threads(ncores);
  #pragma omp parallel for shared(x) reduction(+:sum)
  for (int i=0; i<x.size(); i++){
	// cout << i << ", "omp_get_thread_num() << " of " << omp_get_num_threads() << endl;
    sum += sq(x(i));
  }
  
  return sum;
}

/*** R
vecA = rnorm(100000)
*/
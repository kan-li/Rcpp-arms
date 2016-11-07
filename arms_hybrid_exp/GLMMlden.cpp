# include <RcppArmadillo.h>

// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
double logden_beta(vec x, vec myu_long, mat design_matrix, vec y){
  vec temp = design_matrix*x + myu_long;
  vec one; one.ones(size(temp));
  vec logden = y%temp - log(one + exp(temp));
  return(sum(logden));
}

// [[Rcpp::export]]
double logden_u(double x, vec mybeta, vec myy, mat mycovariate, double mytau){
  vec x_long(size(myy)); x_long.fill(x);
  vec temp = mycovariate*mybeta + x_long;
  vec one; one.ones(size(temp));
  vec logden = myy%temp - log(one + exp(temp));
  return(sum(logden)-0.5*pow(x,2)*mytau);
}
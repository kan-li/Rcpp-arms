# include <RcppArmadillo.h>

// [[Rcpp::depends( RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

extern "C" {
  #include "arms.c"
}

using namespace Rcpp;
using namespace arma;
using namespace std;


struct logden_beta_parm {
  NumericVector myu;
  NumericVector nobs;
  mat design_matrix;
  vec y;
  vec beta;
  vec dir;
};

struct logden_u_parm {
  vec mybeta;
  vec myy;
  mat mycovariate;
  double mytau;
};

NumericVector rept(NumericVector x, NumericVector y) {
  int n = y.size();
  NumericVector myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
    ind += p;
  }
  return myvector;
}


double logden_beta(double x, void *beta_data){
  logden_beta_parm *d;
  d = static_cast<logden_beta_parm *> (beta_data);
  vec new_beta = d->beta + x*d->dir;
  vec myu_long = as<vec>(rept(d->myu, d->nobs));
  vec temp = d->design_matrix*new_beta + myu_long;
  vec one; one.ones(size(temp));
  vec logden = d->y%temp - log(one + exp(temp));
  return(sum(logden));
}


double logden_u(double x, void *u_data){
  logden_u_parm *d;
  d = static_cast<logden_u_parm *> (u_data);
  vec x_long(size(d->myy)); x_long.fill(x);
  vec temp = d->mycovariate*d->mybeta + x_long;
  vec one; one.ones(size(temp));
  vec logden = d->myy%temp - log(one + exp(temp));
  return(sum(logden)-0.5*pow(x,2)*d->mytau);
}

/* *********************************************************************** */

// [[Rcpp::export]]
List sampling(NumericVector nobs, mat design_matrix, vec y, double myshape, vec temp_beta, vec temp_u, double temp_tau, int n_smp ){

  mat beta_smp(temp_beta.n_elem, n_smp);
  mat u_smp(temp_u.n_elem,n_smp);
  vec tau_smp(n_smp);
  
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval,i, it;
  double xinit[10]={0.0,3.0,17.0,20.0}, xl1 = -100.0, xr1 = 100.0, xl2=-30, xr2=30;
  double xcent[10], qcent[10] = {5., 30., 70., 95.};
  unsigned seed = 44;
  double convex = 1.;
  int dometrop = 0;
  double betasamp[100];
  double usamp[100];
  

  /* initialise random number generator */
  srand(seed);
  
  double zero = 0;
  int I = temp_u.n_elem;

  for(it=0; it<n_smp;it++){
// cout << it << endl;
    
    /* sample beta*/
    /*-----------------------------------------------------------------------*/
    /* set up structures for each density function */
      logden_beta_parm beta_data;
      
      vec dir = as<vec>(rnorm(2,0,1));
      /* initialise data needed by normal density function */
      beta_data.myu = temp_u;
      beta_data.nobs = nobs;
      beta_data.design_matrix = design_matrix;
      beta_data.y = y;
      beta_data.beta = temp_beta;
      beta_data.dir = dir;
      
      err = arms(xinit,ninit,&xl2,&xr2,logden_beta,&beta_data,&convex,
                 npoint,dometrop,&zero,betasamp,nsamp,qcent,xcent,ncent,&neval);
      
      if(err>0){
        cout << err << endl;
        exit(1);
      }
      
      temp_beta = temp_beta+betasamp[0]*dir;
      beta_smp.col(it) = temp_beta;
    
    
    /* sample tau*/
    /*-----------------------------------------------------------------------*/
    double scale = 1/(0.5*sum(square(temp_u)));
    temp_tau = rgamma(1,myshape,scale)[0];
    tau_smp(it) = temp_tau;
    
    
    /* sample u*/
    /*-----------------------------------------------------------------------*/
    int counter = 0;
    int start = 0;
    int end = 0;
    
    for(i=0;i<I;i++){
      
      logden_u_parm u_data;
      
      start = counter + 1;
      end = counter + nobs(i);
//       cout << i << endl;
//       cout << "start "<< start << " end "<< end << endl;
      
      u_data.mybeta = temp_beta;
      u_data.myy = y.subvec((start-1),(end-1));
      u_data.mycovariate = design_matrix.rows((start-1),(end-1));
      u_data.mytau = temp_tau;
      
//       u_data.myy.print();
//       cout << endl;
//       u_data.mycovariate.print();
//       cout << endl;
    
      err = arms(xinit,ninit,&xl1,&xr1,logden_u,&u_data,&convex,
                 npoint,dometrop,&temp_u(i),usamp,nsamp,qcent,xcent,ncent,&neval);
      if(err>0){
        cout << err << endl;
        exit(1);
      }
      
      counter=end;
      
      temp_u(i) = u_smp(i, it) = usamp[0];
//    cout << "sample: " << temp_u(i) << ", " << u_smp(i, it) << endl;
      /*-----------------------------------------------------------------------*/
    }
  }
  
  List ret;
  ret["beta"] = beta_smp;
  ret["u"] = u_smp;
  ret["tau"] = tau_smp;
  
  return (ret);
  
}

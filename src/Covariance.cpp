// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Covariance
//' 
//' Calculates the correlation between two vectors.
//' 
//' @param A NxP matrix.
//' @param B NxQ matrix.
//' @param corMat Return correlation matrix?
//' @export
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP matCov(const arma::mat A, const arma::mat B, const bool corMat=false){
  // Dimensions
  const int n = A.n_rows;
  const int p = A.n_cols;
  const int q = B.n_cols;
  // Initialize
  const arma::colvec u = arma::ones(n);
  arma::mat Za(n,p);
  Za.zeros();
  arma::mat Zb(n,q);
  Zb.zeros();
  arma::mat R(p,q);
  R.zeros();
  // ** Calculation
  // Center
  Za = Za-arma::mean(Za,0);
  Zb = Zb-arma::mean(Zb,0);
  if(corMat){
    Za = arma::normalise(Za);
    Zb = arma::normalise(Zb);
  }
  // Inner product
  if(corMat){
    R = Za.t()*Zb;
  } else {
    R = Za.t()*Zb/(n-1);
  }
  // Output
  return Rcpp::wrap(R);
} 
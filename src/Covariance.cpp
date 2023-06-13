// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout <<  << std::endl; 

// Validated: 20-04-28

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
SEXP MatCov(const arma::mat A, const arma::mat B, const bool corMat=false){
  // Dimensions.
  const int n = A.n_rows;
  const int p = A.n_cols;
  const int q = B.n_cols;
  // Initialize.
  const arma::colvec u = arma::ones(n);
  arma::mat Za = arma::zeros(n, p);
  arma::rowvec ma = arma::zeros(1, p);
  arma::mat Ma = arma::zeros(n, p);
  arma::mat Zb = arma::zeros(n, q);
  arma::rowvec mb = arma::zeros(1, q);
  arma::mat Mb = arma::zeros(n, q);
  arma::mat R = arma::zeros(p, q);
  // ** Calculation
  // Column-wise means.
  ma = arma::mean(A, 0);
  mb = arma::mean(B, 0);
  // Generate centering matrices.
  for(int i=0; i<n; i++){
  	Ma.row(i) = ma;
  	Mb.row(i) = mb;
  }
  // Centered matrices.
  Za = A - Ma;
  Zb = B - Mb;
  // Scale, if correlation matrix. 
  if(corMat){
    Za = arma::normalise(Za);
    Zb = arma::normalise(Zb);
  }
  // // Inner product
  if(corMat){
    R = Za.t() * Zb;
  } else {
    R = (Za.t() * Zb) / (n-1);
  }
  // Output
  return Rcpp::wrap(R);
} 


//' Eigenvalues of Symmetric Matrix. 
//' 
//' Calculates the eigenvalues of a symmetric matrix. 
//' 
//' @param A symmetric matrix. 
//' 
//' @export
//' 
//' @return Numeric vector.
// [[Rcpp::export]]

SEXP EigSym(const arma::mat A) {
  
  arma::vec e;
  e = arma::eig_sym(A);
  
  // Output
  return Rcpp::wrap(e);
} 
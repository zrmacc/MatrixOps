// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Lp Norm
//'
//' Calculates the Lp norm of a vector
//' 
//' @param x Jx1 numeric vector.
//' @param p scalar power. 
//'
//' @return Scalar norm.
//' @export 
// [[Rcpp::export]]
SEXP Norm(const arma::vec x, const int p=2){
	// Norm
	const double L = arma::norm(x,p);
	// Export
	return Rcpp::wrap(L);
}
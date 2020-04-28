// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Validated: 19-11-10

//' Log Sum Exp
//'
//' Calculates the Jx1 log sum exp vector. 
//' 
//' @param x Jx1 numeric vector.
//' @param cum Return the cumulative sum vector?
//'
//' @return Jx1 cumulative log sum exp vector. 
//' @export 
// [[Rcpp::export]]
SEXP logSumExp(const arma::vec x, const bool cum=false){
	// Dimensions
	const int J = x.n_elem;
	const arma::vec y = arma::sort(x,"descend");
	// Base vector
	arma::vec S = y[0]*arma::ones(J);
	for(int j=0;j<(J-1);j++){
		S(j+1)=fmax(y(j+1),S(j))+log1p(exp(-1*abs(y(j+1)-S(j))));
	}
	// Export
	if(cum){
		return Rcpp::wrap(S);
	} else {
		return Rcpp::wrap(S(J-1));
	}
}
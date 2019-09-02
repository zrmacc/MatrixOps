// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Projection Decomposition
//' 
//' Decomposes matrix \eqn{Y} into the projection onto the image of \eqn{X},
//' and the projection onto the orthogonal complement of the image. 
//' 
//' @param X NxP Numeric matrix.
//' @param Y NxQ Numeric matrix.
//' @return List containing the following:
//' \item{Coord}{Coordinates of the projection w.r.t. X.}
//' \item{Para}{Projection onto the image of X.}
//' \item{Ortho}{Projection onto the orthogonal complement.}
// [[Rcpp::export]]

SEXP projDecomp(const arma::mat X, const arma::mat Y){
  const arma::mat B = arma::solve(X.t()*X,X.t()*Y,arma::solve_opts::likely_sympd);
  const arma::mat P = X*B;
  const arma::mat Q = Y-P;
  return Rcpp::List::create(Rcpp::Named("Coord")=B,Rcpp::Named("Para")=P,Rcpp::Named("Ortho")=Q);
}

//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y Nx1 Numeric vector.
//' @param X NxP Numeric matrix.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' @export
// [[Rcpp::export]]
SEXP fitOLS(const arma::colvec y, const arma::mat X){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.n_cols;
  // Information
  const arma::mat A = X.t()*X;
  // Estimate beta
  const arma::vec b = arma::solve(A,X.t()*y,arma::solve_opts::likely_sympd);
  // Calculate residuals
  const arma::vec eps = (y-X*b);
  // Scale
  const double v = arma::as_scalar(eps.t()*eps/(n-p));
  // Information
  const arma::mat Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}

//' Weighted Least Squares
//' 
//' Fits the following weighted least squares model: 
//' \eqn{y_{i}=x_{i}'\beta+\epsilon_{i}}. Here, the subject-specific residual is
//' normally distributed with mean zero and variance \eqn{\sigma^{2}/w_{i}}.
//' \eqn{w_{i}} is a known, subject-specific weight, and \eqn{\sigma} is a
//' common scale parameter.
//' 
//' @param y Nx1 Response vector.
//' @param X NxP Design matrix.
//' @param w Nx1 Weight vector.
//' @export
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//'
// [[Rcpp::export]]
SEXP fitWLS(const arma::vec y, const arma::mat X, const arma::vec w){
  // Dimensions
  const int n = X.n_rows;
  const int p = X.n_cols;
  // Initialize
  double qf = 0;
  arma::mat A(p,p);
  A.zeros();
  arma::vec u(p);
  u.zeros();
  // ** Calculation
  // A = X'WX;
  for(int i=0;i<n;i++){
    A += w(i)*X.row(i)*X.row(i).t();
  }
  // u = X'Wy;
  for(int i=0;i<n;i++){
    u += X.row(i).t()*w(i)*y(i);
  }
  // Estimate beta
  const arma::vec b = arma::solve(A,u,arma::solve_opts::likely_sympd);
  // Calculate residuals
  const arma::vec eps = (y-X*b);
  // Scale
  for(int i=0;i<n;i++){
    qf += eps(i)*w(i)*w(i)*eps(i);
  }
  const double v = qf/(n-p);
  // Information
  const arma::mat Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Linear Projection
//' 
//' Decomposes matrix \eqn{Y} into the projection onto the image of \eqn{X},
//' and the projection onto the orthogonal complement of the image. 
//' 
//' @param X Numeric matrix.
//' @param Y Numeric matrix.
//' @export
//' 
//' @return List containing the following:
//' \item{Coord}{Coordinates of the projection w.r.t. X.}
//' \item{Para}{Projection onto the image of X.}
//' \item{Ortho}{Projection onto the orthogonal complement.}
// [[Rcpp::export]]

SEXP linProj(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Y){
  const Eigen::MatrixXd B=(X.transpose()*X).ldlt().solve(X.transpose()*Y);
  const Eigen::MatrixXd P=X*B;
  const Eigen::MatrixXd Q=Y-P;
  return Rcpp::List::create(Rcpp::Named("Coord")=B,Rcpp::Named("Para")=P,Rcpp::Named("Ortho")=Q);
}

//' Univariate OLS model.
//' 
//' Fits the standard OLS model.
//' 
//' @param y Outcome.
//' @param X Model matrix.
//' @export 
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' 
// [[Rcpp::export]]

SEXP fitOLS(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.cols();
  // Gram matrix
  const Eigen::MatrixXd XtX = X.transpose()*X;
  // Estimate beta
  const Eigen::VectorXd b = (XtX).ldlt().solve(X.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eps = (y-X*b);
  // Scale
  const double qf = (eps.transpose()*eps);
  const double v = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = XtX/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}
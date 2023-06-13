// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y Nx1 Numeric vector.
//' @param X NxP Numeric matrix.
//' 
//' @return List of model components.
//' @export
// [[Rcpp::export]]
SEXP FitOLS(
    const arma::colvec y,
    const arma::mat X
){
  // Observations.
  const int n = y.size();
  
  // Estimated parameters.
  const int p = X.n_cols;
  
  // Information.
  const arma::mat A = X.t() * X;
  
  // Estimate beta.
  const arma::vec b = arma::solve(A, X.t() * y, arma::solve_opts::likely_sympd);
  
  // Calculate residuals.
  const arma::vec yhat = X * b;
  const arma::vec eps = (y - yhat);
  
  // Scale.
  const double v = arma::as_scalar((eps.t() * eps) / (n-p));
  
  // Information.
  const arma::mat Ibb = A / v;
  const arma::mat invIbb = arma::pinv(Ibb);
  const arma::vec se = arma::sqrt(invIbb.diag());

  return Rcpp::List::create(
    Rcpp::Named("beta") = b,
    Rcpp::Named("info") = Ibb,
    Rcpp::Named("resid") = eps,
    Rcpp::Named("se") = se,
    Rcpp::Named("sigma2") = v,
    Rcpp::Named("yhat") = yhat
  );
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
//' 
//' @return List of model components.
//' @export
// [[Rcpp::export]]
SEXP FitWLS(const arma::vec y, const arma::mat X, const arma::vec w){
  // Dimensions
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // Initialize
  double qf = 0;
  arma::mat A = arma::zeros(p, p);
  arma::vec u = arma::zeros(p);
  
  // Calculation
  // A = X'WX;
  for(int i=0; i<n; i++){
    A += w(i) * X.row(i).t() * X.row(i);
  }
  
  // u = X'Wy;
  for(int i=0; i<n; i++){
    u += X.row(i).t() * w(i) * y(i);
  }

  // Estimate beta
  const arma::vec b = arma::solve(A, u, arma::solve_opts::likely_sympd);
  
  // Calculate residuals.
  const arma::vec yhat = X * b;
  const arma::vec eps = (y - yhat);
  
  // Scale
  for(int i=0; i<n; i++){
    qf += eps(i) * w(i) * eps(i);
  }
  const double v = qf / (n-p);
  
  // Information
  const arma::mat Ibb = A / v;
  const arma::mat invIbb = arma::pinv(Ibb);
  const arma::vec se = arma::sqrt(invIbb.diag());

  return Rcpp::List::create(
    Rcpp::Named("beta") = b, 
    Rcpp::Named("info") = Ibb, 
    Rcpp::Named("resid") = eps,
    Rcpp::Named("se") = se,
    Rcpp::Named("sigma2") = v,
    Rcpp::Named("yhat") = yhat
  );
}
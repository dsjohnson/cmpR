#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;


double com_log_sum(double x, double y){
  if(x==R_NegInf){return y;}
  else if(y==R_NegInf){return x;}
  else if(x > y){
    return x + log(1 + exp(y - x));
  } else{
    return y + log(1 + exp(x - y));
  }
}

double compute_log_z(const double& mu, const double& nu,
                     const double& log_eps, const bool& approx) {
  double ln_z;
  if(mu<=0 | nu<0){ln_z = R_NaN;}
  else if(nu==0){
     if(mu<1) ln_z = -log(1-pow(mu,nu));
     else ln_z=R_NaN;
  }
  else{
    if(approx){
      ln_z = nu*mu - (nu-1)*log(mu)/(nu) - (nu-1)*log(2*M_PI)/2 - log(nu)/2;
    } else{
      ln_z = R_NegInf;
      double ln_z_last = 0;
      int j = 0;
      while(std::abs(ln_z - ln_z_last) > log_eps){
        ln_z_last = ln_z;
        ln_z = com_log_sum(ln_z, j*nu*log(mu) - nu*lgamma(j+1));
        j = j+1;
      }
    }
  }
  return ln_z;
}

double log_z(const double& mu, const double& nu){
  double out = compute_log_z(mu, nu, 1.0E-8, true);
  if(out<2000) out = compute_log_z(mu, nu, 1.0E-8, false);
  return out;
}

double cmp_log_pmf(const int& x, const double& mu, const double& nu){
  double out;
    if(x<0){out = R_NaN;}
    else{
      out = x*nu*log(mu) - nu*lgamma(x+1) - log_z(mu,nu);
    }
  return out;
}

// [[Rcpp::export]]
NumericVector ln_dcmp(const NumericVector& x, const NumericVector& mu, const NumericVector& nu){
  NumericVector out(x.size());
  for(int i=0; i<out.size(); i++){
    out(i) = cmp_log_pmf(x(i), mu(i), nu(i));
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector samp_cmp(const int& n, const double& mu, const double& nu){
  IntegerVector out(n);
  double U;
  double ln_cmf;
  double ln_cmf_last;
  double ln_pmf;
  int x;

  for(int i=0; i<n; i++){
    U = Rcpp::as<double>(Rcpp::runif(1));
    ln_cmf = R_NegInf;
    x= -1;
    while(ln_cmf<=log(U)){
      x += 1;
      ln_cmf = com_log_sum(ln_cmf, cmp_log_pmf(x, mu, nu));
    }
    out(i) = x;
  }
  return out;
}

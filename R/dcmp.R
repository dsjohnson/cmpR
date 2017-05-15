#' @title Conway-Maxwell-Poisson Distribution
#' @aliases rcmp
#' @rdname cmp
#' @description Probability mass function and random generation for the
#' Conway-Maxwell-Poisson distribution for given values of the parameters.
#' @param x quantile at which the prob. mass function is evaluated.
#' @param mu location parameter
#' @param nu dispersion parameter
#' @param log Logical. Whether or not the log PMF is returned.
#' @param n number of random draws to return
#'
#' @details Computes PMF and makes random draws from the Conway-Maxwell-Poisson (CMP)
#' distribution. The PMF of the CMP is given by
#' \deqn{f(x) = (1/Z) ((nu*mu)^x)/(x!^nu).}
#' (Guikema and Goffelt 2008).
#'
#' The normalizing constant is calculated using a combination of finite truncation of the infinite sum
#' as well as an approximation provided by Shmueli et al. (2005) for small nu and large mu.
#' @author Devin S. Johnson and Jeffrey Dunn
#' @references Guikema, Seth D., and Jeremy P. Goffelt. "A flexible count data regression model for risk analysis." Risk analysis 28.1 (2008): 213-223.
#'
#' Shmueli, Galit, et al. "A useful distribution for fitting discrete data: revival of the Conway–Maxwell–Poisson distribution." Journal of the Royal Statistical Society: Series C (Applied Statistics) 54.1 (2005): 127-142.
#' @export

dcmp = function(x, mu, nu, log=FALSE){
  if(length(mu)==1){
    mu = rep(mu, length(x))
  } else if(length(mu!=length(x))) stop("Length of mu vector must equal 1 or length(x)")
  if(length(nu)==1){
    nu = rep(nu, length(x))
  } else if(length(nu!=length(x))) stop("Length of nu vector must equal 1 or length(x)")
  ff = ln_dcmp(x, mu, nu)
  if(log) return (ff)
  else return(exp(ff))
}

#' @author Devin S. Johnson and Jeffrey Dunn
#' @export
rcmp = function(n, mu, nu) {
  if (mu < 0 || nu < 0) return(rep(NaN,n))
  else(return(samp_cmp(n, mu, nu)))
}

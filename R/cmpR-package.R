#' @title Functions for Conway-Maxwell-Poisson Distribution
#'
#' @description Distribution and sampling fucntion are provided for the Conway-Maxwell-Poisson distribution.
#' The 'mu' parameterization of Guikema and Goffelt (2008), where floor(mu) is the median of the distribution.
#' The other parameter 'nu' controls over and underdispersion. For 0 < nu < 1, the distribution
#' is overdispersed relative to a Poisson. For nu > 1, the distribution is underdispersed.
#' For nu = 0, the distribution is equal to a geometric (for mu<1) and if nu = 1,
#' the distribution is equal to the Poisson.
#'
#' \tabular{ll}{
#' Package: \tab cmpR\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.1\cr
#' Date: \tab May 15, 2017\cr
#' License: \tab CC0 \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @note This software package is developed and maintained by scientists at the NOAA Fisheries Alaska
#' Fisheries Science Center and should be considered a fundamental research communication.
#' The reccomendations and conclusions presented here are those of
#' the authors and this software should not be construed as official communication by NMFS, NOAA,
#' or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the
#' National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the
#' highest quality, tools such as this are under constant development and are subject to change.
#'
#' @name cmpR-package
#' @aliases cmp-package cmpR
#' @docType package
#' @author Devin S. Johnson
#'
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#' @references Guikema, S.D. and Goffelt, J.P., 2008. A flexible count data regression
#' model for risk analysis. Risk analysis, 28(1), pp.213-223.
#' @useDynLib cmpR
#' @importFrom Rcpp evalCpp
NULL

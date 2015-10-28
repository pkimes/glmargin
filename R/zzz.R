#' KPEGASOS Reference Class
#'
#' Exposed Rcpp module class
#'
#' @field b a numeric vector
#' @field W a numeric vector
#'
#' @section Methods:
#' \describe{
#' \item{\code{new(x, y, py, lambda, max_iter, min_eps, verbose):}}{constructor call with training dataset}
#' \item{\code{Solve():}}{solve the problem}
#' \item{\code{boundaries():}}{return interval boundaries}
#' }
#'
#' @name KPEGASOS
#' @rdname kpegasos-module
#' @author Patrick Kimes
loadModule("kpegasos", TRUE)


#' Rcpp KPEGASOS object
#'
#' S4 class for fitting KPEGASOS learning problem to training data
#'
#' @exportClass Rcpp_KPEGASOS
#' @name Rcpp_KPEGASOS-class
#' @rdname Rcpp_KPEGASOS-class
#' @aliases Rcpp_KPEGASOS-class
#' @author Patrick Kimes
NULL

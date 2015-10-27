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
#' @name kpegasos-module
#' @author Patrick Kimes
loadModule("kpegasos", TRUE)

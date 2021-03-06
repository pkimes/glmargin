% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/zzz.R
\name{KPEGASOS}
\alias{KPEGASOS}
\title{KPEGASOS Reference Class}
\description{
Exposed Rcpp module class
}
\section{Fields}{

\describe{
\item{\code{b}}{a numeric vector}

\item{\code{W}}{a numeric vector}
}}
\section{Methods}{

\describe{
\item{\code{new(x, y, py, lambda, max_iter, min_eps, verbose):}}{constructor call with training dataset}
\item{\code{Solve():}}{solve the problem}
\item{\code{boundaries():}}{return interval boundaries}
}
}
\author{
Patrick Kimes
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ars.R
\name{ars}
\alias{ars}
\title{Documentation for ars function}
\usage{
ars(g, n_samp, x_bound, x_start = c(), sybolic_deriv = FALSE)
}
\arguments{
\item{g}{(Required) Function containing the target probability density function (pdf).}

\item{n_samp}{(Required) An integer representing the desired sample size.}

\item{x_bound}{(Required) A vector containing: (1) the lower bound of the domain of \code{g},
and (2) the upper bound.}

\item{x_start}{(Optional) A vector containing the initial abscissae.
If left empty, two points will be chosen by the \code{ars()} function. Default: empty.}

\item{symbolic_deriv}{(Optional) A logical vector indicating whether to calculate
derivatives symbolically (\code{symbolic_deriv=TRUE}), or numerically (\code{symbolic_deriv=FALSE}).
Default: TRUE.}
}
\value{
A vector of \code{n_samp} accepted points based on \code{g}.
}
\description{
Returns sample points based on advanced rejection sampling technique
}
\examples{
Normal distribution (sample of 50 points, no x_start given):
g <- function(x,mu=0,sigma=1){
  val <- 1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
}
ars(g,n_samp=5,x_bound=c(-Inf,Inf))


Gamma distribution (sample of 100 points, no x_start given):
g <- function(x,alpha=2,beta=1){
  val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
}
ars(g,n_samp=100,x_bound=c(0,Inf))


Beta distribution (sample of 20 points, x_start given):
g <- function(x,alpha=2,beta=1){
  val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
}
ars(g,n_samp=20,x_bound=c(0,1),x_start=c(0.25,0.75))
}

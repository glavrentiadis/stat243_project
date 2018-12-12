#' Documentation for ars function
#'
#' Returns sample points based on advanced rejection sampling technique
#'
#' @param g (Required) Function containing the target probability density function (pdf).
#' @param n_samp (Required) An integer representing the desired sample size.
#' @param x_bound (Required) A vector containing: (1) the lower bound of the domain of \code{g},
#' and (2) the upper bound.
#' @param x_start (Optional) A vector containing the initial abscissae.
#' If left empty, two points will be chosen by the \code{ars()} function. Default: empty.
#' @param symbolic_deriv (Optional) A logical vector indicating whether to calculate
#' derivatives symbolically (\code{symbolic_deriv=TRUE}), or numerically (\code{symbolic_deriv=FALSE}).
#' Default: TRUE.
#' @return A vector of \code{n_samp} accepted points based on \code{g}.
#' @examples
#' Normal distribution (sample of 50 points, no x_start given):
#' g <- function(x,mu=0,sigma=1){
#'   val <- 1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
#' }
#' ars(g,n_samp=5,x_bound=c(-Inf,Inf))
#'
#'
#' Gamma distribution (sample of 100 points, no x_start given):
#' g <- function(x,alpha=2,beta=1){
#'   val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
#' }
#' ars(g,n_samp=100,x_bound=c(0,Inf))
#'
#'
#' Beta distribution (sample of 20 points, x_start given):
#' g <- function(x,alpha=2,beta=1){
#'   val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
#' }
#' ars(g,n_samp=20,x_bound=c(0,Inf),x_start=c(0.25,0.75))
#' @export
ars <- function(g,n_samp,x_bound,x_start = c(),sybolic_deriv = FALSE){

  #load libraries
  library(assertthat)
  library(testthat)
  library(truncnorm)
  #  library(Ryacas)
  library(Deriv)

  #Function handle for h
  eval_H <- EvalH(g)
  #Choose if numeric or symbolic derivative
  if (sybolic_deriv == FALSE){
    eval_H_prime <- NumericEvalHPrime(g)
  } else {
    eval_H_prime <- SymbolicEvalHPrime(g)
  }

  #initialization step
  #chech starting points if given by the user
  if (length(x_start) != 0){
    out <- CheckXStart(x_start = x_start, x_bound = x_bound, eval_h = eval_H, eval_deriv_h = eval_H_prime)
  } else {
    out <- SampleXStart(x_bound = x_bound, eval_h = eval_H, eval_deriv_h = eval_H_prime)
  }
  X <- out$X
  Z <- out$Z
  H <- out$H
  H_prime <- out$H_prime
  H_norm <- out$H_norm
  P_cum <- out$Pcum


  #Sampling and updating when needed
  x_accept <- numeric(n_samp)
  k <- 0
  while (k <= n_samp){
    #sample x_star
    x_star <- SamplePieceExp(X,Z,H_norm,H_prime,P_cum)
    #print(paste('Master: x_star:',x_star))
    #check if x_star will be kept and update X,Z,H,H_prime if needed
    out <- UpdateAccept(x_star = x_star, X = X ,Z = Z, H = H, H_prime = H_prime, H_norm = H_norm, P_cum = P_cum,
                        x_accept = x_accept, length_accept = k, eval_h = eval_H, eval_h_prime = eval_H_prime)
    X <- out$X
    Z <- out$Z
    H <- out$H
    H_prime <- out$H_prime
    H_norm <- out$H_norm
    P_cum <- out$P_cum
    x_accept <- out$x_accept
    k <- out$k
  }

  #return random samples
  return(x_accept)

}


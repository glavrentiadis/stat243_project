#Testing for auxilary functions


library(testthat)

#Generate inputs


context("Testing CalcProbBin auxiliary function with tolerance difference = 0.1")

test_that("CalcProbBin auxiliary function", {
  # set.seed(1)
  X <- rnorm(100)
  X <- sort(X)
  g <- function(x){
    sigma <- 1
    mu <- 0
    val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
  }
  
  library(Deriv)
  g_prime <- Deriv(g)
  update_H <- function(x) {log(g(x))}
  update_H_prime <- function(x) {g_prime(x)/g(x)}
  z <- function(i) {
    (H[i+1]-H[i]-X[i+1]*H_prime[i+1]+X[i]*H_prime[i])/(H_prime[i]-H_prime[i+1])
  }
  
  H <- sapply(X, update_H)
  H_prime <- sapply(X, update_H_prime)
  Z <- c()
  for (z_order in 1:length(X)-1) {
    Z[z_order] <- z(z_order)
  }
  Z <- append(Z,-Inf,0)
  Z <- append(Z,Inf,length(Z))
  
  Pcum <- CalcProbBin(X,Z,H,H_prime)$Pcum
  Pcum_check <- pnorm(Z[-1])
expect_equivalent(Pcum,Pcum_check,tolerance = 0.1)
})




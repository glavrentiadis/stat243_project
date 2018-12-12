#Testing for Ars and auxilary functions

library(testthat)


#testing

context("Testing ars function with Kolmogorov-Smirnov Test N(0,1)")
test_that("Kolmogorov-Smirnov Test (Standard normal distribution)", {
  #distribution to test
  g <- function(x){
    sigma <- 1
    mu <- 0
    val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
  }

  x_bound = c(-Inf,Inf)
  n_samp = 10000
  x_start = c(-5,5)
  check <- rnorm(n_samp)

  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)

  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value


  expect_gt(p_value,0.05)
})


context("Testing ars function with Kolmogorov-Smirnov Test N(7,2)")

test_that("Kolmogorov-Smirnov Test (Normal distribution, mu = 7, sd = 2)", {
  #distribution to test
  g <- function(x){
    sigma <- 2
    mu <- 7
    val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
  }

  x_bound = c(-Inf,Inf)
  n_samp = 10000
  x_start = c(-2,15)
  check <- rnorm(n_samp,mean = 7, sd = 2)

  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)

  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value

  expect_gt(p_value,0.05)
})


context("Testing ars function with Kolmogorov-Smirnov Test Beta(1,3)")

test_that("Kolmogorov-Smirnov Test (beta 1, 3)", {
  #distribution to test
  g <- function(x){
    alpha <- 1
    beta <- 3
    val <- x^(alpha-1)*(1-x)^(beta-1)/(gamma(alpha)*gamma(beta)/gamma(alpha+beta))
  }

  x_bound = c(0,1)
  n_samp = 10000
  x_start = c(0.1,0.6)
  check <- rbeta(n_samp,1,3)

  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)

  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value

  expect_gt(p_value,0.05)
})


context("Testing ars function with Kolmogorov-Smirnov Test Gamma(2,3)")

test_that("Kolmogorov-Smirnov Test (Gamma(2,3))", {
  #distribution to test
  g <- function(x){
    alpha <- 2
    beta <- 3
    val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
  }

  x_bound = c(0,Inf)
  n_samp = 10000
  x_start = c(1,5)
  check <- rgamma(n_samp,shape=2,rate=3)

  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)

  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value

  expect_gt(p_value,0.05)
})


context("Testing ars function with Kolmogorov-Smirnov Test Exponential(lambda = 5)")

test_that("Kolmogorov-Smirnov Test (Exponential(lambda = 5))", {
  #distribution to test
  lambda <- 5
  g <- function(x) {
    pi*lambda*exp(-lambda*x)
  }

  x_bound = c(0,Inf)
  n_samp = 10000
  x_start = c(0.5,3.5)
  check <- rexp(n_samp,rate = lambda)

  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)

  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value

  expect_gt(p_value,0.05)
})


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



#Testing for auxilary functions


library(testthat)


context("Testing UpdateAccept auxiliary function")
test_that("UpdateAccept auxiliary function", {
  #Generate inputs
  X = c(-2,0,2)
  H <- log(dnorm(X)) #h(x) = log(f(x))
  H_prime <- -1*X #h_prime(x) = -x
  Z <- c(-Inf,-1,1,Inf)
  P_cum <- c(1/6,1-1/6,1)
  H_norm <- c(0.04511176, 1/3, 0.04511176)
  x_star <- -2
  x_accept <- c()
  length_accept <- 1
  #standard normal distribution, h and h_prime functions
  eval_h <- function(x){
    return(log(dnorm(x)))
  }
  eval_h_prime <- function(x){
    #h = log(f) = -1/2*log(2*pi) - 1/2*(x)^2
    #h_prime = deriv(log(f)) = - x
    
    return(-1*x)
  }
  
  H_new <- UpdateAccept(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)$H
  H_prime_new <- UpdateAccept(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)$H_prime
  X_new <- UpdateAccept(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)$X
  x_accept_new <- UpdateAccept(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)$x_accept
  expect_equal(H,H_new)
  expect_equal(H_prime,H_prime_new)
  expect_equal(X,X_new)
  expect_equal(x_star,x_accept_new)
})




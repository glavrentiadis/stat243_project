#Testing for Ars and auxilary functions


library(testthat)

#load ars function
#function_file <- 'ars_version3.R'
#source(file.path('..',function_file))

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
  #x_start = c(0,10)
  check <- rexp(n_samp,rate = lambda)
  
  rand_samp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
  
  #perform Kolmogorov-Smirnov test
  p_value <- ks.test(rand_samp,check)$p.value
  
  expect_gt(p_value,0.05)
})



library(testthat)

##Test for N(0,1) distribution

g <- function(x){
  sigma <- 1
  mu <- 0
  val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
}

x_bound = c(-Inf,Inf)
n_samp = 10000
x_start = c(-5,5)
check <- rnorm(n_samp)

temp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
p_value <- ks.test(temp,check)$p.value

if (test_that("Success", {expect_gt(p_value,0.05)})==TRUE) 
  {
  print("Success. The ars function passes the Kolmogorov-Smirnov Test for N(0,1) function")
} else {
  print("Error")
}



##Test for N(7,2) distribution

g <- function(x){
  sigma <- 2
  mu <- 7
  val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
}

x_bound = c(-Inf,Inf)
n_samp = 10000
x_start = c(-2,15)
check <- rnorm(n_samp,mean = 7, sd = 2)

temp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
p_value <- ks.test(temp,check)$p.value

if (test_that("Success", {expect_gt(p_value,0.05)})==TRUE) 
{
  print("Success. The ars function passes the Kolmogorov-Smirnov Test for N(7,2) function")
} else {
  print("Error")
}


#Test for Beta(1,3) distribution
g <- function(x){
  alpha <- 1
  beta <- 3
  val <- x^(alpha-1)*(1-x)^(beta-1)/(gamma(alpha)*gamma(beta)/gamma(alpha+beta))
}

x_bound = c(0,1)
n_samp = 10000
x_start = c(0.1,0.6)
check <- rbeta(n_samp,1,3)

temp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
p_value <- ks.test(temp,check)$p.value

if (test_that("Success", {expect_gt(p_value,0.05)})==TRUE) 
{
  print("Success. The ars function passes the Kolmogorov-Smirnov Test for Beta(1,3) function")
} else {
  print("Error")
}




#Test for Gamma(2,3) distribution
g <- function(x){
  alpha <- 2 
  beta <- 3
  val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
}

x_bound = c(0,Inf)
n_samp = 10000
x_start = c(1,5)
check <- rgamma(n_samp,shape=2,rate=3)

temp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
p_value <- ks.test(temp,check)$p.value

if (test_that("Success", {expect_gt(p_value,0.05)})==TRUE) 
{
  print("Success. The ars function passes the Kolmogorov-Smirnov Test for Gamma(2,3) function")
} else {
  print("Error")
}


#Test for Exponential(lambda = 5) distribution
lambda <- 5
g <- function(x) {
  pi*lambda*exp(-lambda*x)
}

x_bound = c(0,Inf)
n_samp = 10000
#x_start = c(0,10)
check <- rexp(n_samp,rate = lambda)

temp <-ars(g,n_samp,x_bound,x_start,sybolic_deriv = FALSE)
p_value <- ks.test(temp,check)$p.value

if (test_that("Success", {expect_gt(p_value,0.05)})==TRUE) 
{
  print("Success. The ars function passes the Kolmogorov-Smirnov Test for Exponential(5) function")
} else {
  print("Error")
}


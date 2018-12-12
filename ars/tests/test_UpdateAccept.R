#Testing for auxilary functions


library(testthat)

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
context("Testing UpdateAccept auxiliary function")
test_that("UpdateAccept auxiliary function", {
  expect_equal(H,H_new)
  expect_equal(H_prime,H_prime_new)
  expect_equal(X,X_new)
  expect_equal(x_star,x_accept_new)
})




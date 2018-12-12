##Test CalcProbBin function to generate cumulative probabilities
#Generate inputs
X <- rnorm(10)
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
#Z <- append(Z,-Inf,0)
#Z <- append(Z,Inf,length(Z))

#CalcProbBin <- function(X,Z,H,H_prime)
  
  
return(list(Pcum = Pcum, log_s = log_s, Ptot = Ptot))



##Test UpdateAccept function

#Inputs
X = c(-2,0,2)
H <- log(dnorm(X)) #h(x) = log(f(x))
H_prime <- -1*X #h_prime(x) = -x
Z <- c(-Inf,-1,1,Inf)
P_cum <- c(1/6,1-1/6,1)
H_norm <- c(0.04511176, 1/3, 0.04511176)

UpdateAccept <- function(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)

#
return(list(X = X, Z = Z, H = H, H_prime = H_prime, H_norm = H_norm, P_cum = P_cum, x_accept = x_accept, k = length_accept))  
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
Z <- append(Z,-Inf,0)
Z <- append(Z,Inf,length(Z))

Pcum <- CalcProbBin(X,Z,H,H_prime)$Pcum
Pcum_check <- pnorm(Z[1])
all.equal(Pcum,Pcum_check,tolerance = 0.05)
  
  
return(list(Pcum = Pcum, log_s = log_s, Ptot = Ptot))



##Test UpdateAccept function

#Inputs
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

UpdateAccept(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime)

#It is true that x_star is added to x_accept, and no update for X, H, and H_prime
#####################################
#User supplied g-function

#Exponential
gExp <- function(x,lamb=1){
  value <- lam*exp(-lam*x)
}

#Normal
gN <- function(x,mu=0,sigma=1){
  val <- 3/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/(2*sigma^2))
}

#Gamma
gG <- function(x,alpha=2,beta=1){
  val <- beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
}

#Beta
gB <- function(x,alpha=1,beta=3){
  val <- x^(alpha-1)*(1-x)^(beta-1)/(gamma(alpha)*gamma(beta)/gamma(alpha+beta))
}

#####################################
#Other inputs
n <- 10000
x_bound <- c(-Inf,Inf)
#x_bound <- c(0,1)
g <- function(x){
  gN(x,mu=40,sigma=5) #Distribution of g
  #gB(x)
}
#####################################
#Run function

library(stringr)
run <- function(time=TRUE,saveOut=FALSE){
  if (time==TRUE) {
    print(system.time(
      ars(g,n,x_bound)
    ))
  }
  if (saveOut==TRUE){
    return(ars(g,n,x_bound))
  } 
}
output <- run(time==FALSE,saveOut==TRUE)
hist(output,breaks=25)
print(summary(output))
x <- 0.5
print(log(g(x)))
rm(output)

exp(710)

EvalH <- function(g){
  
  evalH <- function(x){
    return(log(g(x)))
  }
  
  return(evalH)
}
eval_H <- EvalH(g)

eval_H_prime <- function(x){
  #numerical step
  dx <- sqrt(.Machine$double.eps)*abs(max(x,1))
  #compute numerical derivative
  dh_dx <- (log(g(x+dx))-log(g(x-dx)))/(2*dx)
  
  return(dh_dx)    
}



X <- c(999,1001)
#X <- c(0.1,0.95)
H <- eval_H(X)
H_d <- eval_H_prime(X)

Z <- (H[-1] - H[-2] - X[-1]*H_d[-1] + X[-2]*H_d[-2])/(H_d[-2] - H_d[-1])
i_con_Hd <- abs(H_d[-1] - H_d[-2]) < .Machine$double.eps*abs(H_d[-2])
Z_med <- 0.5*(X[-1] + X[-2])
Z[i_con_Hd] <- Z_med[i_con_Hd]
#include the bounds in the array Z
Z <- c(x_bound[1],Z,x_bound[2])
out <- CalcProbBin(X,Z,H,H_d)
Pcum <- out[[1]]
H_norm <- out[[2]]

n_bin <- length(Z) -1
Z[-1]
exp(38.2)
exp(-3.85)/2.5*(exp(2.5*4.35)-exp(2.5*-Inf))
#calculate area under each bin (proportional to probabilty)
a_coef <- H - X*H_d #intersect of each line seg. (log space)
b_coef <- H_d #slope of each line seg. (log space)
P <- exp(a_coef)/b_coef*(exp(b_coef*Z[-1]) - exp(b_coef*Z[-(n_bin+1)]))
(exp(b_coef*Z[-1]) - exp(b_coef*Z[-(n_bin+1)]))
(exp(b_coef*Z[-1]))
exp(b_coef*Z[-(n_bin+1)])

#normalized probabilites
Ptot <- sum(P)
P <- P/Ptot
print(P)
exp(1000)

#temp <- data.frame()
#dat <- data.frame()
temp <- data.frame(X[1],X[2],H[1],H[2],H_d[1],H_d[2],Z[2],
                   a_coef[1],a_coef[2],b_coef[1],b_coef[2],P[1],P[2])
#dat <- temp
dat <- rbind(dat,temp)

print(dat)
#cumulative probabilites
Pcum <- cumsum(P) 
#normalized height of distributions
log_s <- H -log(Ptot)


exp(a_coef)/b_coef

exp(b_coef*Z[-1])-exp(b_coef*Z[-(n_bin+1)])  



SamplePieceExp(X,Z,H_norm,H_prime,Pcum)
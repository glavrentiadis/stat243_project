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
n <- 10
#x_bound <- c(-Inf,Inf)
x_bound <- c(-Inf,Inf)
#x_bound <- c(0,1)
g <- function(x){
  gN(x,mu=60,sigma=5) #Distribution of g
  #gB(x)
}
h <- function(x,orig){
  log(g(x))
}
h_d <- function(x){
  eval(Deriv(h(x),"x"))
}
print(eval(h_d(1,g,h)))

  evalH_prime <- function(x){
    #numerical step
    dx <- sqrt(.Machine$double.eps)*abs(max(x,1))
    #compute numerical derivative
    dh_dx <- (log(g(x+dx))-log(g(x-dx)))/(2*dx)
    
    return(dh_dx)    
  }
#####################################
#Run function
?optimize
?optim
#optimize(g,c(0,20),tol=1)
optimize(f=evalH_prime,interval=c(-.Machine$double.xmax,100),maximum=TRUE,tol=1
         #,control=list(fnscale=-1),
      #lower=-.Machine$double.xmax,upper=.Machine$double.max
      )
optim(par=-100,fn=evalH_prime,method="Brent",
        control=list(fnscale=-1),lower=-.Machine$double.xmax,upper=.Machine$double.max
)
log(h_d(-130.6))
evalH_prime(-140.6)
warnings()
prettyNum(-.Machine$double.xmax)
prettyNum(.Machine$double.xmax)
optimize(f=sin,interval=c(-.Machine$double.xmax,0))
optimize(f=Deriv(sin,"x"),interval=c(0,.Machine$double.xmax))
optim(fn=sin,par=-.Machine$double.xmax/10)
optim(fn=sin,par=.Machine$double.xmax/10)
.Machine$double.xmax/2-(-.Machine$double.xmax/2)
for (i in 1:100){
  
}

n <- 1000
  nums <- -.Machine$double.xmax+(.Machine$double.xmax/n)*seq(1:(n/2))*2
  #print(nums[50])
  nums[((n/2)+1):n] <- rev(nums)
  0
  which.max(nums)
  nums[53]
  nums[50]
  
  library(testthat)
  library(assertthat)
  install.packages("testit")
  library(testit)

  
  has_warning(optim(fn=evalH_prime,par=nums[1]))
  expect_error(optim(fn=evalH_prime,par=nums[1]))
  
  
  optim(fn=h,par=0,control=list(fnscale=-1))
  mean <- 50
  optim(fn=dnorm,mean=200,par=10,method="Brent",control=list(fnscale=-1),
        lower=-200,upper=250)
  ?dnorm


options(show.warning.messages = FALSE)
options(show.error.messages = FALSE)
options(error=NULL)
evalH_prime(-1000)
nums <-
#print(nums)
  #optim(fn=evalH_prime,par=seq[i])
  
  print(has_error(optim(fn=evalH_prime,par=nums[285])))
finite <- TRUE

find_Range <- function(seq){
  finite <- TRUE
  ct <- 1
  
  print(g(253))
  evalH_prime(251)
  print(all.equal(g(253),0))
  evalH_prime(1000000)
  while (finite==TRUE & ct <= 500){
    
    result <- has_error(optim(fn=evalH_prime,par=seq[ct]))
    if (result==TRUE){
      print("Heyo!")
      print(ct)
      finite <- FALSE
    }
    if (i%%50 == 0){
      print(result)
    }
    ct <- ct+1
  }
  #print(has_error(optim(fn=evalH_prime,par=nums[i])))
}

find_Range(nums)
warnings()
evalH_prime(252)

x <- 252
prettyNum(log(g(252+dx)),digits=8)
prettyNum(log(g(252-dx)),digits=8)
?deriv
x <- 1
print(deriv(dnorm))
log(g(252-dx))
dx <- sqrt(.Machine$double.eps)*abs(max(x,1))

#compute numerical derivative
dh_dx <- (log(g(x+dx))-log(g(x-dx)))/(2*dx)

print(optim(fn=evalH_prime,par=500))
options(show.warning.messages =TRUE)
options(show.error.messages = TRUE)
result <- has_error(optim(fn=evalH_prime,par=nums[1]))


-.Machine$double.xmax
sin(.Machine$double.xmax)
?sprintf
h(-12)
print(h(0))
print(h(60))
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
output <- run(time=FALSE,saveOut=TRUE)


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



X <- c(39.99,40.01)
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
pnorm(-7.5)

exp(a_coef)/b_coef

exp(b_coef*Z[-1])-exp(b_coef*Z[-(n_bin+1)])  



SamplePieceExp(X,Z,H_norm,H_prime,Pcum)

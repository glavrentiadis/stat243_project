
rm(list = ls())

ars <- function(g,n_samp,x.Bound,k_start = 2){
  
  #load libraries
  library(truncnorm)
  library(Ryacas)
  library(Deriv)
  
  #initialization step
  out <- Initial(k_start, x.Bound, g)
  X <- out[[1]]
  H <- out[[2]]
  H_prime <- out[[3]]
  Z <- out[[4]]
  H_norm <- out$H_norm
  Pcum <- out$Pcum
  Ptot <- out$Ptot

  rm(out)
  

  #Sampling and updating when needed
  x_accept <- numeric(n_samp)
  k <- 0
  while (k <= n_samp){
    #sample x_star
    x_star <- SamplePieceExp(X,Z,H_norm,H_prime,Pcum)
    #print(paste('Master: x_star:',x_star))
    #check if x_star will be kept and update X,Z,H,H_prime if needed
    out <- update_accept(x_star, g, X, H, H_prime, Z,k,x_accept,H_norm,Pcum) 
    X <- out[[1]]
    H <- out[[2]]
    H_prime <- out[[3]]
    Z <- out[[4]]
    x_accept <- out[[5]]
    k <- out[[6]]
    H_norm <- out$H_norm
    Pcum <- out$Pcum

    
  }
  
  #return random samples
  return(x_accept)

}




#Auxilary functions

#function to compute log of g (function h)
update_H <- function(x,g){
  log(g(x))
}

#function to compute the derivative of h
calc_g_prime <- function(X,g){
  return(1/g(X)*eval(Deriv(~g(X),"X")))
}




#Initial function
Initial <- function(k_start,x.Bound,g,QC=FALSE){

  upBd <- x.Bound[length(x.Bound)]
  lowBd <- x.Bound[1]
  if (lowBd != -Inf){
    needPosHprime <- FALSE
  } else{
    needPosHprime <- TRUE
  }
  if (upBd != Inf){
    needNegHprime <- FALSE
  } else{
    needNegHprime <- TRUE
  }
  upBdInit <- upBd #Upper bound for initialization points
  lowBdInit <- lowBd #Lower bound for initialization points
  if (0==1){
    
    if (QC==TRUE){
      #print("Lower & upper bound pairs")
      #print(lowBd)
      #print(upBd)
      
      #print(upBdInit)
      #print(lowBdInit)
    }
  }
  
  
  k <- 0
  X <- c()
  
  
  #print(paste0("Begin: Find initialization points"))
  while (k < k_start){
    samp <- rtruncnorm(1,a=lowBdInit,b=upBdInit)
    hPrime <- update_H_prime(samp,g)
    if (QC==TRUE){
      #print(paste0("Current X vector: ",X))
      #print(paste0("Initialization bounds: ",lowBdInit," ",upBdInit))
      #print(paste0("Need positive / negative h'? ",needPosHprime," ",needNegHprime))
      
      #print(paste0("Sampled point x: ",samp))
      #print(paste0("h'(x): ",update_H_prime(samp,g)))
    }
    
    if (hPrime >0 & needPosHprime==TRUE){
      needPosHprime <- FALSE
      X <- append(X,samp,0)
      k <- k+1
      #print(paste0("Appended X with h'(X) > 0"))
      
      if (needNegHprime==TRUE) {
        lowBdInit <- samp
        #print(paste0("Changed lower initialization bound to sampled point"))
      }
      else {
        lowBdInit <- lowBd
        #print(paste0("Changed lower initialization bound to distribution lower bound"))
      }
    }
    else if (hPrime < 0 & needNegHprime==TRUE) {
      needNegHprime <- FALSE
      X <- append(X,samp,length(X))
      k <- k+1
      #print(paste0("Appended X with h'(X) < 0"))
      
      if (needPosHprime == TRUE) {
        upBdInit <- samp
        #print(paste0("Changed upper initialization bound to sampled point"))
      }
      else {
        upBdInit <- upBd
        #print(paste0("Changed upper initialization bound to distribution upper bound"))
      }
    }
    else if (k < (k_start-needPosHprime-needNegHprime)){
      X <- append(X,samp,length(X[X<samp]))
      k <- k+1
      #print(paste0("Appended X; no h' requirements met"))
    }
    else {
      #print(paste0("Sampled point does not satisfy any initialization requirements. Sample again."))
    }
  }
  #print(paste0("Finished: Find initialization points"))
  H <- update_H(X,g)
  H_d <- update_H_prime(X,g)
  Z <- (H[-1] - H[-k_start] - X[-1]*H_d[-1] + X[-k_start]*H_d[-k_start])/(H_d[-k_start] - H_d[-1])

  #include the bounds in the array Z
  Z <- c(lowBd,Z,upBd)
  
  #compute the probabilites of each piece-wise exponential bin
  out <- CalcProbBin(X,Z,H,H_d)
  Pcum <- out[[1]]
  H_norm <- out[[2]]
  Ptot <- out[[3]]
  
  
  return(list(X = X,H = H,H_prime = H_d,Z = Z,Pcum = Pcum,H_norm = H_norm,Ptot = Ptot))
}

#Piece wise exponential sampling

SamplePieceExp <- function(X,Z,H_norm,H_prime,Pcum){
  #SamplePieceExp draws a sample out of a piece wise exponential distribution
  #Input:
  # X: numeric array with all the values where h and h' have been evaluated
  # Z: numeric array with the ordinates at the intersections of the piece-wise distribution
  # H_norm: numeric array with the normalized elevations of exponential distribution at locations X
  # H_prime:  numeric array with the evaluations of the first derivative of h at the locations X
  # Pcum: numeric array with the cumulative probability of each of each segment 
  #Output:
  # x_star: random sample of the distribution
  
  #sample inverse probability 
  Pinv <- runif(1)
  
  #sample distribution bin
  i_bin <- sum(Pcum < Pinv) + 1 
  
  #sample x based on inverse prob
  a_coef <- H_norm[i_bin] - X[i_bin]*H_prime[i_bin] #intersect of each line seg. (log space)
  b_coef <- H_prime[i_bin] #slope of each line seg. (log space)
  DP <- Pinv - max(Pcum[i_bin-1],0) 
  
  #sample x_star
  x_star <- 1/b_coef*log(DP*b_coef*exp(-a_coef)+exp(b_coef*Z[i_bin]))
  
  return(x_star)
}


#Calculate probability of each segment
CalcProbBin <- function(X,Z,H,H_prime){
  #CalcProbBin calculates the probabilities of each segment of the exponential distribution.
  #Input: 
  # X: numeric array with all the values where h and h' have been evaluated
  # Z: numeric array with the ordinates at the intersections of the piecewise exponential function u
  # H: numeric array with the evaluations of function h at the locations X
  # H_prime:  numeric array with the evaluations of the first derivative of h at locations X
  #Output:
  # Pcum: numeric array with the cumulative probability of each of each segment 
  # log_s: numeric array with the normalized elevations of log(s) at locations X
  
  #number of bins
  n_bin <- length(Z) -1
  
  #calculate area under each bin (proportional to probabilty)
  a_coef <- H - X*H_prime #intersect of each line seg. (log space)
  b_coef <- H_prime #slope of each line seg. (log space)
  P <- exp(a_coef)/b_coef*(exp(b_coef*Z[-1]) - exp(b_coef*Z[-(n_bin+1)]))
  #normalized probabilites
  Ptot <- sum(P)
  P <- P/Ptot
  #cumulative probabilites
  Pcum <- cumsum(P) 
  #normalized height of distributions
  log_s <- H -log(Ptot)
  
  return(list(Pcum = Pcum, log_s = log_s, Ptot = Ptot))
}

#Accept or reject sample 
#Note: added H, H_prime, Z as parameters to function
update_accept <- function(x_star, g, X, H, H_prime, Z,length_accept,x_accept,H_norm,Pcum) {
  
  l <- function(x,i) {
    if (x < X[length(X)] & x >= X[1]){ #x in [X[i],X_[i+1] 
      leval <- ((X[i+1]-x)*H[i] + (x-X[i])*H[i+1])/(X[i+1]-X[i])
    } else if (x == X[length(X)]) { 
      i <- i -1
      leval <- ((X[i+1]-x)*H[i] + (x-X[i])*H[i+1])/(X[i+1]-X[i])
    } else { #x < X[1] or x > X[k]
      leval <- -Inf
    }
    return(leval)
  }
  
  u <- function(x,i) {H[i] + (x-X[i])*H_prime[i]}
  
  z <- function(i) {
    (H[i+1]-H[i]-X[i+1]*H_prime[i+1]+X[i]*H_prime[i])/(H_prime[i]-H_prime[i+1])
  }
  w <- runif(1)
  
  i_star_x <- max(which.max(X[X<=x_star]),0)
  i_star_z <- length(Z[x_star>=Z])
  #print(paste0("x_star is: ",x_star))
  #print(paste0("X :",list(X)))
  #print(paste0("Z is: ",list(Z)))
  #print(paste0("i_star_x is: ", i_star_x))
  #print(paste0("i_star_z is: ", i_star_z))
  #print(paste0("L is: ",l(x_star,i_star_x)))
  #print(paste0("U is: ",u(x_star,i_star_z)))
  
  if(w <= exp(l(x_star,i_star_x) - u(x_star,i_star_z))){
    x_accept[length_accept] <- x_star
    length_accept <- length_accept+1
    #print(paste0("No update. New x_accept"))
    
  } else {
    #Update step
    
    X <- append(X, x_star,i_star_x)
    H <- append(H,update_H(x_star,g),i_star_x)
    H_prime <- append(H_prime,update_H_prime(x_star,g),i_star_x)
    
    if (x_star != X[1] & x_star != X[length(X)]){ #x_star outside other X
      Z[i_star_x+1] <- z(i_star_x) #Overwrite a z
      Z <- append(Z,z(i_star_x+1),i_star_x+1) #Append a new z
    } else {
      Z <- append(Z,z(i_star_z),i_star_z)
    }
    
    #calculate new cumulative probabilities of the piece-wise exponential bins
    out <- CalcProbBin(X,Z,H,H_prime)
    Pcum <- out[[1]]
    H_norm <- out[[2]]
    Ptot <- out[[3]]
    
    
    #print(paste0("Abscissae is updated"))
    #print(paste0("New H is updated"))
    #print(paste0("New H_prime is updated"))
    
    #accept x_star second test
    if (w <= exp(log(g(x_star)) - u(x_star,i_star_z))) {
      x_accept[length_accept] <- x_star
      length_accept <- length_accept+1
      #print(paste0("Update step is done. New x_accept"))
    }
  }

  
  return(list(X = X,H = H,H_prime = H_prime,Z = Z,x_accept = x_accept,k = length_accept,H_norm = H_norm,Pcum = Pcum))
  
}



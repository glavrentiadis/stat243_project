
ars <- function(g,n_samp,x_bound,x_start = c(),sybolic_deriv = FALSE){
  
  #load libraries
  library(assertthat)
  library(testthat)
  library(truncnorm)
#  library(Ryacas)
  library(Deriv)

  #Function handle for h
  eval_H <- EvalH(g)
  #Choose if numeric or symbolic derivative 
  if (sybolic_deriv == FALSE){
    eval_H_prime <- NumericEvalHPrime(g)
  } else {
    eval_H_prime <- SymbolicEvalHPrime(g)
  }
  
  #initialization step
  #chech starting points if given by the user
  if (length(x_start) != 0){
    out <- CheckXStart(x_start = x_start, x_bound = x_bound, eval_h = eval_H, eval_deriv_h = eval_H_prime)
  } else {
    out <- SampleXStart(x_bound = x_bound, eval_h = eval_H, eval_deriv_h = eval_H_prime)
  }
  X <- out$X
  Z <- out$Z
  H <- out$H
  H_prime <- out$H_prime
  H_norm <- out$H_norm
  P_cum <- out$Pcum


  #Sampling and updating when needed
  x_accept <- numeric(n_samp)
  k <- 0
  while (k <= n_samp){
    #sample x_star
    x_star <- SamplePieceExp(X,Z,H_norm,H_prime,P_cum)
    #print(paste('Master: x_star:',x_star))
    #check if x_star will be kept and update X,Z,H,H_prime if needed
    out <- UpdateAccept(x_star = x_star, X = X ,Z = Z, H = H, H_prime = H_prime, H_norm = H_norm, P_cum = P_cum,
                        x_accept = x_accept, length_accept = k, eval_h = eval_H, eval_h_prime = eval_H_prime)
    X <- out$X
    Z <- out$Z
    H <- out$H
    H_prime <- out$H_prime
    H_norm <- out$H_norm
    P_cum <- out$P_cum
    x_accept <- out$x_accept
    k <- out$k
  }
  
  #return random samples
  return(x_accept)

}


# Auxilary functions
# ----------------------------------------------------------------------------

#function to compute log of g (function h)
EvalH <- function(g){
  
  evalH <- function(x){
    return(log(g(x)))
  }
  
  return(evalH)
}

#fuction to compute the numerically the derivative of h
NumericEvalHPrime <- function(g){
  
  evalH_prime <- function(x){
    #numerical step
    dx <- sqrt(.Machine$double.eps)*abs(max(x,1))
    #compute numerical derivative
    dh_dx <- (log(g(x+dx))-log(g(x-dx)))/(2*dx)

#    assert_that(!is.na(dh_dx))
    return(dh_dx)    
  }
  
  return(evalH_prime)
}

#fuction to compute the symbolicaly the derivative of h
SymbolicEvalHPrime <- function(g){
  
  #compute the derivative of g symbolically
  g_prime_expr <- Deriv(~g(x),"x")
  
  evalH_prime <- function(x){
        return(as.numeric(1/g(x)*eval(g_prime_expr)))
  }
  
  return(evalH_prime)
}

#Check staring points
CheckXStart <- function(x_start,x_bound,eval_h,eval_deriv_h){
  #CheckXStart checks that the staring x points satisfy the requirements
  #and creates the inital arrays (X ,H ,H_prime , Z ,Pcum ,H_norm, Ptot)
  
  #functions to test check the inital points
  #first point check for positve first derivative
  first_point_check <- function(H_d){
    if(H_d[1]>0){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  on_failure(first_point_check) <- function(call, env) {
    print("Invalid starting points, for a left unbounded distribution first derivative of first starting point must be positive")
  }
  #last point check for negative first derivative
  last_point_check <- function(H_d){
    k_len <- length(H_d)
    if(H_d[k_len]<0){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  on_failure(last_point_check) <- function(call, env) {
    print("Invalid starting points, for a right unbounded distribution first derivative of last starting point must be negative")
  }
  
  #assert x_bound has to elements
  assert_that(length(x_bound) == 2)
  #assert that at least two starting points are specified
  assert_that(length(x_start) >= 2)
  #assert that x_star is in ascending order
  assert_that(all(diff(x_start)>0))

  #number of starting points
  k_start <- length(x_start)
  
  #compute H and H derivative
  X <- x_start
  H <- eval_h(X)
  H_d <- eval_deriv_h(X)

  #if lower derivative is -Inf, the derivative of the last point must be positive
  if (is.infinite(x_bound[1])){
    assert_that(first_point_check(H_d))
  }
  #if upper derivative is +Inf, the derivative of the last point must be negative
  if (is.infinite(x_bound[2])){
    assert_that(last_point_check(H_d))    
  }
    
  #evaluate the intersection points
  Z <- (H[-1] - H[-k_start] - X[-1]*H_d[-1] + X[-k_start]*H_d[-k_start])/(H_d[-k_start] - H_d[-1])
  #correct points with the same H'
  i_con_Hd <- abs(H_d[-1] - H_d[-k_start]) < .Machine$double.eps*abs(H_d[-k_start])
  #for these points assign Z_i as (X_i+X_{i+1})/2
  Z_med <- 0.5*(X[-1] + X[-k_start])
  Z[i_con_Hd] <- Z_med[i_con_Hd]
  
  #include the bounds in the array Z
  Z <- c(x_bound[1],Z,x_bound[2])
  
  #compute the probabilites of each piece-wise exponential bin
  out <- CalcProbBin(X,Z,H,H_d)
  Pcum <- out[[1]]
  H_norm <- out[[2]]
  Ptot <- out[[3]]
  
  
  return(list(X = X,H = H,H_prime = H_d,Z = Z,Pcum = Pcum,H_norm = H_norm,Ptot = Ptot))
}


#Sample initial points if they are not given 
SampleXStart <- function(x_bound,eval_h,eval_deriv_h){
  #SampleXStart samples two starting points if staring points are not
  #given by the user


  print("Sampling starting points, if it takes long to complete this step please specify x_start to avoid this step")
  
  #assert x_bound has to elements
  assert_that(length(x_bound) == 2)

  #parse lower and upper distribution limits for sampling
  lowBd_samp <- x_bound[1] 
  upBd_samp <- x_bound[2]
  
  #initialize x_start
  x_start <- numeric(2)
#  browser()
  #sample first point, iterate until H_d > 0, if lowBd = -Inf
  flag_bound_req <- FALSE
  mean_samp <- 0
  std_samp <- 1
  
  #distribution boundaries
  lowBd_samp <- x_bound[1] 
  upBd_samp <- x_bound[2]
  
  while (flag_bound_req == FALSE){
    samp <- rtruncnorm(1, a=lowBd_samp, b=upBd_samp, mean = mean_samp, sd = std_samp)
    hPrime <- eval_deriv_h(samp)

    if (!(!is.na(hPrime) & is.finite(hPrime))) {
      lowBd_samp <- samp
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp + offset_val
    } else if (hPrime > 0){
      flag_bound_req <- TRUE #boundary requirement is satisfied
      x_start[1] <- samp 
    } else {
      #modify sampling distribution parameters move to the left and widen
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp - offset_val
      upBd_samp <- samp
    }
    
  }
    
  #sample second point, iterate until H_d < 0, if lowBd = +Inf
  flag_bound_req <- FALSE
  mean_samp <- 0
  std_samp <- 1
  
  #distribution boundaries
  lowBd_samp <- x_start[1]
  upBd_samp <- x_bound[2]
  
  while (flag_bound_req == FALSE){
    samp <- rtruncnorm(1, a=lowBd_samp, b=upBd_samp, mean = mean_samp, sd = std_samp)
    hPrime <- eval_deriv_h(samp)


    if (!(!is.na(hPrime) & is.finite(hPrime))) {
      upBd_samp <- samp
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp - offset_val
    } else if (hPrime < 0){
      flag_bound_req <- TRUE #boundary requirement is satisfied
      x_start[2] <- samp
    } else {
      #modify sampling distribution parameters move to the right and widen
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp + offset_val
      lowBd_samp <- samp
    }
  }
  
  #sort starting points in acceding order
  x_start <- sort(x_start)
  
  #compute starting values for (X ,H ,H_prime , Z ,Pcum ,H_norm, Ptot)
  out <- CheckXStart(x_start,x_bound,eval_h,eval_deriv_h)
  
  print("Sampling of starting points is complete")
  
  return(out)
  
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
  #a_coef <- H_norm[i_bin] - X[i_bin]*H_prime[i_bin] #intersect of each line seg. (log space)
  #b_coef <- H_prime[i_bin] #slope of each line seg. (log space)
  H_norm_st <- H_norm[i_bin] #intersect of each line seg. (log space)
  H_prime_st <- H_prime[i_bin] #slope of each line seg. (log space)
  
  DP <- Pinv - max(Pcum[i_bin-1],0) 


  #sample x_star
  #x_star <- 1/b_coef*log(DP*b_coef*exp(-a_coef)+exp(b_coef*Z[i_bin]))
  if (abs(H_prime_st) >sqrt(.Machine$double.eps)){
    x_star <- 1/H_prime_st*(log(DP*H_prime_st + exp(H_norm_st + H_prime_st*(Z[i_bin]-X[i_bin]))) + H_prime_st*X[i_bin] - H_norm_st)
  } else {
    DP_scl <- DP/Pcum[i_bin]
    x_star <- Z[i_bin] + diff(Z[i_bin:(i_bin+1)])*DP_scl
  }
    

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
  #P <- exp(a_coef)/b_coef*(exp(b_coef*Z[-1]) - exp(b_coef*Z[-(n_bin+1)])) 
  P <-  (exp(H + H_prime*(Z[-1] - X)) - exp(H + H_prime*(Z[-(n_bin+1)] - X)) )/H_prime 
  
  #probabilities for bins with zero tangent (the above equation is not valid for this case)
  bin_zero_Hprime <- abs(H_prime) <sqrt(.Machine$double.eps)
  P[bin_zero_Hprime] <- (Z[-1][bin_zero_Hprime] - Z[-length(Z)][bin_zero_Hprime])*exp(a_coef[bin_zero_Hprime])
  
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
UpdateAccept <- function(x_star, X ,Z , H, H_prime, H_norm, P_cum, x_accept, length_accept ,eval_h ,eval_h_prime){
  

  l <- function(x,i) {
    # print(paste0("test1: ",x < X[length(X)]))
    # print(paste0("test2: ",x >= X[1]))
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
    if (abs(H_prime[i]-H_prime[i+1]) < .Machine$double.eps*max(abs(H_prime[i]),1)){
      return(0.5*(X[i]+X[i+1]))
#      return((H[i+1]-H[i]-X[i+1]*H_prime[i+1]+X[i]*H_prime[i])/(H_prime[i]-H_prime[i+1]))
    } else {
      return((H[i+1]-H[i]-X[i+1]*H_prime[i+1]+X[i]*H_prime[i])/(H_prime[i]-H_prime[i+1]))
    }
  }
  
  #function to test if function is concave 
  concave_check <- function(u_at_xstar,l_at_xstar,tol){
    # assert_that(l_at_xstar - u_at_xstar  > tol)
    tol <- sqrt(.Machine$double.eps)*max(abs(u_at_xstar),1) #tolerance for concavity check
    if(l_at_xstar - u_at_xstar  < tol){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  on_failure(concave_check) <- function(call, env) {
    print(paste0("h(x) doesn't satisfy the concavity requirements at x*: ",x_star,", u(x*): ",u_at_xstar,", l(x*): ",l_at_xstar))
  }


  w <- runif(1)
  
  i_star_x <- max(which.max(X[X<=x_star]),0)
  i_star_z <- length(Z[x_star>=Z])
  # print(paste0("x_star is: ",x_star))
  #  print(paste0("X :",list(X)))
  # print(paste0("Z is: ",list(Z)))
  # print(paste0("i_star_x is: ", i_star_x))
  # print(paste0("i_star_z is: ", i_star_z))
  # print(paste0("L is: ",l(x_star,i_star_x)))
  # print(paste0("U is: ",u(x_star,i_star_z)))
  
  #concavity check
  u_at_xstar <- u(x_star,i_star_z)
  l_at_xstar <- l(x_star,i_star_x)
  assert_that(concave_check(u_at_xstar,l_at_xstar))

  
  if(w <= exp(l_at_xstar - u_at_xstar)){
    x_accept[length_accept] <- x_star
    length_accept <- length_accept+1
    # print(paste0("No update. New x_accept"))
    
  } else {
    #Update step
    X <- append(X, x_star,i_star_x)
    H <- append(H,eval_h(x_star),i_star_x)
    H_prime <- append(H_prime,eval_h_prime(x_star),i_star_x)
    
    if (x_star > X[1] & x_star < X[length(X)]){ #x_star inside other X
      Z[i_star_x+1] <- z(i_star_x) #Overwrite a z
      Z <- append(Z,z(i_star_x+1),i_star_x+1) #Append a new z
    } else {
      # if (sum(Z<=z(i_star_z))!=i_star_z){
      #   browser()
      # }
      Z <- append(Z,z(i_star_z),i_star_z)
    }
    # if (any(diff(Z) <0)) {
    #   browser()
    # }
    # k_temp <- length(X)
    # Z_temp <- (H[-1] - H[-k_temp] - X[-1]*H_prime[-1] + X[-k_temp]*H_prime[-k_temp])/(H_prime[-k_temp] - H_prime[-1])
    # if (any(abs(Z_temp - Z[2:(k_temp)]) >1e-5)){
    #   browser()
    #   
    # }
    
    
    #calculate new cumulative probabilities of the piece-wise exponential bins
    out <- CalcProbBin(X,Z,H,H_prime)
    P_cum <- out[[1]]
    H_norm <- out[[2]]
    Ptot <- out[[3]]
    
    
   # print(paste0("Abscissae is updated"))
   # print(paste0("New H is updated"))
   # print(paste0("New H_prime is updated"))
   # print(paste0("H is: ",cat(H)))
   # print(paste0("H_prime is: ",cat(H_prime)))
    
    #accept x_star second test
    if (w <= exp(eval_h(x_star) - u_at_xstar)) {
      x_accept[length_accept] <- x_star
      length_accept <- length_accept+1
      # print(paste0("Update step is done. New x_accept"))
    } else {
      # print(paste0("Update step is done. No new x_accept"))
    }
  }

  return(list(X = X, Z = Z, H = H, H_prime = H_prime, H_norm = H_norm, P_cum = P_cum, x_accept = x_accept, k = length_accept))
  
}



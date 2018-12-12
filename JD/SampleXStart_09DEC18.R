#' Initialize .
#' 

#' add(10, 1)

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
  
  #sample first point, iterate until H_d > 0, if lowBd = -Inf
  flag_bound_req <- FALSE
  mean_samp <- 0
  std_samp <- 1
  while (flag_bound_req == FALSE){
    samp <- rtruncnorm(1, a=lowBd_samp, b=upBd_samp, mean = mean_samp, sd = std_samp)
    hPrime <- eval_deriv_h(samp)
    
    print(paste0("samp: ",samp))
    print(paste0("hPrime: ",hPrime))
    if (eval_h(samp)==Inf || is.nan(eval_deriv_h(samp))){
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp + offset_val
      std_samp <- std_samp*offset_val
      lowBd_samp <- samp
    } else if (!is.infinite(x_bound[1]) || hPrime > 0){
      flag_bound_req <- TRUE #boundary requirement is satisfied
      x_start[1] <- samp
      print("Got first point")
    } else {
      #modify sampling distribution parameters move to the left and widen
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp - offset_val
      std_samp <- std_samp*offset_val
      upBd_samp <- samp
    }
  }
  
  #sample second point, iterate until H_d < 0, if lowBd = +Inf
  flag_bound_req <- FALSE
  mean_samp <- 0
  std_samp <- 1
  while (flag_bound_req == FALSE){
    samp <- rtruncnorm(1, a=lowBd_samp, b=x_bound[2], mean = mean_samp, sd = std_samp)
    hPrime <- eval_deriv_h(samp)
    print(paste0("samp: ",samp))
    print(paste0("hPrime: ",hPrime))
    print(typeof(hPrime))
    
    
    if (eval_h(samp)==Inf || is.nan(eval_deriv_h(samp))){
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp - offset_val
      std_samp <- std_samp*offset_val
      upBd_samp <- samp
    } else if (!is.infinite(x_bound[2]) || hPrime < 0){
      flag_bound_req <- TRUE #boundary requirement is satisfied
      x_start[2] <- samp
    } else {
      #modify sampling distribution parameters move to the right and widen
      print("Before mod:")
      print(paste0(c("mean_samp: ",mean_samp)))
      print(paste0(c("std_samp: ",std_samp)))
      offset_val <- max(abs(mean_samp - samp),2)
      mean_samp <- mean_samp + offset_val
      std_samp <- min(10,std_samp*offset_val)
      lowBd_samp <- samp
      print("After mod:")
      print(paste0(c("mean_samp: ",mean_samp)))
      print(paste0(c("std_samp: ",std_samp)))
    }
  }
  
  #sort starting points in acceding order
  x_start <- sort(x_start)
  
  #compute starting values for (X ,H ,H_prime , Z ,Pcum ,H_norm, Ptot)
  out <- CheckXStart(x_start,x_bound,eval_h,eval_deriv_h)
  
  print("Sampling of starting points is complete")
  
  return(out)
  
}


library(stringr)
str_detect
str_count


g <- function(x){
  mu <- 5
  std <- 4
  
  f_pdf <- 1/sqrt(2*pi*std^2)*exp(-(x-mu)^2/(2*std^2))
}

n_samp <- 100000
x_b <- c(-Inf,Inf)

#system.time(
temp <- ars(g, n_samp, x_bound = x_b,  sybolic_deriv = FALSE)
#)



#g <- function(x){
  
  #eval <- exp(x^2)
  
#}





mean(temp)
sd(temp)

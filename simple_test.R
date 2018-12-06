

g <- function(x){
  
  mu <- 0
  std <- 1
  
  f_pdf <- 1/sqrt(2*pi*std^2)*exp(-(x-mu)^2/(2*std^2))
  
}

log(g(2))

n_samp <- 10000
x_b <- c(-Inf,Inf)

system.time(
temp <- ars(g,n_samp,x.Bound = x_b ,k_start = 2)
)


mean(temp)
sd(temp)

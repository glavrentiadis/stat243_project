getwd()

library(devtools)

setwd("./ars")
devtools::document()



setwd("..")
install("ars")

library(ars)

?ars
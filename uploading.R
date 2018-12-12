getwd()

setwd("./ars")
document()


library(devtools)
devtools::document()

setwd("..")
install("ars")

library(ars)

?ars


#install.packages("devtools")
library(devtools)
setwd("/Users/jbunker7/Documents/GitHub/stat243_project")
usethis::create_package("Package")



setwd("/Users/jbunker7/Documents/GitHub/stat243_project/Package")
#?Package
devtools::document()
setwd("..")
getwd()
?help
?.
#install("Package")

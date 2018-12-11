# Mehrnoosh Oghbaie
# 12/10/2018
# Repository for all the functions
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  if (length(bioconductor.packages) > 0) {
    source("http://bioconductor.org/biocLite.R")
  }
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      biocLite(p) 
      library(p, character.only=T)
    }
  }
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
      library(p, character.only=T)  	
    }	
  }
}

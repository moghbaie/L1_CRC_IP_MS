## Mehrnoosh Oghbaie
## 12/11/2018
## Integrate normalized Values

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2")
bioconductor.packages <- c("biomaRt")
source("All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Normalization consists of two stages:
##  1. Normalize intensities by gfp intensity
##  2. Take advantage of normalized intensity on all cases

####################################################################################################################
## Save the output in "../Result/5. Integrated_normalized_data/"
####################################################################################################################

df <- data.frame(matrix(NA, ncol=3))
write.csv(df,file="../Result/5. Integrated_normalized_data/test.csv")

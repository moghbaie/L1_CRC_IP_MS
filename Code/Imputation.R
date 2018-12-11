## Mehrnoosh Oghbaie
## 12/11/2018
## Imputate missing values

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

## Imputation consists of three stages:
##  1. Remove proteins not identofoed in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity



####################################################################################################################
## Save the output in "../Result/3. After_imputation/"
####################################################################################################################

df <- data.frame(matrix(NA, ncol=3))
write.csv(df,file="../Result/3. After_imputation/test.csv")

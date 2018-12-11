## Mehrnoosh Oghbaie
## 12/11/2018
## Preparing data from MaxQuant or ...

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

## Data preparation consists of four stages:
##  1. Calculate iBAQ intensity
##  2. Remove contaminants and reverse proteins
##  3. Log transformation
##  4. Separate different experiment



####################################################################################################################
## Save the output in "../Result/2. After_preparation/"
####################################################################################################################

df <- data.frame(matrix(NA, ncol=3))
write.csv(df,file="../Result/2. After_preparation/test.csv")

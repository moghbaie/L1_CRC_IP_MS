## Mehrnoosh Oghbaie
## 12/11/2018
## Quality control test

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","magrittr","sqldf","PTXQC","corrplot","ggplot2")
bioconductor.packages <- c("biomaRt")
source("All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## QC steps



####################################################################################################################
## QC plot "../Image/QC_plot/"
####################################################################################################################



####################################################################################################################
## Save the output in "../Result/1. QC/"
####################################################################################################################



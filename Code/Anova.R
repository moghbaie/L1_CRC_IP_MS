## Mehrnoosh Oghbaie
## 12/11/2018
## ANOVA test between cases and controls

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

## Variance analysis consists of four stages:
##  1. Filter proteins with less than two peptides
##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)
##  5. Draw Volcano plot



####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################



####################################################################################################################
## Save the output in "../Result/4. ANOVA_result/"
####################################################################################################################

write.csv(df,file="../Result/4. ANOVA_result/Significant_table.RData")

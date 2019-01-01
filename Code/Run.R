## Mehrnoosh Oghbaie
## 12/31/2018
## Define class and run the project

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

CRAN.packages <- c("readr","readxl","data.table","reshape2","dplyr","magrittr",
                   "igraph","sqldf","stringr","corrplot","ggplot2","R6","ggridges",
                   "gridExtra","ggrepel","rgl")
bioconductor.packages <- c("biomaRt","limma")

source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("Data_preparation.R")
source("Imputation.R")
source("Anova.R")
source("Normalize.R")

run_order <- rbind(cbind(rep("Ovary",2),c("Igg","Tumor")),
                   cbind(rep("Liver",3),c("Igg","Normal","Tumor")),
                   cbind(rep("Colon",3),c("Igg","Normal","Tumor")))

input.info.dir <- "../Input_data/Input.info"

CRC <- Template$new()
CRC$removeContaminant()
CRC$logTransformation()
CRC$separatedList()
CRC$removeAllZeros()
CRC$imputeAll(run_order)
CRC$drawComparisonAfterImputation()
CRC$drawAverageImputated()
CRC$anovaAnalysis(run_order)
CRC$drawVenDiagramSignificant()
CRC$drawVolcanoPlot() #or run the following commands at the end (it might crash)

#draw_volcanoplot(data=CRC[["experimentImputed"]],condition="Ovary_Tumor_Igg")
#draw_volcanoplot(data=CRC[["experimentImputed"]], condition="Liver_Tumor_Igg")
#draw_volcanoplot(data=CRC[["experimentImputed"]], condition="Liver_Tumor_Normal")
#draw_volcanoplot(data=CRC[["experimentImputed"]], condition="Colon_Tumor_Igg")
#draw_volcanoplot(data=CRC[["experimentImputed"]], condition="Colon_Tumor_Normal")

CRC$calculateAverageNormalized(run_order)
CRC$drawCommonVenndiagram()
CRC$drawHeatmap()

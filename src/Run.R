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
                   "gridExtra","ggrepel","rgl","venn")
bioconductor.packages <- c("biomaRt","limma","qvalue","msa","ape","seqinr","ggseqlogo")


############################################################################################
## Check the values in input.info file

source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("Data_preparation.R")
source("Imputation.R")
source("Anova.R")
source("Normalize.R")
source("Mutation.R")

eLife_list <- read_excel("../Input_data/eLife_list.xlsx", 
                         col_names = FALSE)
colnames(eLife_list) <- c("geneName","uniprotID","Description","Condition")

run_order2 <- rbind(cbind(rep("Ovary",4),c("144T_IgG_2","144T_IgG_10","144T_ORF1_1","144T_ORF1_9")),
                   cbind(rep("Liver",3),c("159T_IgG_4","159N_ORF1_5","159T_ORF1_3")),
                   cbind(rep("Colon",3),c("163T_IgG_7","163N_ORF1_8","163T_ORF1_6")))
run_order2 <- cbind(run_order2,c(1,2,1,2,1,1,1,1,1,1))

## Data_preparation.R
CRC <- Template$new()
CRC$removeContaminant()
CRC$logTransformation()
CRC$separatedList()

## Imputation.R
CRC$removeAllZeros()
CRC$imputeAll(run_order2)
CRC$imputePerseusMethod(run_order2)
CRC$drawComparisonAfterImputation()
CRC$drawAverageImputated(run_order2)

## Anova.R
CRC$anovaAnalysisPerseus(run_order2)
CRC$anovaAnalysis(run_order2)
CRC$anovaAnalysisPreImpute(run_order2)
CRC$drawVenDiagramSignificant()
#CRC$drawVolcanoPlot() #or run the following commands at the end (it might crash)

draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Ovary_Tumor_144T_IgG_10", eLife_list)
draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Liver_Tumor_159T_IgG_4", eLife_list)
draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Liver_Tumor_159N_ORF1_5", eLife_list)
draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Colon_Tumor_163T_IgG_7", eLife_list)
draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Ovary_Tumor_144T_IgG_2", eLife_list)
draw_volcanoplot(data=CRC[["experimentImputateSignificant"]], condition="Colon_Tumor_163N_ORF1_8", eLife_list)


CRC$calculateAverageNormalized(run_order2)
CRC$drawCommonVenndiagram()
CRC$drawHeatmap()
CRC$drawHeatmapUnnormalized()
CRC$drawMutation()

#save(CRC,file = "../backup.RData")
save(CRC,file = "../backup_v4.RData")

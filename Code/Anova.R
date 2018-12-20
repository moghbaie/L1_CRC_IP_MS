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
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Variance analysis consists of four stages:
##  1. Filter proteins with less than two peptides
##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)
##  5. Draw Volcano plot

load(file="../Result/3. After_imputation/after_imputation.RData")

##  1. Filter proteins with less than two peptides
# I'll skip this one for now


##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)

run_order <- rbind(cbind(rep("Ovary",2),c("Igg","Tumor")),
                   cbind(rep("Liver",3),c("Igg","Normal","Tumor")),
                   cbind(rep("Colon",3),c("Igg","Normal","Tumor")))

for( i in unique(run_order[,1])){
  y <- list_experiment[[i]]
  z <- run_order[run_order[,1]==i,2]
  for(j in z[z!="Tumor"]){
    
      dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
      dt[["ID"]] <- y[["id"]]
      dt[["uniprotID"]] <- y[["uniprotID"]]
      dt[["p.value"]] <- NA
      dt[["logfold"]] <- NA
      dt[["Significant"]] <- NA
      for(k in 1:dim(dt)[1]){
        dt[[k,"p.value"]] <- t.test(unname(unlist(y[k,grepl("Tumor",colnames(y))])),unname(unlist(y[k,grepl(j,colnames(y))])))$p.value
        dt[[k,"logfold"]] <- mean(unname(unlist(y[k,grepl("Tumor",colnames(y))])))-mean(unname(unlist(y[k,grepl(j,colnames(y))])))
      }
      dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
      dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&dt[["logfold"]]>1,"Yes","No")
      list_experiment[[paste0(i,"_","Tumor","_",j)]] <-dt
  }
}


Significant_proteins <- c()
for(i in names(list_experiment)[grep("_",names(list_experiment))]){
  Significant_proteins <- unique(c(Significant_proteins,list_experiment[[i]][list_experiment[[i]][["Significant"]]=="Yes","uniprotID"]))
}

####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################



####################################################################################################################
## Save the output in "../Result/4. ANOVA_result/"
####################################################################################################################

save(Significant_proteins, list_experiment, file="../Result/4. ANOVA_result/Significant.RData")

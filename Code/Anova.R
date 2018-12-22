## Mehrnoosh Oghbaie
## 12/11/2018
## ANOVA test between cases and controls

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2","readxl","limma","ggrepel")
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


####################################################################################################################
## Venn diagram
Significant_proteins <- data.frame(matrix(NA, nrow=length(KM050217_merged$uniprotID), ncol= 2+length( names(list_experiment)[grep("_",names(list_experiment))])))
colnames(Significant_proteins) <- c("id","uniprotID", names(list_experiment)[grep("_",names(list_experiment))])
Significant_proteins$uniprotID <- KM050217_merged$uniprotID
Significant_proteins$id <- KM050217_merged$id


for(i in names(list_experiment)[grep("_",names(list_experiment))]){
  Significant_proteins[[i]] <- list_experiment[[i]][,"Significant"][match(Significant_proteins$uniprotID,list_experiment[[i]][,"uniprotID"])]
  Significant_proteins[[i]][is.na(Significant_proteins[[i]])] <- "NO"
  Significant_proteins[[i]] <- ifelse(Significant_proteins[[i]]=="Yes",1,0)
}

a <- vennCounts(Significant_proteins[,-(1:2)])
vennDiagram(a,circle.col = c("darkmagenta", "darkblue", "pink1", "skyblue", "orange"), main="Comparison significant proteins")

Significant_list <- Significant_proteins%>%
  filter(.[3]|.[4]|.[5]|.[6]|.[7]==1)%>%
  dplyr::select(id,uniprotID)


####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################


draw_volcanoplot <- function(data=list_experiment, condition="Liver_Tumor_Igg"){
  ds <- data[[condition]]
  fold_cutoff = 1
  pvalue_cutoff = 0.05
  png(paste0("../Image/Volcano_plot/Volcano_plot_",condition,".png"),width = 800, height = 600)
  p <- ggplot(ds, aes(logfold, -log10(p.adj), label = ifelse(Significant=="Yes"|ID=="LORF1",ID,""))) +
    geom_point(color = ifelse(ds$Significant == "Yes"|ds$ID=="LORF1", "red", "grey50")) +
    geom_text_repel()+
    geom_vline(xintercept = fold_cutoff, col = "blue")+
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
    ggtitle(condition)
  print(p)
  dev.off()
}


names(list_experiment)[grep("_",names(list_experiment))]
draw_volcanoplot(condition="Ovary_Tumor_Igg")
draw_volcanoplot(condition="Liver_Tumor_Igg")
draw_volcanoplot(condition="Liver_Tumor_Normal")
draw_volcanoplot(condition="Colon_Tumor_Igg")
draw_volcanoplot(condition="Colon_Tumor_Normal")


####################################################################################################################
## Save the output in "../Result/4. ANOVA_result/"
####################################################################################################################

save(KM050217_merged,Significant_proteins,Significant_list, list_experiment, file="../Result/4. ANOVA_result/Significant.RData")


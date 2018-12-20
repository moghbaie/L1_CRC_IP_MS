## Mehrnoosh Oghbaie
## 12/11/2018
## Preparing data from MaxQuant or ...

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","readxl","data.table","reshape2","dplyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2")
bioconductor.packages <- c("biomaRt")
source("Function/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Data preparation consists of four stages:
##  1. Calculate iBAQ intensity
##  2. Remove contaminants and reverse proteins
##  3. Log transformation
##  4. Separate different experiment


KM050217_merged <- read_excel("../Input_data/Data_from_Kelly/KM050217_merged.xlsx")
dim(KM050217_merged)
cols <- colnames(KM050217_merged)

## 1. iBAQ intensity can be extracted
KM050217_merged[,cols[grepl("iBAQ",cols)][-1]]

## 2. Contaminants can be removed using contaminats txt file
contaminants <- read.table("F:/Line1_CRC/Input_data/Data_from_Kelly/contaminants.fasta", sep=";", quote="\"")
contaminant_list <- unlist(lapply(as.character(contaminants[,1])[grepl(">", as.character(contaminants[,1]))], function(x) strsplit(x," |>")[[1]][2]))

# extract ID from fasta folder
KM050217_merged$uniprotID <- unlist(lapply(KM050217_merged$`Fasta headers`, function(x) strsplit(as.character(x),"\\|")[[1]][2]))

# Remove contaminants from the dataframe
df <- KM050217_merged[!KM050217_merged$uniprotID %in% contaminant_list,c("id","uniprotID",cols[grepl("iBAQ",cols)][-1])]
dim(df)

##  3. log transformation
tp <- lapply(df, class)
df_log <- df


for(i in colnames(df)[tp =="numeric"]){
  df_log[i] <- log(unlist(unname(df_log[i])))
  df_log[is.infinite(unname(unlist(df_log[i]))),i] <- NA
}

##  4. Separate different experiment
colnames(df_log) <- c("id","uniprotID",unlist(lapply(colnames(df_log)[-(1:2)], function(x) strsplit(x,"_")[[1]][2])))

list_experiment <- list()
list_experiment[["Ovary"]] <- df_log[, grepl("1|2|uniprotID|id",colnames(df_log))]
colnames(list_experiment[["Ovary"]]) <- paste0(colnames(list_experiment[["Ovary"]]),c(rep("",2),rep("_Tumor",3),rep("_Igg",3)))


list_experiment[["Liver"]] <- df_log[, grepl("3|4|5|uniprotID|id",colnames(df_log))]
colnames(list_experiment[["Liver"]]) <- paste0(colnames(list_experiment[["Liver"]]),c(rep("",2),rep("_Normal",3),rep("_Tumor",3),rep("_Igg",3)))

list_experiment[["Colon"]] <- df_log[, grepl("6|7|8|uniprotID|id",colnames(df_log))]
colnames(list_experiment[["Colon"]]) <- paste0(colnames(list_experiment[["Colon"]]),c(rep("",2),rep("_Normal",3),rep("_Tumor",3),rep("_Igg",3)))

####################################################################################################################
## Save the output in "../Result/2. After_preparation/"
####################################################################################################################

save(list_experiment,file="../Result/2. After_preparation/log_transformed.RData")


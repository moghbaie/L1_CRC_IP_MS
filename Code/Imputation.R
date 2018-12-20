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
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Imputation consists of three stages:
##  1. Remove proteins not identified in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity

load(file="../Result/2. After_preparation/log_transformed.RData")

#####################################################################################################
##  1. Remove proteins not identified in both cases and controls
# Ovary
list_experiment$Ovary <- list_experiment$Ovary[rowSums(list_experiment$Ovary[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Ovary)

# Liver
list_experiment$Liver <- list_experiment$Liver[rowSums(list_experiment$Liver[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Liver)

# Colon
list_experiment$Colon <- list_experiment$Colon[rowSums(list_experiment$Colon[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Colon)

#####################################################################################################
##  2. Impute values for proteins with zero intensity (cases and controls)

### function : finding column with minimum zeros
count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_mutate <- function(na_zeros , mu,sd){
  matrix(runif(3*sum(na_zeros), mu-3*sd, mu-2*sd),ncol=3)
}

### function: mutate all NA values for specific condition( control or case)
mutate_na_zeros <- function(data=list_experiment$Liver,condition="Normal"){
  na_zeros <- rowSums(data[,grepl(condition,colnames(data))],na.rm=T)==0
  min_zero_col <- min_col(apply(data[,grepl(condition,colnames(data))],2,count_zeros))
  sample <- unname(unlist(data[,min_zero_col]))
  mu <- mean(sample, na.rm=T)
  sd <- sd(sample, na.rm=T)
  data[na_zeros,grepl(condition,colnames(data))] <- na_zeros_mutate(na_zeros , mu,sd)
  return(data)
}

run_order <- rbind(cbind(rep("Ovary",2),c("Igg","Tumor")),
                   cbind(rep("Liver",3),c("Igg","Normal","Tumor")),
                   cbind(rep("Colon",3),c("Igg","Normal","Tumor")))

for(i in 1:dim(run_order)[1]){
  list_experiment[[run_order[i,1]]] <- mutate_na_zeros(list_experiment[[run_order[i,1]]],run_order[i,2])
}



#####################################################################################################
##  3. Impute values for proteins that have at least one non-zero intensity

### Imputate partial rows

imutate_partial_na <- function(data=list_experiment$Colon,condition="Tumor"){
  count_na <- apply(data[,grepl(condition,colnames(data))],1, count_zeros)
  y <- data[count_na>0,grepl(condition,colnames(data))]
  for(i in 1:dim(y)[1]){
  col_NAs <- names(y[i,is.na(unlist(y[i,]))])
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)))
      std <- sd(unlist(unname(delta)))
      cor <- cor(data[count_na==0,grepl(condition,colnames(data))])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,grepl(condition,colnames(data))] <- y
 return(data)
} 
  
  
for(i in 1:dim(run_order)[1]){
  list_experiment[[run_order[i,1]]] <- imutate_partial_na(list_experiment[[run_order[i,1]]],run_order[i,2])
}


####################################################################################################################
## Save the output in "../Result/3. After_imputation/"
####################################################################################################################


save(list_experiment,file="../Result/3. After_imputation/after_imputation.RData")

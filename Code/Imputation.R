## Mehrnoosh Oghbaie
## 12/11/2018
## Imputate missing values

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","plyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2","ggridges","gridExtra")
bioconductor.packages <- c("biomaRt")
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Imputation consists of three stages:
##  1. Remove proteins not identified in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity

load(file="../Result/2. After_preparation/log_transformed.RData")
set.seed(123)
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

list_experiment[["PreImputation"]] <- list_experiment
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
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}


rnorm(sum(na_zeros),mu - 1.8*sd,0.3*sd)
### function: mutate all NA values for specific condition( control or case)
mutate_na_zeros <- function(data=list_experiment$Liver,condition="Normal"){
  na_zeros <- rowSums(data[,grepl(condition,colnames(data))],na.rm=T)==0
  min_zero_col <- min_col(apply(data[,grepl(condition,colnames(data))],2,count_zeros))
  sample <- unname(unlist(data[,min_zero_col]))
  mu <- mean(sample, na.rm=T)
  sd <- sd(sample, na.rm=T)
  data[na_zeros,min_zero_col] <- na_zeros_mutate(na_zeros , mu,sd)
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

###############################################################################################
## Compare distribution of data before and after imputation.


plot_density<- function(data=list_experiment, tissue="Ovary"){
    g<-dim(data[[tissue]][,-(1:2)])
    my.env <- new.env()
    for(i in 1:g[2]){
      df <- data.frame(rbind(
        cbind(unlist(unname(data[["PreImputation"]][[tissue]][,-(1:2)][,i])),rep("Before_imputation", g[1])),
        cbind(unlist(unname(data[[tissue]][,-(1:2)][,i])),rep("After_imputation", g[1]))
      ))
      colnames(df) <- c("value","group")
      df$value<- as.numeric(as.character(df$value))
      mu <- ddply(df, "group", summarise, grp.mean=mean(value,na.rm=TRUE))
      assign(paste0("p",i),ggplot(df, aes(x=value, y=group, fill=factor(..quantile..))) +
               stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.25, 0.75)) +
               scale_fill_manual(
                 name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                 labels = c("(0, 0.25]", "(0.25, 0.75]", "(0.75, 1]")
               )+labs(title = paste("Log intensity",names(data[[tissue]][,-(1:2)][,i]))),envir = my.env )
      
    }
  return(my.env)
  }
  
  
my.env <- plot_density(data=list_experiment, tissue ="Colon")

p <- grid.arrange(grobs=list(my.env$p1,my.env$p2,my.env$p3,
                        my.env$p4,my.env$p5,my.env$p6,
                        my.env$p7,my.env$p8,my.env$p9),top = "Colon tissue")
p
 

my.env <-  plot_density(data=list_experiment, tissue ="Liver")
p2<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                 my.env$p4,my.env$p5,my.env$p6,
                 my.env$p7,my.env$p8,my.env$p9)
grid.arrange(p2,top = "Liver tissue")
p2

my.env <- plot_density(data=list_experiment, tissue ="Ovary")
p3<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                my.env$p4,my.env$p5,my.env$p6)
grid.arrange(p3,top = "Ovary tissue")
p3



##########################################################################
## Plot average value


draw_mean_impute <- function(data, tissue="Colon"){
  df <- data[[tissue]][,-(1:2)]
  df_average<- data.frame(matrix(NA, ncol=2, nrow=0))
  for (k in run_order[run_order[,1]==tissue,2]){
    df_average<- rbind(df_average,
                       cbind(apply(df[,grepl(k,colnames(df))],1,mean),rep(k,dim(df)[1])))
  }
  df_average$V1 <- as.numeric(as.character(df_average$V1))
  colnames(df_average) <- c("Log_intensity","Condition")
  
  p <- ggplot(df_average, aes(x = Log_intensity, y = Condition)) +
    geom_density_ridges(
      jittered_points = TRUE, quantile_lines = TRUE, scale = 0.6, alpha = 0.9,
      vline_size = 1.2, vline_color = "red", vline_alpha = 1,
      point_size = 0.8, point_alpha = 1,
      position = position_raincloud(adjust_vlines = TRUE)
    )+labs(title = paste("Average Log intensity",tissue))
  return(p)
}
q1 <- draw_mean_impute(data=list_experiment,tissue="Colon")
q2 <- draw_mean_impute(data=list_experiment,tissue="Liver")
q3 <- draw_mean_impute(data=list_experiment,tissue="Ovary")
q <- grid.arrange(grobs=list(q1,q2,q3),ncol=3,top="After imputation")
q

r1 <- draw_mean_impute(data=list_experiment$PreImputation ,tissue="Colon")
r2 <- draw_mean_impute(data=list_experiment$PreImputation ,tissue="Liver")
r3 <- draw_mean_impute(data=list_experiment$PreImputation , tissue="Ovary")
r <- grid.arrange(grobs=list(r1,r2,r3),ncol=3, top="Before imputation")
r
####################################################################################################################
## Save the output in "../Result/3. After_imputation/"
####################################################################################################################


save(KM050217_merged,list_experiment,file="../Result/3. After_imputation/after_imputation.RData")

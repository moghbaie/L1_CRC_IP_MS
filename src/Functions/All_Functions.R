# Mehrnoosh Oghbaie
# 12/10/2018
# Repository for all the functions
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  #if (length(bioconductor.packages) > 0) {
  #  source("http://bioconductor.org/biocLite.R")
  #}
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p,version = "3.8") 
    }
    library(p,lib.loc="~/R/win-library/3.5",character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
    }	
    #library(p, lib.loc="~/R/win-library/3.5",character.only=T) 
    library(p, lib.loc="~/R/win-library/3.5",character.only=T)
  }
}


###################################################################################
### Imputation functions
##################################################################################

### function : finding column with minimum zeros
count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

perseus_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(rnorm(sum(na_zeros), mean = mu-1.8*sd, sd = 0.3*sd),ncol=1))
}

#####################################################################################################
### Imputation method 6_2

#Method_2:
#  For each protein, which has zeros in each replica, replica with the least number of zeros is chosen.
#  For it mean and standard deviation of non-zero proteins intensities are counted. 
#  New intensity for this protein in this replica is sampled from uniform distribution with parameters: 
#  start = mu - 3*sd, end = mu - 2*sd. Then other replicas for this protein are imputed as in 6th step
#Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_6:
# Impute values for proteins, which have zero intensities only in some replicates:
#  Build distribution of deltas for all non zero proteins, where 
# delta =(Intrep1-Intrep2) mean(Intrep1,Intrep2)
# Calculate mudelta , sddelta
# Calculate new delta and new Intensity:
# deltanew=rnorm(mu=mudelta, sd=sddelta*sqrt(2)*mean(correlations))
# Inew=mean(Intother)*abs(1+deltanew)



### Method2
impute_all_zeros_2 <- function(data=list_experiment$Liver,condition="Normal"){
  condition <- paste0(condition, collapse="|")
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  na_zeros <- rowSums(data[,colnames],na.rm=T)==0
  min_zero_col <- min_col(apply(data[,colnames],2,count_zeros))
  sample <- unname(unlist(data[,min_zero_col]))
  mu <- mean(sample, na.rm=T)
  sd <- sd(sample, na.rm=T)
  data[na_zeros,colnames] <- cbind(na_zeros_impute(na_zeros , mu,sd),na_zeros_impute(na_zeros , mu,sd),na_zeros_impute(na_zeros , mu,sd))
  return(data)
}

### Funtion that Imputate partial rows

### Method6

impute_partial_zero_6_2 <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition,collapse="|")
  }
  data <- impute_all_zeros_2(data,condition)
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  count_na <- apply(data[,colnames],1, count_zeros)
  y <- data[count_na>0,colnames]
  print(paste("There are ",dim(y)[1]," records missing. \n And ", (dim(data)[1]-dim(y)[1]), " records with full values in ",condition))
  for(i in 1:dim(y)[1]){
    col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)), na.rm=T)
      std <- sd(unlist(unname(delta)),na.rm=T)
      cor <- cor(data[count_na==0,colnames])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,colnames] <- y
  return(data)
} 


#Method_1:
# For each replica mean and standard deviation of non-zero proteins intensities are counted. 
# New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
#  start = mu - 3*sd, end = mu - 2*sd
# Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_7:
# For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
# For other proteins uses method6

### function: mutate all NA values for specific condition( control or case)
### Method1
impute_all_zeros_1 <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  less1zero <- rowSums(data[,colnames],na.rm=T)<=1
  for(i in 1:length(colnames)){
    zero_col <- is.na(data[,colnames[i]])
    sample <- data[!less1zero,colnames[i]]
    sample <- sample[complete.cases(sample)]
    print(paste("There are ", length(sample), " complete records in ", colnames[i], ".\n There are ",sum(unname(less1zero&zero_col)), "zero records to fill"))
    #OutVals = boxplot(sample)$out
    #sample <- sample[!sample %in% OutVals]
    print(paste( length(sample), "records are used for imputing zero values in ", colnames[i]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[less1zero&zero_col,colnames[i]] <- na_zeros_impute(less1zero&zero_col, mu,sd)
    #data[less1zero&zero_col,colnames[i]] <- perseus_zeros_impute(less1zero&zero_col, mu,sd)
  }
  return(data)
}


### Method7_1
impute_partial_zero_7_1 <- function(data=list_experiment$Colon,condition="Tumor"){
  print(names(data))
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  count_na_initial <- apply(data[,colnames],1, count_zeros)
  data <- impute_all_zeros_1(data,condition)
  count_na <- apply(data[,colnames],1, count_zeros)
  y <- data[count_na>0,colnames]
  print(paste("There are ",dim(y)[1]," records missing. \n And ", (dim(data)[1]-dim(y)[1]), " records with full values in ",condition))
  for(i in 1:dim(y)[1]){
    col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)), na.rm=T)
      std <- sd(unlist(unname(delta)),na.rm=T)
      cor <- cor(data[count_na==0,colnames])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,colnames] <- y
  return(data)
} 

### Perseus imputation


impute_Perseus <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  for(i in 1:length(colnames)){
    zero_col <- is.na(data[,colnames[i]])
    sample <- data[,colnames[i]]
    sample <- sample[complete.cases(sample)]
    print(paste("There are ", length(sample), " complete records in ", colnames[i], ".\n There are ",sum(unname(zero_col)), "zero records to fill"))
    #OutVals = boxplot(sample)$out
    #sample <- sample[!sample %in% OutVals]
    print(paste(length(sample), "records are used for imputing zero values in ", colnames[i]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[zero_col,colnames[i]] <- perseus_zeros_impute(zero_col, mu,sd)
  }
  return(data)
}



### Function that plot intensity comparison before and after imputation

plot_density<- function(data , tissue="Ovary"){
  g<-dim(data[["PreImputation"]][[tissue]][,-(1:2)])
  my.env <- new.env()
  for(i in 1:g[2]){
    df <- data.frame(rbind(
      cbind(unlist(unname(data[["PreImputation"]][[tissue]][,-(1:2)][,i])),rep("Before_imputation", g[1])),
      cbind(unlist(unname(data[["experimentImputed"]][[tissue]][,-(1:2)][,i])),rep("After_imputation", g[1]))
    ))
    colnames(df) <- c("value","group")
    df$value<- as.numeric(as.character(df$value))
    mu <- plyr::ddply(df, "group", summarise, grp.mean=mean(value,na.rm=TRUE))
    assign(paste0("p",i),ggplot(df, aes(x=value, y=group, fill=factor(..quantile..))) +
             stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.25, 0.75)) +
             scale_fill_manual(
               name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
               labels = c("(0, 0.25]", "(0.25, 0.75]", "(0.75, 1]")
             )+labs(title = paste("Log intensity",names(data[[tissue]][,-(1:2)][,i]))),envir = my.env )
    
  }
  return(my.env)
}


### Function that plot average intensity  after imputation
draw_mean_impute <- function(data, tissue="Colon",run_order){

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


###################################################################################
### Anova functions
##################################################################################
#### Volcano plot

draw_volcanoplot <- function(data, condition, eLife_list){
  outliers <- CRC[["experimentPreImputateSignificant"]][[condition]]$Case>1
  ds <- data[[condition]]
  ds[["outliers"]] <- outliers
  ds <- ds[complete.cases(ds),]
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,";")[[1]][1]))
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,"-")[[1]][1]))
  fold_cutoff = 1
  pvalue_cutoff = 0.05
  gene_protein <- read.delim("protein.tab")
  gene_protein <- gene_protein %>% dplyr::mutate(Gene.names = as.character(Gene.names),
                                                 Entry = as.character(Entry))
  gene_protein$Gene.names  <- apply(gene_protein,1, function(x) strsplit(x[["Gene.names"]], " ")[[1]][1])
  ds$ID <- gene_protein$Gene.names[match(ds$uniprotID,gene_protein$Entry)]
  ds$ID <- ifelse(is.na(ds$ID),ds$uniprotID,ds$ID)
  ds$ID[ds$ID=="Q9UN81"] <- "L1RE1"
  col <- ifelse((ds$Significant == "Yes"& ds$outliers)&(ds$ID=="L1RE1"|grepl("ORF1P",ds$ID)), "red3", ifelse(!ds$outliers&(ds$Significant == "Yes")&(ds$ID=="L1RE1"|grepl("ORF1P",ds$ID)),"grey",ifelse(ds$Significant == "Yes"& ds$outliers,"black",  "grey")))
  group <- ifelse((ds$Significant == "Yes"& ds$outliers)&(ds$ID=="L1RE1"|grepl("ORF1P",ds$ID)), "Strong significant ORF1s",
                  ifelse(!ds$outliers&(ds$Significant == "Yes")&(ds$ID=="L1RE1"|grepl("ORF1P",ds$ID)),"Not significant",
                         ifelse(ds$Significant == "Yes"& ds$outliers,"Significant", "Not significant")))
  size <- ifelse((ds$uniprotID %in% eLife_list$uniprotID)&(ds$Significant == "Yes")&ds$outliers,5,1)
  colour <- ifelse((ds$Significant == "Yes"& ds$outliers)&(grepl("ORF1P", ds$ID)|ds$ID=="L1RE1"), "red3",ifelse((ds$Significant == "Yes"& !ds$outliers)&(grepl("ORF1P", ds$ID)|ds$ID=="L1RE1"),"grey","black") )
  lab <- ifelse((ds$Significant=="Yes"& ds$outliers)|ds$ID=="L1RE1",ds$ID,"")
  
  cols <- c("Strong significant ORF1s"="red3",  "Significant"="black",  "Not significant"="grey")
  
  pdf(paste0("../Image/Volcano_plot/Volcano_plot_",condition,".pdf"),width = 12, height = 18)
  p <- ggplot(ds,aes(logfold, -log10(p.adj),label=ID)) +
    geom_point(col= "orange",size=size)+
    geom_point(aes(colour=group),fill = col, size=2) +
    geom_vline(xintercept = fold_cutoff, col = "blue")+
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
    geom_text_repel(data  = subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))), 
                    colour = ifelse(grepl("L1RE1|ORF1P",subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID)))$ID),"red","black"), segment.size  = 0.2, 
                    segment.alpha =0.35,
                    #nudge_y =0.15,
                    box.padding = unit(0.55, "lines"), 
                    point.padding = unit(0.55, "lines"),
                    size=6
                    )+
    #geom_text_repel(data  = subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold<4.2)), 
    #                colour = ifelse(grepl("L1RE1|ORF1P",subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold<4.2))$ID),"red","black"), segment.size  = 0.2, 
    #                segment.alpha =0.3,
    #                nudge_x =  0.04 -subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold<4.2))$logfold, 
    #                box.padding = unit(0.55, "lines"), 
    #                point.padding = unit(0.55, "lines"),
    #                size=6
    #)+
    #geom_text_repel(data  = subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold>=1)), 
                    #colour = ifelse(grepl("L1RE1|ORF1P",subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold>=1))$ID),"red","black"), segment.size  = 0.2, 
                    #segment.alpha =0.3,
                    #nudge_x = 12- subset(ds, ((Significant=="Yes"& outliers)|ID=="L1RE1")&(!grepl("ORF1P",ID))&(logfold>=1))$logfold, 
                    #box.padding = unit(1, "lines"), 
                    #point.padding = unit(1, "lines"),
                    #segment.size = 0.2,
                   # segment.color = "grey50",
                    #size=6, 
                    #direction="x"
    #)+
    ggtitle(condition)+ scale_colour_manual(values=cols, aesthetics = c("fill","colour"))+
    theme_minimal()+ 
   
    #scale_x_continuous(limits = c(-3, 12))+ #144
    #scale_y_continuous(limits = c(0, 3.8))
    #scale_x_continuous(limits = c(-5, 11))+ #159
    #scale_y_continuous(limits = c(0, 4.5))
    #scale_x_continuous(limits = c(-7, 9.5))+ #163
    #scale_y_continuous(limits = c(0, 3.8))
    scale_x_continuous(limits = c(-7, 12))+
    scale_y_continuous(limits = c(0, 4.5))
  
  print(p)
  dev.off()
}



#######################################################################################
#### Normalize functions
######################################################################################

numerize<- function(x) ifelse(x=="Yes",1,0)

#######################################################################################
#### Normalize before imputation
#######################################################################################

normalize <- function(dx){
  x1 <- unlist(unname(apply(dx[,-1:-2],2, function(x) max(x, na.rm=T))))
  x2 <- unlist(unname(dx[dx$uniprotID=="Q9UN81",-1:-2]))
  x <- ifelse(!is.na(x2),x2,x1)
  y <- t(t(dx[,-1:-2])/x)
  return(y)
}


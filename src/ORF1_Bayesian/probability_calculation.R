library(cleaver)
library(ape)
library(seqinr)

orf1 <- readAAStringSet("F:/Line1/one_file_per_ptne/orf1.fasta")
names(orf1)
lapply(names(orf1), function(x) strsplit(x, "\\|")[[1]][3])

cleave(orf1[names(orf1)==names(orf1)[1]], enzym="trypsin")

ortho1 <- names(orf1)[unlist(lapply(names(orf1), function(x) strsplit(x, "\\|")[[1]][2])) %in% Unique_protein]

cleave(x, enzym="trypsin")

peptides <- read.delim("F:/Line1_CRC/Input_data/Third_Run_03142019/txt/allPeptides.txt", header=FALSE)
dim(peptides)
pep <- peptides[grepl("ORF1P|Q9UN81|L1RE1",peptides$V36),c("V1","V2","V3","V34","V35","V36","V37","V38")]
cpep <- pep[,c("V1","V35")]
cpep$V1 <- as.character(cpep$V1)
cpep$V35 <- as.character(cpep$V35)

cpep2 <- data.frame(do.call("rbind", apply(cpep, 1,function(x) cbind(rep(x[1], length(strsplit(x[2], ";|,")[[1]])), strsplit(x[2], ";")[[1]]))))
cpep2$X1 <- as.character(cpep2$X1)
cpep2$X2 <- as.character(cpep2$X2)
cpep2 <- cpep2[!(grepl("L1_ORF1",cpep2$X2)|(!grepl("ORF1|Q9UN81",cpep2$X2))),]

cpep3 <- dcast(cpep2, X1~X2)
rownames(cpep3) <-cpep3$X1
cpep3 <- cpep3[,-1]
cpep3[!is.na(cpep3)] <- 1.0
cpep3[is.na(cpep3)] <- 0.0
for(i in 1: dim(cpep3)[2]){
  cpep3[,i] <- as.numeric(cpep3[,i])
}

###Remove LIHS-ORF1 peptides from the matrix
LIHS_peps <- rownames(cpep3)[which(cpep3[,"Q9UN81"]==1)]
### Remove LORF1 peptides
cpep4 <- cpep3[!rownames(cpep3) %in% LIHS_peps, colnames(cpep3)!= "Q9UN81"]


dM <- t(Matrix(matrix(t(cpep4),nrow=99, ncol=93)))
str(cpep4)
colnames(dM) <- colnames(cpep4)
rownames(dM) <- rownames(cpep4)
dM2 <- dM[,colnames(dM) %in% names(colSums(dM)[colSums(dM)!=0])]

colSums(dM2)
dd <- t(t(dM2)/colSums(dM2))

################################################################
peptides1 <- read.delim("F:/Line1_CRC/Input_data/Third_Run_03142019/txt/peptides.txt", header=TRUE)
dim(peptides1)
selected_col  <- colnames(peptides1)[grepl("Sequence|Protein|PEP|LFQ.intensity|Start.position|End.position|Amino.acid.before|Amino.acid.after|Missed.cleavage",colnames(peptides1))]
pe <- peptides1[grepl("ORF1P|Q9UN81|L1RE1",as.character(peptides1$Proteins)),c(selected_col) ]
cpe <- pe[pe$Missed.cleavages==0&(!pe$Sequence %in% LIHS_peps),]

run_order2 <- rbind(cbind(rep("Ovary",4),c("144T_IgG_2","144T_IgG_10","144T_ORF1_1","144T_ORF1_9")),
                    cbind(rep("Liver",3),c("159T_IgG_4","159N_ORF1_5","159T_ORF1_3")),
                    cbind(rep("Colon",3),c("163T_IgG_7","163N_ORF1_8","163T_ORF1_6")))

run_order2 <- cbind(run_order2,c(1,2,1,2,1,1,1,1,1,1))

con <- list()
con[["Ovary"]] <- cpe[,selected_col[grepl("144|Sequence|Protein|PEP",selected_col)]]
con[["Liver"]] <- cpe[,selected_col[grepl("159|Sequence|Protein|PEP",selected_col)]]
con[["Colon"]] <- cpe[,selected_col[grepl("163|Sequence|Protein|PEP",selected_col)]]


##Calculate list of significant non-consensus peptides.

anovaAnalysis <- function(run_order2,con){
  anova_result <- list()
  for( i in unique(run_order2[,1])){
    y <- con[[i]]
    x <- unique(run_order2[run_order2[,1]==i,3])
    for(l in x){
      z <- unique(run_order2[(run_order2[,1]==i)&(run_order2[,3]==l),2])
      print(z)
      for(j in z[!grepl("T_ORF1",z)]){
        dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
        dt[["Proteins"]] <- y[["Proteins"]]
        dt[["Sequence"]] <- y[["Sequence"]]
        dt[["PEP"]] <- y[["PEP"]]
        dt[["p.value"]] <- NA
        dt[["logfold"]] <- NA
        dt[["averageCase"]] <- NA
        dt[["averageControl"]] <- NA
        dt[["Significant"]] <- NA
        dt[["Control"]] <- apply(y[,grepl(z[!grepl("T_ORF1",z)],colnames(y))],1, function(x) sum(!is.na(unname(unlist(x)))))
        dt[["Case"]] <- apply(y[,grepl(z[grepl("T_ORF1",z)],colnames(y))],1, function(x) sum(!is.na(unname(unlist(x)))))
        #dt <- dt %>% dplyr::filter(!uniprotID %in% self[["non_unique_ORF"]])
        
        w <- y[,grepl(paste0(c(j,z[grepl("T_ORF1",z)]),collapse="|"),colnames(y))]
        
        w[is.na(w)] <- 0
        for(k in 1:dim(dt)[1]){
          print(k)
          dt[[k,"p.value"]] <- t.test(unname(unlist(w[k,grepl("T_ORF1",colnames(w))])),unname(unlist(w[k,grepl(j,colnames(w))])))$p.value
          dt[[k,"averageCase"]] <- mean(unname(unlist(w[k,grepl("T_ORF1",colnames(w))])))
          dt[[k,"averageControl"]] <- mean(unname(unlist(w[k,grepl(j,colnames(w))])))
          dt[[k,"logfold"]] <- mean(unname(unlist(w[k,grepl("T_ORF1",colnames(w))])))-mean(unname(unlist(w[k,grepl(j,colnames(w))])))
        }
        
        dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
        dt[["Significant"]] <- ifelse((dt[["p.adj"]]<0.05|dt[["averageControl"]]==0)&dt[["logfold"]]>1,"Yes","No")
        dt <- dt[dt$Significant=="Yes",]
        anova_result[[paste0(i,"_","Tumor","_",j)]] <-dt
      }
    }
  }
  return(anova_result)
}


anova_result <- anovaAnalysis(run_order2,con)

#### Calculating the percentage of captured variability
Var <- list()
for(name in names(anova_result)){
  dk <- anova_result[[name]]
  Proteins <- unique(unlist(lapply(dk$Proteins, function(x) strsplit(as.character(x), ";")[[1]])))
  
  df <- data.frame(Protein= Proteins)
  df$Protein <- as.character(df$Protein)
  df[["Sequences"]] <- NA
  #df[["PEP"]] <- NA
  df[["Captured_Variability_ratio"]] <- NA
  df[["Num_changes_consensus"]] <- NA
  for(j in 1:length(Proteins)){
    df[j,"Protein"]
    seq <- dk %>% dplyr::filter(grepl(paste0(df[j,"Protein"],"\\b"), as.character(dk$Proteins)))%>% dplyr::select(Sequence) %>% .$Sequence
    intensity <- dk %>% dplyr::filter(grepl(paste0(df[j,"Protein"],"\\b"), as.character(dk$Proteins)))%>% dplyr::select(Sequence, averageCase)
    
    df[j, "Captured_Variability_ratio"] <- round(sum(dd[rownames(dd) %in% as.character(seq),colnames(dd)==df[j,"Protein"]]),2)
    
    df[j, "Num_changes_consensus"] <- length(unname(unlist(dd[dd[,colnames(dd)==df[j,"Protein"]]!=0 ,colnames(dd)==df[j,"Protein"]])))
    df[j,"Sequences"] <- paste(as.character(seq), collapse="|")
    df[j,"Pep_Intensity"] <- paste(apply(intensity, 1,function(x) paste0(c(x[["Sequence"]],x[["averageCase"]]), collapse="|")),collapse=";")
  }
  Var[[name]] <- df
}

for(name in names(Var)){
  dy <- anova_result[[name]]
  dy[["Proteins_Percentage_captured"]] <- NA
  dy[["Num_changes_consensus"]] <- NA
  for(j in 1:dim(dy)[1]){
    prots <- unlist(strsplit(as.character(dy[j,"Proteins"]),";")[[1]])
    
    varProt <- Var[[name]][Var[[name]]$Protein %in% prots,c("Protein","Captured_Variability_ratio")]
    varProt <- varProt[order(-varProt$Captured_Variability_ratio),]
    
    varPossible <- Var[[name]][Var[[name]]$Protein %in% prots,c("Protein","Num_changes_consensus","Captured_Variability_ratio")]
    varPossible <- varPossible[order(-varPossible$Captured_Variability_ratio),]
    dy[j,"Proteins_Percentage_captured"] <- paste(apply(varProt, 1,function(x) paste(c(x[["Protein"]],x[["Captured_Variability_ratio"]]), collapse="|")),collapse=";")
    dy[j,"Num_changes_consensus"] <- paste(apply(varPossible, 1,function(x) paste(c(x[["Protein"]],x[["Num_changes_consensus"]]), collapse="|")),collapse=";")
  }
  anova_result[[name]] <- dy
}


proteinGroups <- read.delim("F:/Line1_CRC/Input_data/Third_Run_03142019/txt/proteinGroups.txt")

colnames(proteinGroups)
#View(proteinGroups[grepl("ORF1|Q9UN81",proteinGroups$Protein.IDs),])
orf1_significant <- lapply(as.character(c(109,60,118,111,77,110,122,60,55,56,102,35,98,76,45,29,79,18,69,96,13,46,9)),function(x) paste0("ORF1P.",x))
dd <- proteinGroups[grepl("ORF1|Q9UN81",proteinGroups$Protein.IDs),c("Protein.IDs","Gene.names","Fasta.headers","Peptide.counts..all.","Peptides","Peptide.counts..unique.","Peptide.counts..razor.unique.")][unlist(lapply(lapply(proteinGroups[grepl("ORF1|Q9UN81",proteinGroups$Protein.IDs),c("Protein.IDs")],function(x) strsplit(as.character(x),";")[[1]]), function(x) any(x %in% orf1_significant))),]


xlsx.writeMultipleData <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
}


c(names(Var))
file <- "C:/Users/ProteHomics/Desktop/protein_var.xlsx"
write.xlsx(Var$Colon_Tumor_163N_ORF1_8, file, sheetName = "Colon_Tumor_163N_ORF1_8",append=TRUE)

xlsx.writeMultipleData("C:/Users/ProteHomics/Desktop/protein_var.xlsx",
                       Var$Ovary_Tumor_144T_IgG_2,
                       Var$Ovary_Tumor_144T_IgG_10,
                       Var$Liver_Tumor_159T_IgG_4,
                       Var$Liver_Tumor_159N_ORF1_5,
                       Var$Colon_Tumor_163T_IgG_7,
                       Var$Colon_Tumor_163N_ORF1_8)

msms <- read.delim("F:/Line1_CRC/Input_data/Fourth_Run_04202019/txt/msms.txt")
write.csv(msms[grepl("ORF1P|Q9UN81",msms$Proteins),], file="C:/Users/ProteHomics/Desktop/msms_orf1_orthogonal.csv")

msms1 <- read.delim("F:/Line1_CRC/Input_data/Third_Run_03142019/txt/msms.txt")
write.csv(msms1[grepl("ORF1P|Q9UN81",msms1$Proteins),], file="C:/Users/ProteHomics/Desktop/msms_orf1.csv")

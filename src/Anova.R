## Mehrnoosh Oghbaie
## 12/11/2018
## ANOVA test between cases and controls

##  1. Filter proteins with less than two peptides (I'll skip this one for now)
##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)
Template$set("public","experimentPreImputateSignificant", list())
Template$set("public","experimentImputateSignificant", list())
Template$set("public","experimentImputateSignificantPerseus", list())



Template$set("public","anovaAnalysisPerseus", function(run_order2){
  for( i in unique(run_order2[,1])){
    y <- self[["experimentImputedPerseus"]][[i]]
      x <- unique(run_order2[run_order2[,1]==i,3])
      for(l in x){
        z <- unique(run_order2[(run_order2[,1]==i)&(run_order2[,3]==l),2])
        print(z)
        for(j in z[!grepl("T_ORF1",z)]){
          dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
          dt[["ID"]] <- y[["id"]]
          dt[["uniprotID"]] <- y[["uniprotID"]]
          dt[["p.value"]] <- NA
          dt[["logfold"]] <- NA
          dt[["Significant"]] <- NA
          dt <- dt %>% dplyr::filter(!uniprotID %in% self[["non_unique_ORF"]])
          for(k in 1:dim(dt)[1]){
            print(k)
            dt[[k,"p.value"]] <- t.test(unname(unlist(y[k,grepl(z[grepl("T_ORF1",z)],colnames(y))])),unname(unlist(y[k,grepl(j,colnames(y))])))$p.value
            dt[[k,"logfold"]] <- mean(unname(unlist(y[k,grepl(z[grepl("T_ORF1",z)],colnames(y))])),rm.na=T)-mean(unname(unlist(y[k,grepl(j,colnames(y))])),rm.na=T)
          }
          dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
          dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&dt[["logfold"]]>1,"Yes","No")
          dt[["Significant_fdr_0.05"]] <- qvalue(p = dt[["p.value"]],lfdr=TRUE, fdr.level = 0.05)$significant
          dt[["Significant_fdr_0.05_logfold"]] <- ifelse(dt[["Significant_fdr_0.05"]]&dt[["logfold"]]>1,"Yes","No")
          dt[["lfdr"]] <- lfdr(p = dt[["p.value"]])
          self[["experimentImputateSignificantPerseus"]][[paste0(i,"_","Tumor","_",j)]] <-dt
        }
      }
  }
})


Template$set("public","anovaAnalysisPreImpute", function(run_order2){
  for( i in unique(run_order2[,1])){
    y <- self[["PreImputation"]][[i]]
    x <- unique(run_order2[run_order2[,1]==i,3])
    for(l in x){
      z <- unique(run_order2[(run_order2[,1]==i)&(run_order2[,3]==l),2])
      print(z)
      for(j in z[!grepl("T_ORF1",z)]){
        dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
        dt[["ID"]] <- y[["id"]]
        dt[["uniprotID"]] <- y[["uniprotID"]]
        dt[["p.value"]] <- NA
        dt[["logfold"]] <- NA
        dt[["Significant"]] <- NA
        dt[["Control"]] <- apply(y[,grepl(z[!grepl("T_ORF1",z)],colnames(y))],1, function(x) sum(!is.na(unname(unlist(x)))))
        dt[["Case"]] <- apply(y[,grepl(z[grepl("T_ORF1",z)],colnames(y))],1, function(x) sum(!is.na(unname(unlist(x)))))
        dt <- dt %>% dplyr::filter(!uniprotID %in% self[["non_unique_ORF"]])
        
        w <- y[,grepl(paste0(c(j,z[grepl("T_ORF1",z)]),collapse="|"),colnames(y))]
     
        w[is.na(w)] <- 0
        for(k in 1:dim(dt)[1]){
          print(k)
            dt[[k,"p.value"]] <- t.test(unname(unlist(w[k,grepl("T_ORF1",colnames(w))])),unname(unlist(w[k,grepl(j,colnames(w))])))$p.value
            dt[[k,"logfold"]] <- mean(unname(unlist(w[k,grepl("T_ORF1",colnames(w))])))-mean(unname(unlist(w[k,grepl(j,colnames(w))])))
        }
       
        dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
        dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&dt[["logfold"]]>1,"Yes","No")
        dt[!is.na(dt[["p.value"]]),"Significant_fdr_0.05"] <- qvalue(p = dt[["p.value"]][!is.na(dt[["p.value"]])], lfdr=TRUE,fdr.level = 0.05)$significant
        dt[["Significant_fdr_0.05_logfold"]] <- ifelse(dt[["Significant_fdr_0.05"]]&dt[["logfold"]]>1,"Yes","No")
        dt[["lfdr"]] <- lfdr(p = dt[["p.value"]])
        self[["experimentPreImputateSignificant"]][[paste0(i,"_","Tumor","_",j)]] <-dt
      }
    }
  }
})


Template$set("public","anovaAnalysis", function(run_order2){
for( i in unique(run_order2[,1])){
  y <- self[["experimentImputed"]][[i]]
    x <- unique(run_order2[run_order2[,1]==i,3])
    for(l in x){
      z <- unique(run_order2[(run_order2[,1]==i)&(run_order2[,3]==l),2])
      print(z)
      for(j in z[!grepl("T_ORF1",z)]){
        dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
        dt[["ID"]] <- y[["id"]]
        dt[["uniprotID"]] <- y[["uniprotID"]]
        dt[["p.value"]] <- NA
        dt[["logfold"]] <- NA
        dt[["Significant"]] <- NA
        dt <- dt %>% dplyr::filter(!uniprotID %in% self[["non_unique_ORF"]])
        for(k in 1:dim(dt)[1]){
            print(k)
            dt[[k,"p.value"]] <- t.test(unname(unlist(y[k,grepl(z[grepl("T_ORF1",z)],colnames(y))])),unname(unlist(y[k,grepl(j,colnames(y))])))$p.value
            dt[[k,"logfold"]] <- mean(unname(unlist(y[k,grepl(z[grepl("T_ORF1",z)],colnames(y))])),rm.na=T)-mean(unname(unlist(y[k,grepl(j,colnames(y))])),rm.na=T)
        }
       
        dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
        dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&dt[["logfold"]]>1,"Yes","No")
        dt[["Significant_fdr_0.05"]] <- qvalue(p = dt[["p.value"]],lfdr=TRUE,  fdr.level = 0.05)$significant
        dt[["Significant_fdr_0.05_logfold"]] <- ifelse(dt[["Significant_fdr_0.05"]]&dt[["logfold"]]>1,"Yes","No")
        dt[["lfdr"]] <- lfdr(p = dt[["p.value"]])
        
        self[["experimentImputateSignificant"]][[paste0(i,"_","Tumor","_",j)]] <-dt
      }
    }

  }
})


####################################################################################################################
## Venn diagram
Template$set("public","Significant_list", NA)
Template$set("public","Significant_proteins", NA)

Template$set("public","drawVenDiagramSignificant", function(){
  Significant_proteins <- data.frame(matrix(NA, nrow=length(self[["input_merged"]]$uniprotID), ncol= 2+length( names(self[["experimentImputateSignificant"]]))))
  colnames(Significant_proteins) <- c("id","uniprotID", names(self[["experimentImputateSignificant"]]))
  Significant_proteins$uniprotID <- self[["input_merged"]]$uniprotID
  #Significant_proteins$id <- self[["input_merged"]]$id
  
  for(i in names(self[["experimentImputateSignificant"]])){
    case1more <- self[["experimentPreImputateSignificant"]][[i]][["Case"]][match(Significant_proteins$uniprotID,self[["experimentPreImputateSignificant"]][[i]][,"uniprotID"])]
    Significant_proteins[[i]] <- ifelse(case1more>1, 
                                        self[["experimentImputateSignificant"]][[i]][,"Significant"][match(Significant_proteins$uniprotID,self[["experimentImputateSignificant"]][[i]][,"uniprotID"])],
                                        "No")
    
    Significant_proteins[[i]][is.na(Significant_proteins[[i]])] <- "NO"
    Significant_proteins[[i]] <- ifelse(Significant_proteins[[i]]=="Yes"&(self[["experimentPreImputateSignificant"]][[i]][,"Case"]>1),1,0)
  }
self[["Significant_proteins"]] <- Significant_proteins

Significant_list <- Significant_proteins%>%
  dplyr::filter(.[[3]]|.[[4]]|.[[5]]|.[[6]]|.[[7]]|.[[8]]==1)%>%
  dplyr::select(id,uniprotID)
self[["Significant_list"]] <- Significant_list
png(paste0("../Image/Volcano_plot/Significant_Venndiagram.png"),width = 800, height = 600)
venn(Significant_proteins[,-c(1:2)],zcolor = "style")
dev.off()
})

####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################
Template$set("public","drawVolcanoPlot", function(){
  for(condition in names(self[["experimentImputed"]])[grep("_",names(self[["experimentImputed"]]))]){
    draw_volcanoplot(data=self[["experimentImputed"]],condition=condition)
    Sys.sleep(2)
  }
})




## Mehrnoosh Oghbaie
## 05/21/2019
## Draw mutation from consensus in ORF1s sequences


Template$set("public","drawMutation", function(){
  detected_orf1 <- CRC$Significant_list$uniprotID[grepl("ORF1P",CRC$Significant_list$uniprotID)]
  
  
  peptides <- read.delim(paste0("../Input_data/",strsplit(as.character(CRC$input$V2[1]),"/")[[1]][1],"/txt/peptides.txt"), header=FALSE)
  pep <- peptides[peptides$V36 %in% detected_orf1,c("V1","V2","V3","V34","V35","V36","V37","V38")]

  orf1_dir <- "../Input_data/orf1.fasta"
  mySequences <- readAAStringSet(orf1_dir)
  name_selected <- names(mySequences)[unlist(lapply(names(mySequences), function(x) strsplit(x,"\\|")[[1]][2])) %in% detected_orf1]
  myseq <- mySequences[name_selected]
  
  myseq1 <- myseq
  names(myseq1) <- paste(lapply(names(myseq1), function(x) strsplit(x,"\\|")[[1]][2]),lapply(names(myseq1), function(x) strsplit(x,"\\.")[[1]][5]), sep="-")
  
  orf1_align <- msa(myseq1)
  orf1_align2 <- msaConvert(orf1_align, type="seqinr::alignment")
  
  d <- dist.alignment(orf1_align2, "identity")
  orf1Tree <- nj(d)
  png(filename="../Image/Phylogenetic_tree/tree.png",width = 1000, height = 900)
  plot(orf1Tree, main="ORF1_mutations")
  dev.off()
  
  names(myseq) <- unlist(lapply(names(myseq), function(x) strsplit(x,"\\|")[[1]][2]))
  
  for(name in names(CRC$experimentImputateSignificant)){
    print(name)
   
    case1more <- CRC[["experimentPreImputateSignificant"]][[name]][["Case"]][match(CRC$experimentImputateSignificant[[name]][["uniprotID"]],CRC[["experimentPreImputateSignificant"]][[name]][,"uniprotID"])]
    significant_name <- CRC$experimentImputateSignificant[[name]][["uniprotID"]][CRC$experimentImputateSignificant[[name]][["Significant"]]=="Yes"&case1more>1]

    significant_name_ORF <- significant_name[grepl("ORF1",significant_name)]
    myseq_Tumor_name <- myseq[names(myseq) %in% significant_name_ORF]
    if(length(myseq_Tumor_name)>1){
      orf1_align_name <- msa(myseq_Tumor_name)
      orf1_align_name2 <- msaConvert(orf1_align_name, type="seqinr::alignment")
      
      aln = data.frame(
        letter=strsplit(paste0(orf1_align_name2$seq, collapse=""), "")[[1]], 
        mutation = rep(orf1_align_name2$nam, each=338),
        x       = rep(1:338, orf1_align_name2$nb)
      )
      
      pep_orf1 <- pep[pep$V36 %in% orf1_align_name2$nam,]
      aln$mut = 'no'
      for(i in as.character(pep_orf1$V1)){
        aln[(aln$x %in% c(as.numeric(as.character(pep_orf1[as.character(pep_orf1$V1)==i, "V37"])):as.numeric(as.character(pep_orf1[as.character(pep_orf1$V1)==i, "V38"]))))& as.character(aln$mutation)== as.character(pep_orf1[as.character(pep_orf1$V1)==i, "V36"]),"mut"] <- "yes"
      }
      
      aln$consensus <- "yes"
      Lorf1 <- strsplit("MGKKQNRKTGNSKTQSASPPPKERSSSPATEQSWMENDFDELREEGFRRSNYSELREDIQTKGKEVENFEKNLEECITRITNTEKCLKELMELKTKARELREECRSLRSRCDQLEERVSAMEDEMNEMKREGKFREKRIKRNEQSLQEIWDYVKRPNLRLIGVPESDVENGTKLENTLQDIIQENFPNLARQANVQIQEIQRTPQRYSSRRATPRHIIVRFTKVEMKEKMLRAAREKGRVTLKGKPIRLTADLSAETLQARREWGPIFNILKEKNFQPRISYPAKLSFISEGEIKYFIDKQMLRDFVTTRPALKELLKEALNMERNNRYQPLQNHAKM","")[[1]]
      aln$LORF1 <- NA
      for( i in as.character(aln$mutation)){
        aln$LORF1[aln$mutation==i] <- Lorf1[aln$x[aln$mutation==i]]
      }
      aln$consensus <- ifelse(as.character(aln$letter)==as.character(aln$LORF1),"yes","no")
      
      dt <- data.frame(cbind(Lorf1,rep("L1RE1",338),c(1:338),rep("yes",338),rep("yes",338),Lorf1))
      colnames(dt) <- colnames(aln)
      dt$x <- as.integer(as.character(dt$x))
      dt$mut <- as.character(dt$mut)
      dt$consensus <- as.character(dt$consensus)
      dt$LORF1 <- as.character(dt$LORF1)
      
      aln <- rbind(aln,dt)
      #unique(aln$mut)
      
      png(filename=paste0("../Image/Phylogenetic_tree/",name,".png"),width = 2000, height =(orf1_align_name2$nb+1)*30)
      # Generate the sequence alignment
      p2 = ggplot(aln, aes(x, mutation)) +
        geom_text(aes(label=letter, color=consensus, size=mut)) + 
        scale_x_continuous(breaks=1:10, expand = c(0.0105, 0)) + xlab('') + 
        scale_color_manual(values=c('red', "black")) + 
        scale_size_manual(values=c(1.6, 4)) + 
        theme_logo() + 
        theme(legend.position = 'none', axis.text.x = element_blank())
      
      print(p2)
      dev.off()
    }
  }
})
             
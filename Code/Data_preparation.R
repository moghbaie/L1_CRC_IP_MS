## Mehrnoosh Oghbaie
## 12/11/2018
## Preparing data from MaxQuant or ...

## Data preparation consists of four stages:
##  1. Calculate iBAQ intensity
##  2. Remove contaminants and reverse proteins
##  3. Log transformation
##  4. Separate different experiment


Template <- R6Class("Template",
                   private = list(
                     input.info.dir = "../Input_data/Input.info",
                     df = NA,
                     df_log=NA
                   ),
                   active = list(
                     changeInfoFile = function(new_info_file){
                       if(missing(new_info_file)) return(private$input.info.dir)
                       else private$input.info.dir <- new_info_file
                     },
                     input_dir = function() {
                       dirname(private$input.info.dir)
                     }, 
                     input = function(){
                       read.delim(private$input.info.dir, header=FALSE)
                     },
                     data_type = function(){
                       self$input %>%
                         dplyr::filter(V1=="data_type") %>%
                         dplyr::mutate(V2=as.character(V2))%>% 
                         dplyr::select(V2) %>% 
                         .$V2
                     }, 
                     input_merged = function(){
                       tbl <- read_excel(paste(self$input_dir,self$input %>%
                                          dplyr::filter(V1=="MaxQuant_Output") %>%
                                          dplyr::mutate(V2=as.character(V2))%>%
                                          dplyr::select(V2) %>%
                                          .$V2, sep="/"))
                       
                       tbl$uniprotID = unlist(lapply(tbl$`Fasta headers`, function(x) strsplit(as.character(x),"\\|")[[1]][2]))
                       return(tbl)
                     },
                     contaminant_list =function(){
                       contaminants_tbl <- read.table(paste(self$input_dir,
                                                                       self$input %>%
                                                                         dplyr::filter(V1=="Contaminants") %>%
                                                                         dplyr::mutate(V2=as.character(V2))%>% 
                                                                         dplyr::select(V2) %>% 
                                                                         .$V2, sep="/"),
                                                                 sep=";", quote="\"")
                       contaminant_list <- unlist(lapply(as.character(contaminants_tbl[,1])[grepl(">", as.character(contaminants_tbl[,1]))], function(x) strsplit(x," |>")[[1]][2]))
                       return(contaminant_list)
                       } 
                      )
                   )


Template$set("public","removeContaminant", function(){
  cols <- colnames(self$input_merged)
  private$df <- self$input_merged[!self$input_merged$uniprotID %in% self$contaminant_list,c("id","uniprotID",cols[grepl(self$data_type,cols)][-1])]
}
)

Template$set("public","logTransformation", function(){
  tp <- lapply(private$df, class)
  df_log <- private$df
  for(i in colnames(private$df)[tp =="numeric"]){
    df_log[i] <- log(unlist(unname(df_log[i])))
    df_log[is.infinite(unname(unlist(df_log[i]))),i] <- NA
  }
  colnames(df_log) <- c("id","uniprotID",unlist(lapply(colnames(df_log)[-(1:2)], function(x) strsplit(x,"_")[[1]][2])))
  private$df_log <- df_log
}
)

Template$set("public","experiment", list())
                      
Template$set("public","separatedList", function(){
  self$experiment[["Ovary"]] <- private$df_log[, grepl("1|2|uniprotID|id",colnames(private$df_log))]
  colnames(self$experiment[["Ovary"]]) <- paste0(colnames(self$experiment[["Ovary"]]),c(rep("",2),rep("_Tumor",3),rep("_Igg",3)))
  
  
  self$experiment[["Liver"]] <- private$df_log[, grepl("3|4|5|uniprotID|id",colnames(private$df_log))]
  colnames(self$experiment[["Liver"]]) <- paste0(colnames(self$experiment[["Liver"]]),c(rep("",2),rep("_Normal",3),rep("_Tumor",3),rep("_Igg",3)))
  
  self$experiment[["Colon"]] <- private$df_log[, grepl("6|7|8|uniprotID|id",colnames(private$df_log))]
  colnames(self$experiment[["Colon"]]) <- paste0(colnames(self$experiment[["Colon"]]),c(rep("",2),rep("_Normal",3),rep("_Tumor",3),rep("_Igg",3)))
}
)


## Mehrnoosh Oghbaie
## 12/11/2018
## Quality control test

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","yaml","reshape2","dplyr","magrittr","sqldf","PTXQC","corrplot","ggplot2","methods")
bioconductor.packages <- c("biomaRt")
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## QC steps


####################################################################################################################
## QC plot "../Image/QC_plot/"
####################################################################################################################

txt_folder <-  file.path(dirname(getwd()),"Input_data","Fourth_Run_04202019/txt")

r = createReport(txt_folder)

cat(paste0("\nReport generated as '", r$report_file, "'\n\n"))


####################################################################################################################
## Save the output in "../Result/1. QC/"
####################################################################################################################

fh_out = getReportFilenames(txt_folder)
if (file.exists(fh_out$yaml_file))
{
  cat("\nUsing YAML config already present in target directory ...\n")
  yaml_config = yaml.load_file(input = fh_out$yaml_file)
} else {
  cat("\nYAML config not found in folder '", txt_folder, "'. The first run of PTXQC will create one for you.", sep="")
  yaml_config = list()
}

r = createReport(txt_folder, yaml_config)

cat(paste0("\nReport generated as '", r$report_file, "'\n\n"))

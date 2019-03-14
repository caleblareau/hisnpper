options(warn=-1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# This function / script for processed chromosome files to produce consensus barcode doublets
# as well as summary and QC metrics and visuals

# TO DO:
# actually call doublets -- make another additional text file
# produce QC plot
# make barcode key/value

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

outdir <- args[i+1] # directory of .rds files
file <- args[i+2] # file path to the number of barcodes for each observed barcode

importDT <- function(file){
  if(tools::file_ext(file) == "gz"){
    cov <- data.table::fread(cmd = paste0("zcat < ", file), stringsAsFactors = FALSE, header = TRUE)
  } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
    cov <- data.table::fread(paste0(file), stringsAsFactors = FALSE, header = TRUE)
  } else{
    stop("Provide a valid file format for the  file (.gz, .txt, .csv, or .tsv)")
  }
}

df <- importDT(file)
#df <- df[df[,3] != df[,4],]
listy <- split( df , f = df[,1])
lapply(1:length(listy), function(i){
  chr <- names(listy)[i]
  write.table(listy[[i]],
            file = paste0(outdir, "/", "SNPs_", chr, ".tsv"),
           quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  chr
}) -> token

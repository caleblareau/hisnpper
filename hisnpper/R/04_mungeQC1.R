options(warn=-1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(stringr)))

"%ni%" <- Negate("%in%")

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

directory <- args[i+1] # file path for the files
outfile <- args[i+2] # final resting places

if(FALSE){
  directory <- "/data/aryee/caleb/mESC_allele/biorad-bams/hs_Exp75-S15/logs/readstats/"
  outfile <- "o.txt"
}

lf <- list.files(directory, pattern = "^quant", full.names = TRUE)

lapply(lf, read.table, col.names = c("group", "count")) %>% rbindlist() %>% data.frame() -> wholeDF
wholeDF %>% group_by(group) %>% summarise(total = sum(count)) %>% data.frame()-> sumDF

write.table(sumDF, file = outfile, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

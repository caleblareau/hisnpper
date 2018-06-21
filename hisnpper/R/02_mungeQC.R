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

outfile <- args[i+1] # final resting places
directory <- args[i+2] # file path for the files

lf <- list.files(directory, pattern = "_stats.txt$", full.names = TRUE)

lapply(lf, read.table) %>% rbindlist() %>% data.frame() -> wholeDF
splitDF <- data.frame(str_split_fixed(as.character(wholeDF[,1]), "_", 3), stringsAsFactors = FALSE)
colnames(splitDF) <- c("chromosome", "what", "count")
splitDF$count <- as.numeric(as.character(splitDF$count))

splitDF%>% group_by(what) %>% summarise(total = sum(count)) %>% data.frame()-> sumDF

write.table(sumDF, file = outfile, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

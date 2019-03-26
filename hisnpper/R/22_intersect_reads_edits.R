suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
edit_desired <- args[i+1]
input_fromAwk_file <- args[i+2]
input_read_barcode_file <- args[i+3]
keep_positions_file  <- args[i+4]
remove_positions_file  <- args[i+5]
out_filtered_HQ_file <- args[i+6]
out_stats_file <- args[i+7]

if(FALSE){
  base <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/hisnpper/tests/edits"
  edit_desired <- "N_G"
  input_fromAwk_file <- paste0(base, "/","temp/02_frombam/rawFromBam_chr1.txt")
  input_read_barcode_file <- paste0(base, "/","temp/01_split/splitBam.chr1.read_barcode.tsv.gz")
  keep_positions_file <- "none"
  remove_positions_file <- "none"
  out_filtered_HQ_file <- paste0(base, "/","temp/03_processed/merged_chr1.txt")
  out_stats_file <- paste0(base, "/","logs/readstats/quant1.chr1.txt")
}


# Import chromosome for analysis
reads_from_awk <- fread(input_fromAwk_file)
colnames(reads_from_awk) <- c("read_id", "chr", "pos", "edit", "BQletter")

# Import and filter for specific files
if(keep_positions_file != "none" & file.exists(keep_positions_file)){
  if(file.size(keep_positions_file) > 0){
    filt_dt <- fread(keep_positions_file)
    boo <- reads_from_awk$pos %in% filt_dt[[2]]
    reads_from_awk <- reads_from_awk[boo]
  }
}

if(remove_positions_file != "none"& file.exists(remove_positions_file)){
  if(file.size(remove_positions_file) > 0){
    filt_dt <- fread(remove_positions_file)
    boo <- reads_from_awk$pos %ni% filt_dt[[2]]
    reads_from_awk <- reads_from_awk[boo]
  }
}

# Filter out for the edit
reads_from_awk <-  reads_from_awk[edit == edit_desired,]

# Numeric conversion of BQ
BQvec <- 0:40
names(BQvec) <-  c("!",'"',"#","$","%","&","'","(",")","*","+",",","-",".","/",
                   "0","1","2","3","4","5","6","7","8","9",":",";","<","=",">",
                   "?","@","A","B","C","D","E","F","G","H","I")
reads_from_awk[,BQ := BQvec[BQletter]]

# Merge read_id with barcode
read_barcode <- fread(cmd = paste0("zcat < ", input_read_barcode_file), col.names = c("read_id", "barcode_id"))

# Remove duplicates from read_barcode
read_barcode <- unique(read_barcode)
N_reads <- dim(read_barcode)[1]

# Remove reads with an NA barcode
read_barcode <- read_barcode[complete.cases(read_barcode)]
N_reads_bc <- dim(read_barcode)[1]

# Now merge
if(all(read_barcode$barcode_id == "none")) {
  mdf <- reads_from_awk[,c("chr", "pos", "edit", "BQ", "read_id")]
  
} else {
  mdf <- merge(reads_from_awk, read_barcode)
  mdf <- mdf[,c("chr", "pos", "edit", "BQ", "barcode_id", "read_id")]
}



# Write out polished, processed data
# TO DO: update BQ to BQ_score
write.table(mdf, file = out_filtered_HQ_file, 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Write out initial QC stats
mdf %>% group_by(read_id) %>%
  summarize(count_n = n()) %>%
  group_by(count_n) %>%
  summarize(n_reads = n()) -> count_df

n_0 <- N_reads - sum(count_df$n_reads)
qc_df <- rbind(data.frame(count_n = 0, n_reads = n_0), count_df)
write.table(qc_df, file = out_stats_file,
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

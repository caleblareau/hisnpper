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
input_snp_file <- args[i+1]
input_fromAwk_file <- args[i+2]
input_read_barcode_file <- args[i+3]
out_filtered_HQ_file <- args[i+4]
out_stats_file <- args[i+5]

if(FALSE){
  base <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/hisnpper/tests"
  input_snp_file <- paste0(base, "/","out/temp/01_split/SNPs_chr1.tsv")
  input_fromAwk_file <- paste0(base, "/","out/temp/02_frombam/rawFromBam_chr1.txt")
  input_read_barcode_file <- paste0(base, "/","out/temp/01_split/splitBam.chr1.read_barcode.tsv.gz")
  out_filtered_HQ_file <- paste0(base, "/","out/temp/03_processed/merged_chr1.txt")
  out_stats_file <- paste0(base, "/","out/logs/readstats/quant1.chr1.tsv")
  
}

# Import chromosome for analysis
reads_from_awk <- fread(input_fromAwk_file)
colnames(reads_from_awk) <- c("chr", "bp", "pos_in_read", "BQletter", "base", "read_id")
reads_from_awk[,pos := bp + pos_in_read -1]

# Import and filter for SNPs
snp_dt <- fread(input_snp_file)
boo <- reads_from_awk$pos %in% snp_dt[["V2"]]
boo <- TRUE # for testing UPDATE
filt_snps <- reads_from_awk[boo]

# Numeric conversion of BQ
BQvec <- 0:40
names(BQvec) <-  c("!",'"',"#","$","%","&","'","(",")","*","+",",","-",".","/",
                "0","1","2","3","4","5","6","7","8","9",":",";","<","=",">",
                "?","@","A","B","C","D","E","F","G","H","I")
filt_snps[,BQ := BQvec[BQletter]]

# Merge read_id with barcode
read_barcode <- fread(cmd = paste0("zcat < ", input_read_barcode_file), col.names = c("read_id", "barcode_id"))

# Remove duplicates from read_barcode
read_barcode <- unique(read_barcode)
N_reads <- dim(read_barcode)[1]

# Remove reads with an NA barcode
read_barcode <- read_barcode[complete.cases(read_barcode)]
N_reads_bc <- dim(read_barcode)[1]

# Now merge
mdf <- merge(filt_snps, read_barcode)
N_reads_bc_SNP <- length(unique(mdf[["read_id"]]))

# Write out polished, processed data
# TO DO: update BQ to BQ_score
write.table(mdf[,c("chr", "pos", "base", "BQ", "barcode_id", "read_id")], file = out_filtered_HQ_file, 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Write out initial QC stats
qc_df <- data.frame(
  what = c("N_reads", "N_reads_bc", "N_reads_bc_SNP"),
  howmany = c(N_reads, N_reads_bc, N_reads_bc_SNP)
)

write.table(qc_df, file = out_stats_file,
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

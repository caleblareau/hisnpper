suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(reshape2)))

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
input_snp_file <- args[i+1]
input_filtered_HQ_file <- args[i+2]
out_conf_file <- args[i+3]

if(FALSE){
  base <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/hisnpper/tests"
  input_snp_file <- paste0(base, "/","out/temp/01_split/SNPs_chr1.tsv")
  input_filtered_HQ_file <- paste0(base, "/","out/temp/03_processed/merged_chr1.txt")
  out_conf_file <- paste0(base, "/","out/temp/04_assign/assigned_chr1.txt")
}

# Import SNPs parsed
reads_dt <- fread(input_filtered_HQ_file, col.names = c("chr", "pos", "base", "BQ", "barcode", "read"))

# Annotate SNP read with haplotype
snp_dt <- fread(input_snp_file)
m_snp_dt <- data.table(reshape2::melt(data.frame(snp_dt), id.vars = c("chr", "start")))
colnames(m_snp_dt) <- c("chr", "pos", "haplotype", "base")
merge_dt <- merge(reads_dt, m_snp_dt, by = c("chr", "pos", "base"))
merge_dt <- merge_dt[complete.cases(merge_dt)]

# Sum up BQ for pertinent reads
#agg_dt <- merge_dt[, .(BQ_transform = (1 - 10^(-1*sum(BQ/10)))), by = list(read, haplotype)] 
agg_dt <- merge_dt[, .(BQ_transform = sum(BQ)), by = list(read, haplotype)] 
dcast_dt <- dcast(agg_dt, read~haplotype, value.var = "BQ_transform", fill = 0)

# Polish up and assign
assign <- colnames(dcast_dt[,2:3])[max.col(dcast_dt[,2:3],ties.method="random")]
bq_ratio <- pmax(dcast_dt[[2]], dcast_dt[[3]]) / (dcast_dt[[2]] + dcast_dt[[3]])

out_df <- data.frame(read_id = dcast_dt[[1]], assign, round(bq_ratio, 2))
write.table(out_df, file = out_conf_file,
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


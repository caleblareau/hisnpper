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
input_qc1 <- args[i+1]
input_conf <- args[i+2]
cutoff <- as.numeric(args[i+3])

if(FALSE){
  base <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/hisnpper/tests"
  input_qc1 <- paste0(base, "/","out/mESC_test.qcQuant.tsv")
  input_conf <- paste0(base, "/","out/mESC_test.haplo_conf.tsv")
  cutoff <- 0.9
  
}

# Import existing QC
qc_df <- read.table(input_qc1); colnames(qc_df) <- c("call", "count")

# Import overall confidence calls
conf_dt <- fread(input_conf)
conf_count_df <- conf_dt %>% mutate(call = ifelse(V3 > cutoff, V2, "ambiguous")) %>%
  group_by(call) %>% summarize(count = n())

# append percentages and confidence amounts
maxx <- max(qc_df$count)
odf <- rbind(qc_df, conf_count_df) %>% arrange(desc(count)) %>%
  mutate(prop = round((count / maxx)*100,2))
colnames(odf) <- c("Category", "count", "Percentage")

write.table(odf, file = input_qc1, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


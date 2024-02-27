library(ape)
library(ggplot2)

## plot number of peaks per sample and peak selection method ##

sample_id_vec <- c("Bak13","Bak8","De_Donno","Filmore","LA_Y3C","LM10","OC8","Oak_35874","Red_Oak2","Riv16","Riv19","Temecula_1")
peak_vec <- c(
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20",
  "peaks_uBH_0.001"
)

num_peaks_list <- list()

for(sample_id in sample_id_vec) {
  for(peak_id in peak_vec) {
    sample_name <- paste(sample_id, peak_id, sep = "_")
    input_peaks <- paste0("/data/samples/",sample_id,"/peaks/",sample_name,".fasta")
    
    input_peak_DNAbin <- read.dna(file = input_peaks, format = "fasta")
    num_peaks_list[[sample_name]] <- data.frame("sample" = sample_id, "peak_id" = peak_id, "num_peaks" = nrow(input_peak_DNAbin))
  }
}

num_peaks_df <- do.call(rbind, num_peaks_list)

input_file_map <- c(
  "peaks_uBH_0.001" = "naive peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20" = "advanced peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25" = "advanced peak calling with k-mer"
)
num_peaks_df$peak_calling_method <- input_file_map[num_peaks_df$peak_id]
num_peaks_df$peak_calling_method <- factor(num_peaks_df$peak_calling_method, levels = c("naive peak calling","advanced peak calling","advanced peak calling with k-mer"))

ggplot(num_peaks_df, aes(x = sample, y = num_peaks, fill = peak_calling_method)) +
  geom_bar(position = "dodge2", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"), name = "peak calling method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("number of peaks")
ggsave("/data/output/results/number_of_peaks_barplot.png", width = 10, height = 5)

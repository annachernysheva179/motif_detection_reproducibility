library(stringdist)
library(ggplot2)
library(pheatmap)
library(stringr)
library(ape)

## read function:
source("/data/motif_detection_reproducibility/helper_functions//lib.R")

## iterate through all samples:
sample_id_vec <- c("Bak13","Bak8","De_Donno","Filmore","LA_Y3C","LM10","OC8","Oak_35874","Red_Oak2","Riv16","Riv19","Temecula_1")
peak_vec <- c(
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20",
  "peaks_uBH_0.001"
)
motif_map <- list(
  "Bak13" = c("GACNNNNNGTG","GATGNNNNNNNTCC","ATCGAT","TTCGAA","GCCAGG"),
  "Bak8" = c("GACNNNNNGTG","GATGNNNNNNNTCC","TTCGAA","GCCAGG"),
  "De_Donno" = c("ACNNNNNNNNGT","CCANNNNNNTACC","GCCAGG"),
  "Filmore" = c("GACNNNNNNNGTC","GAGNNNNNTGA","CCSGG","CGTACG"),
  "LA_Y3C" = c("CCANNNNNNTGG","GAANNNNNNRTGT","CGTACG"),
  "LM10" = c("GACNNNNNNNGTC","GAGNNNNNTGA","CCSGG","CGTACG","GCCGGC"),
  "OC8" = c("AGACC"),
  "Oak_35874" = c("GAGNNNNNTGA"),
  "Red_Oak2" = c("GACNNNNNNNGTC","GAANNNNNNNCTC","CCSGG","CGTACG"),
  "Riv16" = c("GACNNNNNNGTC","TTCGAA"),
  "Riv19" = c("TTCGAA"),
  "Temecula_1" = c("GACNNNNNGTG","GATGNNNNNNNTCC","ATCGAT","TTCGAA","GCCAGG")
)

benchmark_plot_list <- list()
benchmark_res_list <- list()
for(sample_id in sample_id_vec) {
  for(peak_id in peak_vec) {
    sample_name <- paste(sample_id, peak_id, sep = "_")
    print(paste0("Running benchmark for sample: ",sample_name))
    
    benchmark_res <- benchmark_one_sample(diff_path = paste0("/data/samples/",sample_id,"/diff/",sample_id,"_difference_FN.RDS"),
                         input_peaks = paste0("/data/samples/",sample_id,"/peaks/",sample_name,".fasta"),
                         true_motifs = motif_map[[sample_id]],
                         cisfinder_path = paste0("/data/output/results/cisFinder_res/cisFinder_",sample_name,".csv"),
                         vcnn_path = paste0("/data/output/results/vCNN_res/vCNN_",sample_name,".csv"),
                         meme_chip_path = paste0("/data/output/results/memechip_res/memechip_",sample_name,".csv"),
                         comos_path = paste0("/data/output/results/CoMos/CoMoS_",sample_id,"_sel4.0_amb2.0.csv"),
                         meme_path = paste0("/data/output/results/meme_res/meme_",sample_name,".csv"),
                         plot_path = "/data/output/benchmark/",
                         sample_name = sample_name)
    benchmark_plot_df <- get_plot_df(benchmark_res, sample_id = sample_id)
    ggplot(benchmark_plot_df, aes(tool, value)) +
      geom_bar(stat = "identity") +
      facet_wrap(~metric, ncol = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("/data/output/benchmark/",sample_name,"_metric_barplot.png"), width = 7, height = 5)
    benchmark_plot_df$input_file <- peak_id
    benchmark_plot_list[[sample_name]] <- benchmark_plot_df
    benchmark_res_list[[sample_name]] <- benchmark_res
  }
}

saveRDS(benchmark_plot_list, "/data/output/benchmark/benchmark_plot_list_complete.rds")
saveRDS(benchmark_res_list, "/data/output/benchmark/benchmark_res_list_complete.rds")

################################################################################################################################


## combine all samples:
combined_plot_df <- do.call(rbind, benchmark_plot_list)
input_file_map <- c(
  "peaks_uBH_0.001" = "naive peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20" = "advanced peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25" = "advanced peak calling with k-mer"
)
combined_plot_df$peak_calling_method <- input_file_map[combined_plot_df$input_file]
combined_plot_df$peak_calling_method <- factor(combined_plot_df$peak_calling_method, levels = c("naive peak calling","advanced peak calling","advanced peak calling with k-mer"))

ggplot(combined_plot_df, aes(x = tool, y = value, fill = peak_calling_method)) +
  geom_boxplot(position = "dodge2") +
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"), name = "peak calling method") +
  facet_wrap(~metric, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("/data/output/benchmark/all_samples_tool_comparison_metric_scatter.png", width = 10, height = 5)


################################################################################################################################


## benchmark refined motifs ##
refined_bm_plot_list <- list()
refined_bm_res_list <- list()

for(sample_id in sample_id_vec) {
  for(peak_id in peak_vec) {
    sample_name <- paste(sample_id, peak_id, sep = "_")
    print(paste0("Running refined benchmark for sample: ",sample_name))
    
    benchmark_res <- try(benchmark_one_sample(diff_path = paste0("/data/samples/",sample_id,"/diff/",sample_id,"_difference_FN.RDS"),
                                          input_peaks = paste0("/data/samples/",sample_id,"/peaks/",sample_name,".fasta"),
                                          true_motifs = motif_map[[sample_id]],
                                          cisfinder_path = paste0("/data/output/benchmark/postprocess/cisFinder_res/",sample_id,"/",peak_id,"/Motifs_classification_cisFinder_res_",sample_id,"_nn_model.tsv"),
                                          vcnn_path = paste0("/data/output/benchmark/postprocess/vCNN_res/",sample_id,"/",peak_id,"/Motifs_classification_vCNN_res_",sample_id,"_nn_model.tsv"),
                                          meme_path = paste0("/data/output/benchmark/postprocess/meme_res/",sample_id,"/",peak_id,"/Motifs_classification_meme_res_",sample_id,"_nn_model.tsv"),
                                          comos_path = paste0("/data/output/benchmark/postprocess/CoMos/",sample_id,"/sel4.0_amb2.0/Motifs_classification_CoMos_",sample_id,"_nn_model.tsv"),
                                          plot_path = "/data/output/benchmark/postprocess/plots/",
                                          sample_name = sample_name,
                                          input_sep = "\t"))
    if(inherits(benchmark_res, "try-error")) {
      ## if error: skip
      next
    }
    benchmark_plot_df <- get_plot_df(benchmark_res, sample_id = sample_id)
    ggplot(benchmark_plot_df, aes(tool, value)) +
      geom_bar(stat = "identity") +
      facet_wrap(~metric, ncol = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("/data/output/benchmark/postprocess/plots/",sample_name,"_metric_barplot.png"), width = 7, height = 5)
    benchmark_plot_df$input_file <- peak_id
    refined_bm_plot_list[[sample_name]] <- benchmark_plot_df
    refined_bm_res_list[[sample_name]] <- benchmark_res
  }
}

saveRDS(refined_bm_plot_list, "/data/output/benchmark/postprocess/refined_benchmark_plot_list.rds")
saveRDS(refined_bm_res_list, "/data/output/benchmark/postprocess/refined_benchmark_result_list.rds")



########################################################################
#### compare benchmark before and after nanodisco motif refinement #####

benchmark_plot_list <- readRDS("/data/output/benchmark/benchmark_plot_list_complete.rds")
benchmark_plot_list <- benchmark_plot_list[names(benchmark_plot_list) %in% names(refined_bm_plot_list)]

combined_plot_df <- do.call(rbind, benchmark_plot_list)
refined_combined_plot_df <- do.call(rbind, refined_bm_plot_list)

diff_combined_plot_df <- combined_plot_df
diff_combined_plot_df$value_difference <- refined_combined_plot_df$value - combined_plot_df$value

input_file_map <- c(
  "peaks_uBH_0.001" = "naive peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20" = "advanced peak calling",
  "peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20_k_3_kmer_quantile_0.25" = "advanced peak calling with k-mer"
)
diff_combined_plot_df$peak_calling_method <- input_file_map[diff_combined_plot_df$input_file]
diff_combined_plot_df$peak_calling_method <- factor(diff_combined_plot_df$peak_calling_method, levels = c("naive peak calling","advanced peak calling","advanced peak calling with k-mer"))

ggplot(diff_combined_plot_df, aes(x = tool, y = value_difference, fill = peak_calling_method)) +
  geom_boxplot(position = "dodge2") +
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"), name = "peak calling method") +
  facet_wrap(~metric, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("/data/output/benchmark/postprocess/plots/all_samples_tool_comparison_difference_metric_scatter.png", width = 10, height = 5)

## also plot absolute value for refined motifs:
refined_combined_plot_df$peak_calling_method <- input_file_map[refined_combined_plot_df$input_file]
refined_combined_plot_df$peak_calling_method <- factor(refined_combined_plot_df$peak_calling_method, levels = c("naive peak calling","advanced peak calling","advanced peak calling with k-mer"))

ggplot(refined_combined_plot_df, aes(x = tool, y = value, fill = peak_calling_method)) +
  geom_boxplot(position = "dodge2") +
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"), name = "peak calling method") +
  facet_wrap(~metric, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("/data/output/benchmark/postprocess/plots/all_samples_tool_comparison_metric_scatter.png", width = 10, height = 5)



########################################################################
############################ runtime ###################################

## runtime in minutes (for simplicity)
runtime_df <- data.frame(
  "tool" = rep(c("meme-chip","vCNN","cisFinder","CoMoS","meme"), each = 3),
  "sample" = rep(c("Riv19","Temecula_1","LA_Y3C"), 5),
  "runtime" = c(12.45833,22.01403,15.22482,0.54755,4.37351,2.54337,0.03407,0.07318,0.05759,0.39758,0.39656,0.44798,42.00140292,42.56187098,42.51947737)
)

ggplot(runtime_df, aes(x = tool, y = runtime, color = sample)) +
  geom_jitter(width = 0.3, size = 2) +
  theme_bw() +
  ylab("runtime (mins)")
ggsave("/data/output/benchmark/runtime_scatter.png", width = 8, height = 4)



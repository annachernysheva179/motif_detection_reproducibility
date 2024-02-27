## pre-process diff files ##
library(ape)
library(GenomicRanges)
library(stringr)


source("/data/motif_detection_reproducibility/helper_functions//lib.R")

pre_process_diff <- function(diff_path, ref_path, out_path) {
  diff <- readRDS(diff_path)
  
  ## remove NA rows
  diff <- diff[!is.na(diff$mean_diff), ]
  
  ## read fasta file:
  ref_seq <- read.dna(file = ref_path, format = "fasta")
  ref_seq <- as.character.DNAbin(ref_seq)
  
  ## fix index:
  genomic_index <- rep(NA, nrow(diff))
  genomic_index[diff$dir == "fwd"] <- diff$position[diff$dir == "fwd"] + 3
  genomic_index[diff$dir == "rev"] <- diff$position[diff$dir == "rev"] + 4
  
  diff$genomic_index <- genomic_index
  
  print("Fetching genomic sequence..")
  ### also add sequence to position ###
  diff$genomic_sequence <- sapply(1:nrow(diff), function(idx) {
    if(diff$dir[idx] == "fwd") {
      toupper(paste(ref_seq[1, max(1,diff$genomic_index[idx]-22):min(diff$genomic_index[idx]+22, ncol(ref_seq))], collapse = ""))
    } else {
      toupper(paste(reverse_complement(ref_seq[1, max(1,diff$genomic_index[idx]-22):min(diff$genomic_index[idx]+22, ncol(ref_seq))]), collapse = ""))
    }
  })
  
  print("Doing multiple testing correction..")
  ## add BH corrected p-value
  diff$u_test_pval_BH <- p.adjust(diff$u_test_pval, method = "BH")
  
  print("Saving object..")
  saveRDS(diff, out_path)
  print("Done!")
}

## pre-process all samples:
pre_process_diff(diff_path = "/data/samples/Riv19/diff/Riv19_difference.RDS",
                 ref_path = "/data/samples/Riv19/ref_genome/Riv19.fasta",
                 out_path = "/data/samples/Riv19/diff/Riv19_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Filmore/diff/Filmore_difference.RDS",
                 ref_path = "/data/samples/Filmore/ref_genome/Filmore.fasta",
                 out_path = "/data/samples/Filmore/diff/Filmore_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/LA_Y3C//diff/LA_Y3C_difference.RDS",
                 ref_path = "/data/samples/LA_Y3C/ref_genome/LA_Y3C.fasta",
                 out_path = "/data/samples/LA_Y3C/diff/LA_Y3C_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Bak13/diff/Bak13_difference.RDS",
                 ref_path = "/data/samples/Bak13/ref_genome/Bak13.fasta",
                 out_path = "/data/samples/Bak13/diff/Bak13_difference_FN.RDS")


pre_process_diff(diff_path = "/data/samples/Bak8/diff/Bak8_difference.RDS",
                 ref_path = "/data/samples/Bak8/ref_genome/Bak8.fasta",
                 out_path = "/data/samples/Bak8/diff/Bak8_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/De_Donno/diff/De_Donno_difference.RDS",
                 ref_path = "/data/samples/De_Donno/ref_genome/De_Donno.fasta",
                 out_path = "/data/samples/De_Donno/diff/De_Donno_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/LM10/diff/LM10_difference.RDS",
                 ref_path = "/data/samples/LM10/ref_genome/LM10.fasta",
                 out_path = "/data/samples/LM10/diff/LM10_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/OC8/diff/OC8_difference.RDS",
                 ref_path = "/data/samples/OC8/ref_genome/OC8.fasta",
                 out_path = "/data/samples/OC8/diff/OC8_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Oak_35874/diff/Oak_35874_difference.RDS",
                 ref_path = "/data/samples/Oak_35874/ref_genome/Oak_35874.fasta",
                 out_path = "/data/samples/Oak_35874/diff/Oak_35874_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Red_Oak2/diff/Red_Oak2_difference.RDS",
                 ref_path = "/data/samples/Red_Oak2/ref_genome/Red_Oak2.fasta",
                 out_path = "/data/samples/Red_Oak2/diff/Red_Oak2_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Riv16/diff/Riv16_difference.RDS",
                 ref_path = "/data/samples/Riv16/ref_genome/Riv16.fasta",
                 out_path = "/data/samples/Riv16/diff/Riv16_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Riv25/diff/Riv25_difference.RDS",
                 ref_path = "/data/samples/Riv25/ref_genome/Riv25.fasta",
                 out_path = "/data/samples/Riv25/diff/Riv25_difference_FN.RDS")

pre_process_diff(diff_path = "/data/samples/Temecula_1/diff/Temecula_1_difference.RDS",
                 ref_path = "/data/samples/Temecula_1/ref_genome/Temecula_1.fasta",
                 out_path = "/data/samples/Temecula_1/diff/Temecula_1_difference_FN.RDS")




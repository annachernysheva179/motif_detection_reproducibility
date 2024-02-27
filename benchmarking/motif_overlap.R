## motif overlap ##
library(ggvenn)


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

found_motif_list <- list()
for(sample_id in sample_id_vec) {
  for(peak_id in peak_vec) {
    sample_name <- paste(sample_id, peak_id, sep = "_")
    
    cisfinder_path = paste0("/data/output/results/cisFinder_res/cisFinder_",sample_name,".csv")
    vcnn_path = paste0("/data/output/results/vCNN_res/vCNN_",sample_name,".csv")
    meme_chip_path = paste0("/data/output/results/memechip_res/memechip_",sample_name,".csv")
    comos_path = paste0("/data/output/results/CoMos/CoMoS_",sample_id,"_sel4.0_amb2.0.csv")
    meme_path = paste0("/data/output/results/meme_res/meme_",sample_name,".csv")
    
    input_sep <- ","
    cisfinder_res <- read.table(cisfinder_path, h=T, sep = input_sep)
    vCNN_res <- read.table(vcnn_path, h=T, sep = input_sep, row.names = 1)
    meme_chip_res <- read.table(meme_chip_path, h=T, sep = input_sep, row.names = 1)
    comos_res <- read.table(comos_path, h=T, sep = input_sep, row.names = 1)
    meme_res <- read.table(meme_path, h=T, sep = input_sep, row.names = 1)
    
    found_motif_list[[sample_id]][[peak_id]] <- list(
      "cisfinder"=cisfinder_res$motif, "vcnn" = vCNN_res$motif,
      "meme_chip"=meme_chip_res$motif, "comos"=comos_res$motif, "meme"=meme_res$motif
    )
  }
}

## how often are comos motifs supported
for(sample_id in sample_id_vec) {
  
  unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["comos"]]})))
  
  print(paste0("Sample: ", sample_id))
  comos_motifs <- unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["comos"]]})))
  print(paste0("Motifs found by CoMoS: ", paste(comos_motifs, collapse = ", ")))
  
  cisfinder_motifs <- unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["cisfinder"]]})))
  print(paste0("Motifs supported by cisfinder: ", paste(unique(cisfinder_motifs[cisfinder_motifs %in% comos_motifs]), collapse = ", ")))
  
  vcnn_motifs <- unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["vcnn"]]})))
  print(paste0("Motifs supported by vCNN: ", paste(unique(vcnn_motifs[vcnn_motifs %in% comos_motifs]), collapse = ", ")))
  
  meme_chip_motifs <- unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["meme_chip"]]})))
  print(paste0("Motifs supported by MEME-ChIP: ", paste(unique(meme_chip_motifs[meme_chip_motifs %in% comos_motifs]), collapse = ", ")))
  
  meme_motifs <- unique(unlist(lapply(found_motif_list[[sample_id]], function(el) {el[["meme"]]})))
  print(paste0("Motifs supported by MEME: ", paste(unique(meme_motifs[meme_motifs %in% comos_motifs]), collapse = ", ")))
  
  print("")
  print("")
}

### try for venn diagrams ##
for(sample_id in sample_id_vec) {
  found_list_sub <- found_motif_list[[sample_id]][["peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20"]]
  comp_list <- list(
    "CoMoS" = found_list_sub$comos,
    "MEME ChIP" = found_list_sub$meme_chip,
    "vCNN" = found_list_sub$vcnn,
    "cisFinder" = found_list_sub$cisfinder
  )
  ggvenn(comp_list, show_elements = F, show_percentage = T, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
         stroke_size = 0.5, set_name_size = 4) +
    ggtitle(paste0("Venn diagram for ", sample_id))
  ggsave(paste0("/data/output/benchmark/motif_overlap/",sample_id, "_motif_overlap_venn.png"), width = 6, height = 5)
}


##########################################################################################
####################### same for refined motifs ##########################################

found_refined_motif_list <- list()
for(sample_id in sample_id_vec) {
  for(peak_id in peak_vec) {
    sample_name <- paste(sample_id, peak_id, sep = "_")
    
    cisfinder_path = paste0("/data/output/benchmark/postprocess/cisFinder_res/",sample_id,"/",peak_id,"/Motifs_classification_cisFinder_res_",sample_id,"_nn_model.tsv")
    vcnn_path = paste0("/data/output/benchmark/postprocess/vCNN_res/",sample_id,"/",peak_id,"/Motifs_classification_vCNN_res_",sample_id,"_nn_model.tsv")
    meme_path = paste0("/data/output/benchmark/postprocess/meme_res/",sample_id,"/",peak_id,"/Motifs_classification_meme_res_",sample_id,"_nn_model.tsv")
    comos_path = paste0("/data/output/benchmark/postprocess/CoMos/",sample_id,"/sel4.0_amb2.0/Motifs_classification_CoMos_",sample_id,"_nn_model.tsv")
    
    input_sep <- "\t"
    if(file.exists(cisfinder_path)) {
      cisfinder_res <- read.table(cisfinder_path, h=T, sep = input_sep)
      if(nrow(cisfinder_res) == 0) {cisfinder_res <- NULL}
    } else {
      print(paste0("Invalid cisfinder path: ", cisfinder_path))
      cisfinder_res <- NULL
    }
    if(file.exists(vcnn_path)) {
      vCNN_res <- read.table(vcnn_path, h=T, sep = input_sep, row.names = 1)
      if(nrow(vCNN_res) == 0) {vCNN_res <- NULL}
    } else {
      print(paste0("Invalid vCNN path: ", vcnn_path))
      vCNN_res <- NULL
    }
    if(file.exists(meme_path)) {
      meme_chip_res <- read.table(meme_path, h=T, sep = input_sep, row.names = 1)
      if(nrow(meme_chip_res) == 0) {meme_chip_res <- NULL}
    } else {
      print(paste0("Invalid meme-chip path: ", meme_path))
      meme_chip_res <- NULL
    }
    if(file.exists(comos_path)) {
      comos_res <- read.table(comos_path, h=T, sep = input_sep, row.names = 1)
      if(nrow(comos_res) == 0) {comos_res <- NULL}
    } else {
      print(paste0("Invalid CoMoS path: ", comos_path))
      comos_res <- NULL
    }
    
    found_refined_motif_list[[sample_id]][[peak_id]] <- list(
      "cisfinder"=cisfinder_res$Motif, "vcnn" = vCNN_res$Motif,
      "meme_chip"=meme_chip_res$Motif, "comos"=comos_res$Motif
    )
  }
}

for(sample_id in sample_id_vec) {
  found__refined_list_sub <- found_refined_motif_list[[sample_id]][["peaks_uBH_0.001_peak_dist_2_min_cov_20_min_dist_20"]]
  comp_list <- list(
    "CoMoS" = found__refined_list_sub$comos,
    "MEME ChIP" = found__refined_list_sub$meme_chip,
    "vCNN" = found__refined_list_sub$vcnn,
    "cisFinder" = found__refined_list_sub$cisfinder
  )
  ggvenn(comp_list, show_elements = F, show_percentage = T, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
         stroke_size = 0.5, set_name_size = 4) +
    ggtitle(paste0("Venn diagram using refined motifs for: ", sample_id))
  ggsave(paste0("/data/output/benchmark/motif_overlap/",sample_id, "_refined_motif_overlap_venn.png"), width = 6, height = 5)
}


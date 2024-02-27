## helper functions ##


## BENCHMARKING ##
convert_motif <- function(motif) {
  map <- c(
    "R"="[GA]",
    "Y"="[CT]",
    "M"="[AC]",
    "K"="[GT]",
    "S"="[GC]",
    "W"="[AT]",
    "B"="[CGT]",
    "D"="[AGT]",
    "H"="[ACT]",
    "V"="[ACG]",
    "N"="[ACGT]",
    "A"="A",
    "C"="C",
    "G"="G",
    "T"="T"
  )
  motif_split <- strsplit(motif, "")[[1]]
  motif_long <- map[motif_split]
  motif_long <- paste(motif_long, collapse = "")
  return(motif_long)
}

get_reverse_complement <- function(motif_vec) {
  compl_map <- c(
    "A"="T","C"="G","G"="C","T"="A",
    "R"="Y","Y"="R","M"="K","K"="M","S"="S","W"="W",
    "B"="V","D"="H","H"="D","V"="B","N"="N"
  )
  rc_vec <- sapply(motif_vec, function(el) {
    el_split <- strsplit(el,"")[[1]]
    el_split_rev <- el_split[length(el_split):1]
    el_split_rev_comp <- compl_map[el_split_rev]
    paste(el_split_rev_comp, collapse = "")
  })
  return(rc_vec)
}

get_motif_similarity <- function(res, true_motifs, out_file, sim_method = "osa") {
  ##account for null res:
  if(is.null(res)) {
    return(NA)
  }
  
  sim_res <- stringsimmatrix(true_motifs, res$motif, method = sim_method)
  rownames(sim_res) <- true_motifs; colnames(sim_res) <- res$motif
  sim_res <- as.data.frame(sim_res)
  
  ## check if rownames (tool motifs) contain same motif twice:
  already_visited <- c()
  for(i in 1:ncol(sim_res)) {
    if(i %in% already_visited) {next}
    current_motif <- colnames(sim_res)[i]
    other_motifs <- colnames(sim_res)[colnames(sim_res) != current_motif]
    if(get_reverse_complement(current_motif) %in% other_motifs) {
      rc_idx <-which(colnames(sim_res) == get_reverse_complement(current_motif))
      sim_res[,i] <- apply(sim_res[,c(i,rc_idx)], 1, max)
      colnames(sim_res)[i] <- paste(c(colnames(sim_res)[i], colnames(sim_res)[rc_idx]), collapse = "_")
      already_visited <- c(already_visited, rc_idx)
    }
  }
  all_idx <- 1:ncol(sim_res)
  if(ncol(sim_res) > 1) {
    sim_res <- sim_res[,all_idx[!all_idx %in% already_visited]]
  }
  
  
  if(ncol(sim_res) < 2 & nrow(sim_res) >= 2) {
    pheatmap(sim_res, cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = T, 
             color = colorRampPalette(c("grey", "blue"))(50), 
             filename = out_file, width = 12, height = 8)
  } else if(nrow(sim_res) < 2 & ncol(sim_res) >= 2) {
    pheatmap(sim_res, cluster_rows = F, cluster_cols = T, show_rownames = T, show_colnames = T, 
             color = colorRampPalette(c("grey", "blue"))(50), 
             filename = out_file, width = 12, height = 8)
  } else if(nrow(sim_res) < 2 & ncol(sim_res) < 2) {
    print(paste0("No heatmap plot for ", out_file))
  }else {
    pheatmap(sim_res, cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = T, 
             color = colorRampPalette(c("grey", "blue"))(50), 
             filename = out_file, width = 12, height = 8)
  }
  
  ## aggregate motif similarities
  motif_similarity_index <- mean(apply(sim_res, 2, max))
  return(motif_similarity_index)
}

get_confusion_matrix <- function(diff, res, true_motifs) {
  if(is.null(res)) {
    return(NA)
  }
  sites <- sapply(res$motif, function(el) {
    el <- convert_motif(el)
    str_detect(diff$genomic_sequence, el)
  })
  sites <- apply(sites, 1, sum)
  true_sites <- sapply(true_motifs, function(el) {
    el <- convert_motif(el)
    str_detect(diff$genomic_sequence, el)
  })
  true_sites <- apply(true_sites, 1, sum)
  tp <- sum(true_sites >= 1 & sites >= 1)
  fp <- sum(true_sites < 1 & sites >= 1)
  tn <- sum(true_sites < 1 & sites < 1)
  fn <- sum(true_sites >= 1 & sites < 1)
  
  m <- matrix(c(tp, fp, fn, tn), byrow = F, nrow = 2, 
              dimnames = list(c("site methylated","site unmethylated"), c("positive prediction","negative prediction")))
  return(m)
}

get_tpr <- function(m) {
  if(!is.na(m)) {
    return((m[1,1] / (m[1,1] + m[1,2])))
  } else {
    return(NA)
  }
}

get_tnr <- function(m) {
  if(!is.na(m)) {
    return((m[2,2] / (m[2,2] + m[2,1])))
  } else {
    return(NA)
  }
}

get_acc <- function(m) {
  if(!is.na(m)) {
    return(((m[1,1] + m[2,2]) / sum(m)))
  } else {
    return(NA)
  }
}

get_plot_df <- function(benchmark_res, sample_id) {
  plot_df <- data.frame(
    "metric" = rep(c("motif_similarity","TPR","TNR","ACC"), 5),
    "value" = c(benchmark_res$cisfinder_motif_similarity, benchmark_res$cisfinder_TPR, benchmark_res$cisfinder_TNR, benchmark_res$cisfinder_accuracy,
                benchmark_res$vcnn_motif_similarity, benchmark_res$vcnn_TPR, benchmark_res$vcnn_TNR, benchmark_res$vcnn_accuracy,
                benchmark_res$meme_chip_motif_similarity, benchmark_res$meme_chip_TPR, benchmark_res$meme_chip_TNR, benchmark_res$meme_chip_accuracy,
                benchmark_res$comos_motif_similarity, benchmark_res$comos_TPR, benchmark_res$comos_TNR, benchmark_res$comos_accuracy,
                benchmark_res$meme_motif_similarity, benchmark_res$meme_TPR, benchmark_res$meme_TNR, benchmark_res$meme_accuracy),
    "tool" = rep(c("cisFinder","vCNN","meme_chip", "CoMoS", "meme"), each = 4),
    "sample" = sample_id
  )
  return(plot_df)
}

fraction_explained_sites <- function(diff_obj, motif_vec, detailed_plots = FALSE, plot_path, sample_name, tool_name) {
  sites <- sapply(motif_vec, function(el) {
    el <- convert_motif(el)
    str_detect(diff_obj$genomic_sequence, el)
  })
  diff_value_list <- lapply(1:ncol(sites), function(i) {diff_obj$mean_diff[sites[,i]]})
  
  var_complete <- var(diff_obj$mean_diff)
  varex <- sapply(diff_value_list, var) / var_complete
  names(varex) <- colnames(sites)
  
  frac_sites <- apply(sites, 2, sum) / nrow(sites)
  
  if(detailed_plots) {
    ## fill up diff_value_list with NAs:
    max_el <- max(sapply(diff_value_list, length))
    for(i in 1:length(diff_value_list)) {
      if(!length(diff_value_list[[i]]) == max_el) {
        diff_value_list[[i]] <- c(diff_value_list[[i]], rep(NA, max_el - length(diff_value_list[[i]])))
      }
    }
    diff_value_df <- data.frame(
      "motif" = rep(colnames(sites), each = max_el),
      "diff_value" = unlist(diff_value_list)
    )
    ggplot(diff_value_df, aes(motif, diff_value)) +
      geom_violin() +
      geom_jitter() +
      theme_bw()+
      ggtitle(paste0("mean diff of sites with motifs (",tool_name,")\n in ", sample_name)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0(plot_path, tool_name, "_", sample_name, "_diff_at_motif_sites.png"), width = 8, height = 6)
  }
  
  ## aggregate bound sites:
  sites <- apply(sites, 1, sum)
  
  res <- list("sites_explained" = sum(sites >= 1) / length(sites), "sites_explained_by_motif" = frac_sites, "varience_explained_by_motif" = varex)
  return(res)
}

fix_colnames <- function(res) {
  if("Motif" %in% colnames(res)) {
    colnames(res)[colnames(res) == "Motif"] <- "motif"
    return(res)
  }
  return(res)
}

benchmark_one_sample <- function(diff_path, input_peaks, true_motifs, cisfinder_path, vcnn_path, meme_chip_path, comos_path, meme_path, plot_path, sample_name, input_sep = ",") {
  print("Reading Diff-object..")
  diff_obj <- readRDS(diff_path)
  ## todo: update this accordingly:
  ##diff_obj <- diff_obj[sample(1:nrow(diff_obj), size = 50000), ]
  ##diff_obj <- diff_obj[diff_obj$u_test_pval_BH < 0.01, ]
  #diff_obj <- diff_obj[order(diff_obj$u_test_pval_BH, decreasing = F), ]
  #diff_obj <- diff_obj[1:50000, ]
  
  input_peak_DNAbin <- read.dna(file = input_peaks, format = "fasta")
  peak_labels <- sapply(labels(input_peak_DNAbin), function(el) {
    el_split <- strsplit(el, "_")[[1]]
    paste(el_split[(length(el_split)-1):length(el_split)], collapse = "_")
  })
  diff_obj$position2 <- paste(diff_obj$dir, diff_obj$genomic_index, sep = "_")
  diff_obj <- diff_obj[diff_obj$position2 %in% peak_labels, ]
  ## todo: why not recover all sites?
  
  ## read results from tools:
  print("Reading tool input files..")
  if(file.exists(cisfinder_path)) {
    cisfinder_res <- read.table(cisfinder_path, h=T, sep = input_sep)
    cisfinder_res <- fix_colnames(cisfinder_res)
    if(nrow(cisfinder_res) == 0) {cisfinder_res <- NULL}
  } else {
    print(paste0("Invalid cisfinder path: ", cisfinder_path))
    cisfinder_res <- NULL
  }
  if(file.exists(vcnn_path)) {
    vCNN_res <- read.table(vcnn_path, h=T, sep = input_sep, row.names = 1)
    vCNN_res <- fix_colnames(vCNN_res)
    if(nrow(vCNN_res) == 0) {vCNN_res <- NULL}
  } else {
    print(paste0("Invalid vCNN path: ", vcnn_path))
    vCNN_res <- NULL
  }
  if(file.exists(meme_chip_path)) {
    meme_chip_res <- read.table(meme_chip_path, h=T, sep = input_sep, row.names = 1)
    meme_chip_res <- fix_colnames(meme_chip_res)
    if(nrow(meme_chip_res) == 0) {meme_chip_res <- NULL}
  } else {
    print(paste0("Invalid meme-chip path: ", meme_chip_res))
    meme_chip_res <- NULL
  }
  if(file.exists(comos_path)) {
    comos_res <- read.table(comos_path, h=T, sep = input_sep, row.names = 1)
    comos_res <- fix_colnames(comos_res)
    if(nrow(comos_res) == 0) {comos_res <- NULL}
  } else {
    print(paste0("Invalid CoMoS path: ", comos_path))
    comos_res <- NULL
  }
  if(file.exists(meme_path)) {
    meme_res <- read.table(meme_path, h=T, sep = input_sep, row.names = 1)
    meme_res <- fix_colnames(meme_res)
    if(nrow(meme_res) == 0) {meme_res <- NULL}
  } else {
    print(paste0("Invalid meme path: ", meme_path))
    meme_res <- NULL
  }
  
  ## 1.) motif similarity to TP motifs
  print("Inferring motif similarities..")
  cisfinder_motif_similarity <- get_motif_similarity(cisfinder_res, true_motifs, out_file = paste0(plot_path,"motif_similarity/", sample_name, "_cisFinder_motif_similarity_heat.png"))
  vcnn_motif_similarity <- get_motif_similarity(vCNN_res, true_motifs, out_file = paste0(plot_path,"motif_similarity/", sample_name, "_vCNN_motif_similarity_heat.png"))
  meme_chip_motif_similarity <- get_motif_similarity(meme_chip_res, true_motifs, out_file = paste0(plot_path,"motif_similarity/", sample_name, "_meme_chip_motif_similarity_heat.png"))
  comos_motif_similarity <- get_motif_similarity(comos_res, true_motifs, out_file = paste0(plot_path,"motif_similarity/", sample_name, "_CoMoS_motif_similarity_heat.png"))
  meme_motif_similarity <- get_motif_similarity(meme_res, true_motifs, out_file = paste0(plot_path,"motif_similarity/", sample_name, "_meme_motif_similarity_heat.png"))
  
  
  ## 2.) modified site recovery (TP, FP)
  print("Calculating confusion matrices..")
  cisfinder_confusion_mtx <- get_confusion_matrix(diff_obj, cisfinder_res, true_motifs)
  vcnn_confusion_mtx <- get_confusion_matrix(diff_obj, vCNN_res, true_motifs)
  meme_chip_confusion_mtx <- get_confusion_matrix(diff_obj, meme_chip_res, true_motifs)
  comos_confusion_mtx <- get_confusion_matrix(diff_obj, comos_res, true_motifs)
  meme_confusion_mtx <- get_confusion_matrix(diff_obj, meme_res, true_motifs)
  
  ## 3.) how many sites are explained by each tool?
  
  if(file.exists(cisfinder_path) & file.exists(vcnn_path) & file.exists(meme_chip_path) & file.exists(comos_path) & file.exists(meme_path)) {
    frac_sites_explained <- list(
      "true_motifs" = fraction_explained_sites(diff_obj, motif_vec = true_motifs, detailed_plots = T,
                                               plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "true_motif"),
      "cisfinder" = fraction_explained_sites(diff_obj, motif_vec = cisfinder_res$motif, detailed_plots = T,
                                             plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "cisfinder"),
      "vcnn" = fraction_explained_sites(diff_obj, motif_vec = vCNN_res$motif, detailed_plots = T,
                                        plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "vcnn"),
      "meme_chip" = fraction_explained_sites(diff_obj, motif_vec = meme_chip_res$motif, detailed_plots = T,
                                             plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "meme_chip"),
      "comos" = fraction_explained_sites(diff_obj, motif_vec = comos_res$motif, detailed_plots = T,
                                         plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "comos"),
      "meme" = fraction_explained_sites(diff_obj, motif_vec = meme_res$motif, detailed_plots = T,
                                        plot_path = "/data/output/benchmark/fraction_explained/", sample_name = sample_name, tool_name = "meme")
    )
  } else {
    frac_sites_explained <- NA
  }
  
  
  res <- list(
    "cisfinder_motif_similarity" = cisfinder_motif_similarity,
    "cisfinder_confusion_mtx" = cisfinder_confusion_mtx,
    "cisfinder_TPR" = get_tpr(cisfinder_confusion_mtx), "cisfinder_TNR" = get_tnr(cisfinder_confusion_mtx), "cisfinder_accuracy" = get_acc(cisfinder_confusion_mtx),
    "vcnn_motif_similarity" = vcnn_motif_similarity,
    "vcnn_confusion_mtx" = vcnn_confusion_mtx,
    "vcnn_TPR" = get_tpr(vcnn_confusion_mtx), "vcnn_TNR" = get_tnr(vcnn_confusion_mtx), "vcnn_accuracy" = get_acc(vcnn_confusion_mtx),
    "meme_chip_motif_similarity" = meme_chip_motif_similarity,
    "meme_chip_confusion_mtx" = meme_chip_confusion_mtx,
    "meme_chip_TPR" = get_tpr(meme_chip_confusion_mtx), "meme_chip_TNR" = get_tnr(meme_chip_confusion_mtx), "meme_chip_accuracy" = get_acc(meme_chip_confusion_mtx),
    "comos_motif_similarity" = comos_motif_similarity,
    "comos_confusion_mtx" = comos_confusion_mtx,
    "comos_TPR" = get_tpr(comos_confusion_mtx), "comos_TNR" = get_tnr(comos_confusion_mtx), "comos_accuracy" = get_acc(comos_confusion_mtx),
    "meme_motif_similarity" = meme_motif_similarity,
    "meme_confusion_mtx" = meme_confusion_mtx,
    "meme_TPR" = get_tpr(meme_confusion_mtx), "meme_TNR" = get_tnr(meme_confusion_mtx), "meme_accuracy" = get_acc(meme_confusion_mtx),
    "fraction_sites_explained" = frac_sites_explained
  )
  if(!is.na(frac_sites_explained)) {
    res[["sites_explained_by_tool"]] <- data.frame("cisfinder"=frac_sites_explained$cisfinder$sites_explained,"vcnn"=frac_sites_explained$vcnn$sites_explained,
                                                   "meme_chip"=frac_sites_explained$meme_chip$sites_explained,"comos"=frac_sites_explained$comos$sites_explained,
                                                   "meme"=frac_sites_explained$meme$sites_explained,
                                                   "true" = frac_sites_explained$true_motifs$sites_explained)
  }
  print("Done!")
  print("-----------")
  return(res)
}







######################## NOT USED ANYMORE ###################################################
get_site_distances <- function(diff_obj) {
  diff_obj_fwd <- diff_obj[diff_obj$dir == "fwd", ]
  diff_obj_rev <- diff_obj[diff_obj$dir == "rev", ]
  
  dist_vec_fwd <- diff_obj_fwd$position[2:nrow(diff_obj_fwd)] - diff_obj_fwd$position[1:(nrow(diff_obj_fwd)-1)]
  dist_vec_rev <- diff_obj_rev$position[2:nrow(diff_obj_rev)] - diff_obj_rev$position[1:(nrow(diff_obj_rev)-1)]
  return(list("fwd"=dist_vec_fwd, "rev"=dist_vec_rev))
}
peak_size_runner <- function(diff_obj_sub) {
  res <- c()
  peak_size = 0
  for(i in 2:nrow(diff_obj_sub)) {
    if(diff_obj_sub$position[i] == diff_obj_sub$position[i-1]+1) {
      peak_size <- peak_size + 1
    } else {
      res <- c(res, peak_size)
      peak_size <- 0
    }
  }
  return(res)
}
get_peak_size <- function(diff_obj) {
  ## subset for forward and reverse strand
  diff_obj_fwd <- diff_obj[diff_obj$dir == "fwd", ]
  diff_obj_rev <- diff_obj[diff_obj$dir == "rev", ]
  ## infer peak sizes:
  peak_size_fwd <- peak_size_runner(diff_obj_fwd)
  peak_size_rev <- peak_size_runner(diff_obj_rev)
  ## return peak sizes:
  return(list("fwd"=peak_size_fwd, "rev"=peak_size_rev))
}
reverse_complement <- function(seq) {
  base_map <- c("a"="t", "c"="g","g"="c","t"="a")
  seq <- seq[length(seq):1]
  seq <- base_map[seq]
  return(seq)
}
## EXPORT PEAKS TO FASTA ##
export_fasta <- function(riv_diff_sub, filename, th, riv_ref_seq_mtx) {
  ## open connection:
  f <- file(filename)
  ## iteratively fill fasta file:
  fasta_content <- c()
  for(i in 1:nrow(riv_diff_sub)) {
    ## 1: header:
    h <- paste0(">riv19_diff_peaks_uBH_", th, "_contig_", riv_diff_sub$contig[i], "_", i, "_position_", riv_diff_sub$genomic_index[i], "_direction_", riv_diff_sub$dir[i])
    ## 2: sequence:
    seq_mtx <- paste(riv_ref_seq_mtx[1, max(1,riv_diff_sub$genomic_index[i]-22):min(riv_diff_sub$genomic_index[i]+22, ncol(riv_ref_seq_mtx))], collapse = "")
    ## append to content:
    fasta_content <- c(fasta_content, h)
    #fasta_content <- c(fasta_content, "\n")
    fasta_content <- c(fasta_content, seq_mtx)
    #fasta_content <- c(fasta_content, "\n")
  }
  ## write content
  writeLines(fasta_content, f)
  ## close connection
  close(f)
}
## COUNT KNOWN MOTIFS ##
count_motif_occurence <- function(riv_diff_sub, motifs, riv_ref_seq_mtx) {
  motif_count <- lapply(motifs, function(motif) {
    motif_occur <- str_detect(riv_diff_sub$genomic_sequence, motif)
    sum(motif_occur)
  })
  names(motif_count) <- motifs
  return(motif_count)
}
######################################################################################################


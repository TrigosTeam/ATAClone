fix_peak_range <- function(peak_info, target_width){
  l_width <- peak_info$summit - peak_info$start
  r_width <- peak_info$end - peak_info$summit
  new_start <- round(peak_info$summit - target_width * l_width / (l_width + r_width))
  new_end <- round(peak_info$summit + target_width * r_width / (l_width + r_width))
  peak_info$start <- new_start
  peak_info$end <- new_end
  peak_info
}


get_peaks_altius <- function(altius_peak_mat_path, altius_sample_info_path, altius_peak_info_path, n_systems, fixed_peak_width = NULL){
  require("Matrix")
  require("MatrixGenerics")
  require("GenomicRanges")

  #files from here: https://zenodo.org/records/3838751
  altius_peak_mat <- Matrix::readMM(altius_peak_mat_path)
  #for some reason, line 734 is malformed? It shouldn't be there anyway
  altius_sample_info <- read.delim(altius_sample_info_path, sep = "\t", header = T)[1:733,]
  altius_peak_info <- read.delim(altius_peak_info_path, sep = "\t", header = T)

  mat_list <- list()
  for(organ.system in sort(unique(altius_sample_info$System))){
    mat_list[[organ.system]] <- altius_peak_mat[,altius_sample_info$System == organ.system]
  }
  system_mat <- do.call(cbind, lapply(mat_list, rowSums))
  use_idx <- which(rowSums(system_mat != 0) == n_systems)

  altius_peak_info <- altius_peak_info[use_idx,]
  if (!is.null(fixed_peak_width) & is.integer(fixed_peak_width)){
    altius_peak_info <- fix_peak_range(altius_peak_info, fixed_peak_width)
  }
  GenomicRanges::makeGRangesFromDataFrame(altius_peak_info)
}

make_bins <- function(chr_arm_info, bin_width){
  chr_arm_width <- (chr_arm_info[2] - chr_arm_info[1])
  n_bins <- chr_arm_width / bin_width
  if (n_bins < 1){
    n_bins <- 1
  } else {
    n_bins <- ifelse(abs(chr_arm_width / floor(n_bins) - bin_width) < abs(chr_arm_width / ceiling(n_bins) - bin_width), floor(n_bins), ceiling(n_bins))
  }
  new_bin_width <- chr_arm_width / n_bins
  as.integer(floor(seq(chr_arm_info[1], chr_arm_info[2], new_bin_width)))
}

get_chr_arm_bins <- function(bin_width){
  #https://www.biostars.org/p/383786/
  chr.cyto.dt <- data.table::fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
                                   col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
  chr.arm.dt <- chr.cyto.dt[ , .(length = sum(chromEnd - chromStart)),
                             by = .(chrom, arm = substring(name, 1, 1)) ]
  chr.arm.dt <- chr.arm.dt[chrom %in% paste0("chr", c(1:22, "X", "Y"))]
  #reorder. Causes a warning - fix later
  #chr.arm.dt <- chr.arm.dt[order(as.integer(gsub("chr", "", chr.arm.dt$chrom))),]
  chr.arm.split <- split(chr.arm.dt, chr.arm.dt$chrom)
  chr.arm.coords <- list()
  for (chr in names(chr.arm.split)){
    chr.arm.coords[[chr]] <- list(p = as.integer(c(1, chr.arm.split[[chr]]$length[1])), q = as.integer(c(chr.arm.split[[chr]]$length[1] + 1, chr.arm.split[[chr]]$length[1] + chr.arm.split[[chr]]$length[2])))
  }
  chr.arm.coords <- unlist(chr.arm.coords, recursive = F)
  lapply(chr.arm.coords, make_bins, bin_width)
}

get_peak_overlaps <- function(atac.granges, peak.granges, method){
  if (method == "PIC"){
    frag.starts <- GenomicRanges::resize(atac.granges, 1, fix = "start")
    frag.ends <- GenomicRanges::resize(atac.granges, 1, fix = "end")
    incl.idx2.1 <- GenomicRanges::findOverlaps(frag.starts, peak.granges)
    incl.idx2.2 <- GenomicRanges::findOverlaps(frag.ends, peak.granges)
    atac.peak.granges <- atac.granges[unique(c(incl.idx2.1@from, incl.idx2.2@from))]
  } else {
    peak_overlaps <- GenomicRanges::findOverlaps(atac_fragments_filtered, atac_peaks)
    atac.peak.granges <- atac.granges[unique(peak_overlaps@from)]
  }
  atac.peak.granges
}

get_counts <- function(atac_fragments, bins){
  nbins <- length(bins)
  bin_tabs <- list()
  for (i in 2:nbins){
    bin_tabs[[paste0(bins[i - 1], "-", bins[i])]] <- table(atac_fragments[start >= bins[i - 1] & start < bins[i]]$name)
  }
  bin_tabs
}

get_bin_counts <- function(atac_fragments, chr_arm_bins, barcodes){
  atac_fragments_dt <- data.table::as.data.table(as.data.frame(atac_fragments))
  atac_chr_list <- split(atac_fragments_dt, atac_fragments_dt$seqnames)
  #is this the best place to put this?
  atac_chr_list <- atac_chr_list[names(atac_chr_list) %in% paste0("chr", c(1:22, "X", "Y"))]
  chr_names <- gsub("\\..", "", names(chr_arm_bins))
  chr_bins <- lapply(split(chr_arm_bins, chr_names), unlist)
  keep_idx <- lapply(lapply(chr_bins, names), grep, pattern = "q1$", invert = T)
  chr_bins <- mapply( `[`, chr_bins, keep_idx, SIMPLIFY = F)

  #get arms for each bin
  keep_idx2 <- lapply(lapply(chr_bins, names), grep, pattern = "p1$", invert = T)
  chr_bin_arms <- lapply(lapply(chr_bins, names), gsub, pattern = "[0-9]+$", replacement = "")
  chr_bin_arms <- mapply(`[`, chr_bin_arms, keep_idx2)
  chr_bin_arms <- unlist(chr_bin_arms, use.names = F)

  chr_bin_tabs <- mapply(get_counts, atac_chr_list, chr_bins[names(atac_chr_list)], SIMPLIFY = F)
  chr_bin_tabs <- unlist(chr_bin_tabs, recursive = F)
  chr_bin_tabs <- lapply(chr_bin_tabs, `[`, barcodes)

  count_mat <- do.call(rbind, chr_bin_tabs)
  colnames(count_mat) <- barcodes
  rownames(count_mat) <- paste0(chr_bin_arms, ":", gsub(".*\\.", "", rownames(count_mat)))
  count_mat[is.na(count_mat)] <- 0
  count_mat
}

get_atac_barcode_stats <- function(cr_arc_barcode_stats_path, cr_atac_barcode_stats_path, barcodes){
  arc.barcode.info <- read.delim(cr_arc_barcode_stats_path, header = T, sep = ",")
  atac.barcode.info <- read.delim(cr_atac_barcode_stats_path, header = T, sep = ",")

  rownames(arc.barcode.info) <- arc.barcode.info$barcode
  arc.barcode.info <- arc.barcode.info[barcodes,]

  rownames(atac.barcode.info) <- atac.barcode.info$barcode
  atac.barcode.info <- atac.barcode.info[arc.barcode.info$atac_barcode,]

  atac.barcode.info
}

#normalise_counts <- function(x, ){
#
#}

simulate_counts <- function(x, overdispersion){
  lib.sizes <- colSums(x)
  gene.props <- rowSums(x) / sum(x)
  mu.mat <- t(lib.sizes %*% t(gene.props))
  mu.list <- as.list(as.data.frame(mu.mat))

  sim.list <- list()
  if (overdispersion == 0){
    for (i in seq_along(mu.list)){
      sim.list[[i]] <- lapply(mu.list[[i]], rpois, n = 1)
    }
  } else {
    for (i in seq_along(mu.list)){
      sim.list[[i]] <- lapply((1 / overdispersion) / (mu.list[[i]] + 1 / overdispersion), rnbinom, n = 1, size = 1 / overdispersion)
    }
  }
  sim.mat <- do.call(cbind, lapply(sim.list, as.numeric))
  rownames(sim.mat) <- rownames(x)
  colnames(sim.mat) <- colnames(x)
  sim.mat
}

permute_rows <- function(x){
  x.permute <- do.call(rbind, lapply(as.list(as.data.frame(t(x))), sample))
  rownames(x.permute) <- rownames(x)
  colnames(x.permute) <- colnames(x)
  x.permute
}

choose_pcs <- function(sdev_test, sdev_null, tolerance){
  which(sapply(lapply(sdev_test, `<`, sdev_null), sum) <= tolerance)
}

get_internal_edges <- function(knn_graph, leiden_clusters){
  cluster.ids <- sort(unique(leiden_clusters))
  subgraph.list <- list()
  for (i in cluster.ids){
    subgraph.list[[i]] <- subgraph(knn_graph, which(leiden_clusters == i))
  }
  sum(sapply(subgraph.list, ecount))
}

sweep_leiden_resolution <- function(knn_graph, objective_function, leiden_resolutions, n_reps, n_iter){
  require("igraph")
  leiden_list <- list()
  for (i in seq_along(leiden_resolutions)){
    leiden_list[[i]] <- list()
    for (j in (1:n_reps)){
      leiden_list[[i]][[j]] <- cluster_leiden(knn_graph, objective_function = objective_function, resolution_parameter = leiden_resolutions[i], n_iterations = n_iter)
    }
  }
  member_list <- lapply(leiden_list, lapply, membership)
  edge_list <- lapply(member_list, sapply, get_internal_edges, knn_graph = knn_graph)
  edge_list
}

choose_best_leiden_partition <- function(knn_graph, objective_function, leiden_resolution, n_reps, n_iter){
  require("igraph")
  leiden_list <- list()
  for (i in (1:n_reps)){
    leiden_list[[i]] <- cluster_leiden(knn_graph, objective_function = objective_function, resolution_parameter = leiden_resolution, n_iterations = n_iter)
  }
  leiden_quality <- sapply(leiden_list, function(x){x$quality})
  best_idx <- which(order(leiden_quality, decreasing = T) == 1)
  leiden_clusters <- membership(leiden_list[[best_idx]])
  sorted_tabs <- sort(table(leiden_clusters), decreasing = T)
  new_cluster_ids <- setNames(1:length(sorted_tabs), names(sorted_tabs))
  factor(unname(new_cluster_ids[as.character(leiden_clusters)]), levels = 1:length(sorted_tabs))
}

get_barcode_probability <- function(barcodes, weights = rep(1, length(barcodes))){
  barcode_1_9 <- substr(colnames(peak_counts), 1, 9)
  barcode_1_9_tabs <- table(rep(barcode_1_9, weights))
  barcode_1_9_probs <- (barcode_1_9_tabs / sum(barcode_1_9_tabs))[barcode_1_9]
  barcode_8_16 <- substr(colnames(peak_counts), 8, 16)
  barcode_8_16_tabs <- table(rep(barcode_8_16, weights))
  barcode_8_16_probs <- (barcode_8_16_tabs / sum(barcode_8_16_tabs))[barcode_8_16]
  setNames(barcode_1_9_probs * barcode_8_16_probs, barcodes)
}

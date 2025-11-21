get_copy_number <- function(count_matrix, clusters, ref.cluster, joint, is_ref_female, ploidy.correction.method = "prop", external.ref = NULL){
  #require("MatrixGenerics")
  #removed bias correction - need to investigate whether it improves calling. Scaling to mean in presence of zeroes may be a problem.
  cor.list <- list()
  if (is.null(ref.cluster) & !is.null(external.ref)){
    ref.cluster <- "external"
    if (ploidy.correction.method == "prop"){
      cor.list[["external"]] <- rowMeans(external.ref) / sum(external.ref)
    }
  }
  for (i in levels(clusters)){
    cor.list[[i]] <- rowMeans(count_matrix[,as.integer(clusters) == i])
  }
  cn.mat <- count_matrix
  if (ploidy.correction.method == "prop"){
    for (i in seq_along(cor.list)){
      cor.list[[i]] <- cor.list[[i]] / sum(cor.list[[i]])
    }
    if (!joint){
      cn.mat <- t(t(count_matrix) / colSums(count_matrix))
    }
  }
  if (joint){
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- rep(2 * as.numeric(cor.list[[i]]) / as.numeric(cor.list[[ref.cluster]]), sum(as.integer(clusters) == i))
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- ifelse(is_ref_female, 1, 0.5) * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  } else {
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- 2 * cn.mat[,as.integer(clusters) == i] / as.numeric(cor.list[[ref.cluster]])
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- ifelse(is_ref_female, 1, 0.5) * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  }
  cn.mat[is.infinite(cn.mat)] <- NA
  cn.mat
}

get_new_order <- function(cn.mat, count_matrix_cor, clusters){
  #sort within clusters based on corrected count matrix
  hclust.list <- list()
  for (i in levels(clusters)){
    cell.dist <- dist(t(count_matrix_cor[,as.integer(clusters) == i]))
    hclust.list[[i]] <- colnames(cn.mat)[as.integer(clusters) == i][hclust(cell.dist, method = "ward.D2")$order]
  }
  #perform hierarchical clustering on per-cluster copy-number estimates
  #get first cell for each cluster (they all have the same estimate)
  cluster.idx <- setNames(integer(length(levels(clusters))), levels(clusters))
  for (cluster in levels(clusters)){
    cluster.idx[cluster] <- which(as.integer(clusters) == cluster)[1]
  }
  cluster.dist <- dist(t(cn.mat[,cluster.idx]), method = "manhattan")
  cluster.hclust <- hclust(cluster.dist, method = "complete")
  cluster.order <- names(cluster.idx)[cluster.hclust$order]
  cluster.order <- c(ref.cluster, cluster.order[cluster.order != ref.cluster])
  new.cluster.levels <- factor(clusters[unlist(hclust.list[cluster.order])], levels = cluster.order)

  chr.num <- as.integer(gsub("chr", "", sapply(strsplit(rownames(cn.mat), "\\."), `[`, 1)))
  new.rownames <- c(rownames(cn.mat)[order(chr.num)[!is.na(chr.num)]], rownames(cn.mat)[is.na(chr.num)])

  new.row.levels <- gsub("\\.", "", gsub("chr", "", sapply(strsplit(new.rownames, ":"), `[`, 1)))
  new.row.levels <- factor(new.row.levels, unique(new.row.levels))

  list(new.rownames, unlist(hclust.list), new.cluster.levels, new.row.levels)
}

get_new_feature_names <- function(feature.names){
  chr.num <- as.integer(gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1)))
  new.rownames <- c(feature.names[order(chr.num)[!is.na(chr.num)]], feature.names[is.na(chr.num)])

  new.rownames
}

get_feature_factor <- function(feature.names){
  #chr.num <- as.integer(gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1)))

  #new.row.levels <- gsub("\\.", "", gsub("chr", "", sapply(strsplit(feature.names, ":"), `[`, 1)))
  #new.row.levels <- factor(new.row.levels, unique(new.row.levels))
  #new.row.levels
  new.row.levels <- gsub("chr", "", sapply(strsplit(feature.names, "\\."), `[`, 1))
  factor(new.row.levels, unique(new.row.levels))
}

get_new_cell_order <- function(cn.mat, count_matrix_cor, clusters){
  #sort within clusteres based on corrected count matrix
  hclust.list <- list()
  for (i in levels(clusters)){
    cell.dist <- dist(t(count_matrix_cor[,as.integer(clusters) == i]))
    hclust.list[[i]] <- colnames(cn.mat)[as.integer(clusters) == i][hclust(cell.dist, method = "ward.D2")$order]
  }
  unlist(hclust.list)
}

get_cell_factor <- function(cn.mat, clusters){
  #perform hierarchical clustering on per-cluster copy-number estimates
  #get first cell for each cluster (they all have the same estimate)
  cluster.idx <- setNames(integer(length(levels(clusters))), levels(clusters))
  for (cluster in levels(clusters)){
    cluster.idx[cluster] <- which(as.integer(clusters) == cluster)[1]
  }
  cluster.dist <- dist(t(cn.mat[,cluster.idx]), method = "manhattan")
  cluster.hclust <- hclust(cluster.dist, method = "complete")
  cluster.order <- names(cluster.idx)[cluster.hclust$order]
  #cluster.order <- c(ref.cluster, cluster.order[cluster.order != ref.cluster])
  new.cluster.levels <- factor(clusters, levels = as.integer(cluster.order))
  new.cluster.levels
}

get_chr_arms <- function(bin.names){
  sapply(strsplit(sapply(strsplit(bin.names, ":"), `[`, 1), "\\."), `[`, 2)
}

reorder_matrix <- function(){}

#cell.list <- as.list(as.data.frame(peak_counts_filtered[new_feature_names,]))
#cell.arm.count <- do.call(cbind, lapply(cell.list, get_arm_counts, feature_factor, get_chr_arms(new_feature_names)))
get_arm_counts <- function(x, chr.factor, arm.factor){
  sapply(split(x, list(chr.factor, arm.factor)), sum)
}

get_purple_cn <- function(purple_path, bin.names){
  purple.cn <- read.delim(purple_path, header = T)
  purple.cn.granges <- makeGRangesFromDataFrame(purple.cn)
  ataclone.granges <- GRanges(gsub("\\..", "", bin.names))
  cn.overlaps <- findOverlaps(ataclone.granges, purple.cn.granges)
  cn.overlaps.split <- split(cn.overlaps@to, cn.overlaps@from)
  cn.overlaps.widths <- list()
  for (i in seq_along(cn.overlaps.split)){
    cn.overlaps.widths[[i]] <- numeric(length = length(cn.overlaps.split[[i]]))
    for (j in seq_along(cn.overlaps.widths[[i]])){
      #cn.overlaps.widths[[i]][j] <- width(GenomicRanges::intersect(ataclone.granges[i], purple.cn.granges[cn.overlaps.split[[i]][j]]))
      cn.overlaps.widths[[i]][j] <- min(end(ataclone.granges[i]), end(purple.cn.granges[cn.overlaps.split[[i]][j]])) - max(start(ataclone.granges[i]), start(purple.cn.granges[cn.overlaps.split[[i]][j]])) + 1
    }
  }
  ataclone.purple.estimates <- numeric(length(cn.overlaps.split))
  for (i in seq_along(ataclone.purple.estimates)){
    ataclone.purple.estimates[i] <- sum(purple.cn[cn.overlaps.split[[i]],"copyNumber"] * cn.overlaps.widths[[i]] / sum(cn.overlaps.widths[[i]]))
  }
  ataclone.purple.estimates
}

#messy - need to fix
get_new_hclust <- function(x, new_cell_order, clusters){
  test_hclust <- hclust(dist(t(x[,new_cell_order]), method = "manhattan"), method = "single")
  test_clusters <- cutree(test_hclust, length(levels(clusters)))
  idx_split <- split(1:length(clusters), test_clusters[new_cell_order])

  get_cluster_merges <- function(cluster_indices, cluster_number){
    idx.split <- split(cluster_indices, seq_along(cluster_indices) > 3)
    #js are merged cluster indices
    js <- (cluster_indices[1] - cluster_number + 1):(cluster_indices[length(cluster_indices)] - cluster_number)
    #go along first half in order
    left.idx <- integer(length(js))
    right.idx <- integer(length(js))
    left.idx[1] <- -idx.split[[1]][1]
    right.idx[1] <- -idx.split[[1]][2]
    for (i in 2:(length(idx.split[[1]]) - 1)){
      left.idx[i] <- -(cluster_indices[i] + 1)
      right.idx[i] <- js[i - 1]
    }
    left.final.j <- right.idx[i] + 1
    #go along second half in order
    left.idx[length(idx.split[[1]])] <- -idx.split[[2]][2]
    right.idx[length(idx.split[[1]])] <- -idx.split[[2]][1]
    for (i in (length(idx.split[[1]]) + 1):(length(cluster_indices) - 2)){
      left.idx[i] <- -(cluster_indices[i] + 2)
      right.idx[i] <- js[i - 1]
    }
    right.final.j <- right.idx[i] + 1
    left.idx[i + 1] <- right.final.j
    right.idx[i + 1] <- left.final.j
    cbind(left.idx, right.idx)
  }

  merge_list <- mapply(get_cluster_merges, idx_split, 1:length(levels(clusters)), SIMPLIFY = F)
  #can't remember what this is
  #subgraph_indices <- sapply(merge_list, function(x){x[nrow(x),2] + 1})
  #make graphs of pseudo-bulk
  use.clusters <- 1:length(levels(clusters))
  cluster.cn <- list()
  for (i in seq_along(use.clusters)){
    cluster.cn[[i]] <- x[,clusters == use.clusters[i]][,1]
  }
  cluster.cn <- do.call(cbind, cluster.cn)
  joint_hclust <- hclust(dist(t(cluster.cn), method = "manhattan"), "single")
  get_final_rows <- function(merge_list, joint_hclust){
    subgraph.indices <- sapply(merge_list, max) + 1
    l.indices <- integer(length(subgraph.indices) - 1)
    r.indices <- integer(length(subgraph.indices) - 1)
    merge.idx <- max(subgraph.indices) + 1
    #use to track merge indices as they are created
    new.merge.indices <- integer(length(subgraph.indices) - 1)
    #use j to write new indices
    j <- 1
    #use k to retrieve new indices after writing based off joint hclust
    k <- 1
    for (i in 1:(length(subgraph.indices) - 1)){
      new.merge.indices[j] <- merge.idx
      l.idx <- joint_hclust$merge[i,1]
      r.idx <- joint_hclust$merge[i,2]
      if (l.idx < 0){
        l.indices[i] <- subgraph.indices[-l.idx]
      } else {
        l.indices[i] <- new.merge.indices[k]
        k <- k + 1
        #merge.idx <- merge.idx + 1
      }
      if (r.idx < 0){
        r.indices[i] <- subgraph.indices[-r.idx]
      } else {
        r.indices[i] <- new.merge.indices[k]
        k <- k + 1
        #merge.idx <- merge.idx + 1
      }
      j <- j + 1
      merge.idx <- merge.idx + 1
    }
    cbind(l.indices, r.indices)
  }
  #final.rows <- cbind(c(386,509,596,233,558,137,650,646,311), c(448,622,645,648,649,647,651,652,653))
  #new_merge2 <- rbind(do.call(rbind, merge_list), final.rows)
  new_merge2 <- rbind(do.call(rbind, merge_list), get_final_rows(merge_list, joint_hclust))
  test_hclust$merge <- new_merge2
  test_hclust
}

#' @export
plot_copy_number <- function(x, external_ref, clusters, pca_obj, discard_pcs, is_ref_female){
  use_pcs <- 1:ncol(pca_obj$x)
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]

  single_cn_estimates <- get_copy_number(x, clusters, ref.cluster = NULL, joint = F, external.ref = external_ref, is_ref_female = is_ref_female)
  joint_cn_estimates <- get_copy_number(x, clusters, ref.cluster = NULL, joint = T, external.ref = external_ref, is_ref_female = is_ref_female)
  #new_orders <- get_new_order(joint_cn_estimates, stable_counts_filtered_norm_cor, leiden_clusters)
  new_feature_names <- get_new_feature_names(rownames(joint_cn_estimates))
  feature_factor <- get_feature_factor(new_feature_names)

  x_norm_cor_denoised <- t(pca_obj$x[,use_pcs] %*% t(pca_obj$rotation[,use_pcs]))
  new_cell_order <- get_new_cell_order(joint_cn_estimates, x_norm_cor_denoised, leiden_clusters)
  cell_factor <- get_cell_factor(joint_cn_estimates[,new_cell_order], clusters[new_cell_order])
  chr_arm <- ComplexHeatmap::HeatmapAnnotation(chr_arm = get_chr_arms(new_feature_names), col = list(chr_arm = c(p = "black", q = "grey")), show_legend = F)

  color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))

  use.clusters <- 1:length(levels(clusters))
  cluster.cn <- list()
  for (i in seq_along(use.clusters)){
    cluster.cn[[i]] <- joint_cn_estimates[,clusters == use.clusters[i]][,1]
  }
  cluster.cn <- do.call(cbind, cluster.cn)
  test_hclust <- get_new_hclust(joint_cn_estimates, new_cell_order, clusters)

  ComplexHeatmap::Heatmap(t(as.matrix(single_cn_estimates[new_feature_names,][,new_cell_order])), cluster_rows = test_hclust, cluster_columns = F, column_split = feature_factor, row_labels = rep('', ncol(joint_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, col = color_fun, top_annotation = chr_arm, row_split = length(levels(leiden_clusters)), row_title = hclust(dist(t(cluster.cn), method = "manhattan"), "single")$order, heatmap_legend_param = list(title = "copy_number"))
}

#' @export
get_absolute_copy_number <- function(x, clusters, ref_cluster, is_outlier, is_doublet, is_ref_female){
  new_feature_names <- get_new_feature_names(rownames(x))
  cn.list <- list()
  for (i in levels(clusters)[levels(clusters) != ref_cluster]){
    cn <- rowMeans(x[new_feature_names, clusters == i & !is_outlier & !is_doublet]) /
      rowMeans(x[new_feature_names, clusters == ref_cluster & !is_outlier & !is_doublet])
    chr <- gsub("chr", "", sapply(strsplit(new_feature_names, "\\."), `[`, 1))
    chr_arm <- get_chr_arms(new_feature_names)
    chr <- factor(chr, unique(chr))
    cn <- ifelse(grepl("X|Y", chr) & !is_ref_female, 1, 2) * cn
    cn.list[[i]] <- cn
  }
  cn <- rowMeans(x[new_feature_names, clusters != ref_cluster & !is_outlier & !is_doublet]) /
    rowMeans(x[new_feature_names, clusters == ref_cluster & !is_outlier & !is_doublet])
  chr <- factor(chr, unique(chr))
  cn <- ifelse(grepl("X|Y", chr) & !is_ref_female, 1, 2) * cn
  cn.list[["all"]] <- cn
  cn.list
}

#' @export
plot_absolute_cn <- function(x, cluster_name){
  chr <- gsub("chr", "", sapply(strsplit(names(x), "\\."), `[`, 1))
  chr <- factor(chr, unique(chr))
  new_xs <- unlist(sapply(split(names(x), chr), function(x){1:length(x)}))
  plot(ggplot(data.frame(cn = x, chr, new_xs), aes(x =new_xs, y = cn, color = chr)) + geom_point() + facet_grid(cols = vars(chr), space = "free", scales = "free") + theme(panel.spacing = unit(0, "lines"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.background = element_rect(color = "white"), strip.background=element_rect(color = "white"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none") + scale_x_continuous(expand = expansion(add = 2)) + ylim(c(0,8)) + ggtitle(paste0("Absolute copy number (Cluster = ", cluster_name, ")")))
}

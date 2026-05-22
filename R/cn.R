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
  new_cell_order <- get_new_cell_order(joint_cn_estimates, x_norm_cor_denoised, clusters)
  cell_factor <- get_cell_factor(joint_cn_estimates[,new_cell_order], clusters[new_cell_order])
  chr_arm <- ComplexHeatmap::HeatmapAnnotation(chr_arm = get_chr_arms(new_feature_names), col = list(chr_arm = c(p = "black", q = "grey")), show_legend = F)

  color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))

  #overwrite ComplexHeatmap::cluster_between_groups to allow different distances and linkages
  cluster_between_groups <- function(mat, factor, dist.method, linkage)
  {
    if (!is.factor(factor)) {
      factor = factor(factor, levels = unique(factor))
    }
    dend_list = list()
    order_list = list()
    for (le in unique(levels(factor))) {
      m = mat[, factor == le, drop = FALSE]
      if (ncol(m) == 1) {
        order_list[[le]] = which(factor == le)
        dend_list[[le]] = structure(which(factor == le),
                                    class = "dendrogram", leaf = TRUE, height = 0,
                                    label = 1, members = 1)
      }
      else if (ncol(m) > 1) {
        hc1 = hclust(dist(1:ncol(m)))
        dend_list[[le]] = reorder(as.dendrogram(hc1), wts = 1:ncol(m),
                                  agglo.FUN = mean)
        order_list[[le]] = which(factor == le)[order.dendrogram(dend_list[[le]])]
        dend_list[[le]] <- ComplexHeatmap:::`order.dendrogram<-`(dend_list[[le]], order_list[[le]])
      }
      attr(dend_list[[le]], ".class_label") = le
    }
    parent = as.dendrogram(hclust(dist(t(sapply(order_list, function(x) rowMeans(mat[,
                                                                                     x, drop = FALSE]))), dist.method), linkage))
    dend_list = lapply(dend_list, function(dend) dendrapply(dend,
                                                            function(node) {
                                                              attr(node, "height") = 0
                                                              node
                                                            }))
    dend = ComplexHeatmap::merge_dendrogram(parent, dend_list)
    dend <- ComplexHeatmap:::`order.dendrogram<-`(dend, unlist(order_list[order.dendrogram(parent)]))
    return(dend)
  }

  #https://github.com/jokergoo/ComplexHeatmap/issues/694
  dend <- cluster_between_groups(as.matrix(single_cn_estimates[new_feature_names,][,new_cell_order]), clusters[new_cell_order], "manhattan", "single")

  ComplexHeatmap::Heatmap(t(as.matrix(single_cn_estimates[new_feature_names,][,new_cell_order])), row_split = length(levels(clusters)), cluster_columns = F, column_split = feature_factor, row_labels = rep('', ncol(single_cn_estimates)), column_labels = rep('', nrow(joint_cn_estimates)), border = T, col = color_fun, top_annotation = chr_arm, heatmap_legend_param = list(title = "copy_number"), cluster_rows = dend)
}

#' @export
get_absolute_copy_number <- function(x, clusters, ref_cluster, scale_factors, is_excluded, is_ref_female){
  new_feature_names <- get_new_feature_names(rownames(x))
  x <- x[,!is_excluded]
  #hacky way to add support for multiple reference clusters
  if (length(ref_cluster) > 1){
    clusters[clusters %in% ref_cluster & !(clusters == ref_cluster[1])] <- ref_cluster[1]
    ref_cluster <- ref_cluster[1]
  }
  clusters <- factor(clusters[!is_excluded])
  cn.list <- list()
  for (i in levels(clusters)[levels(clusters) != ref_cluster]){
    sf <- scale_factors[paste0("cluster_", i)]
    cn <- sf * rowMeans(x[new_feature_names, clusters == i]) /
      rowMeans(x[new_feature_names, clusters == ref_cluster])
    chr <- gsub("chr", "", sapply(strsplit(new_feature_names, "\\."), `[`, 1))
    chr_arm <- get_chr_arms(new_feature_names)
    chr <- factor(chr, unique(chr))
    cn <- ifelse(grepl("X|Y", chr) & !is_ref_female, 1, 2) * cn
    cn.list[[i]] <- cn
  }
  cn <- colMeans(scale_factors[paste0("cluster_", clusters[clusters != ref_cluster])] * t(x[new_feature_names, clusters != ref_cluster])) /
    rowMeans(x[new_feature_names, clusters == ref_cluster])
  chr <- gsub("chr", "", sapply(strsplit(new_feature_names, "\\."), `[`, 1))
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

#' @export
plot_absolute_cn2 <- function(absolute_cn_list, max_dif = 1.5, y_lim = c(0,8), exclude_clusters = c("all")){
  is.empty <- sapply(lapply(absolute_cn_list, is.na), sum) == lengths(absolute_cn_list)
  absolute_cn_list <- absolute_cn_list[!is.empty]
  cn_mat <- do.call(cbind, absolute_cn_list)
  test_plt <- reshape2::melt(cn_mat)
  colnames(test_plt) <- c("bin", "clone", "cn")
  #or, can pick a reference cluster
  test_plt$cn_resid <- as.vector(cn_mat - absolute_cn_list$all)
  test_plt$cn_resid[test_plt$cn_resid < -max_dif] <- -max_dif
  test_plt$cn_resid[test_plt$cn_resid > max_dif] <- max_dif
  test_plt <- test_plt[!test_plt$clone %in% exclude_clusters,]
  test_plt$idx <- setNames(1:313, names(absolute_cn_list$`1`))[test_plt$bin]
  chr <- gsub("chr", "", sapply(strsplit(as.character(test_plt$bin), "\\."), `[`, 1))
  chr <- factor(chr, unique(chr))
  test_plt$chr <- chr
  plot(ggplot(test_plt, aes(x = idx, y = cn, color = cn_resid)) + geom_point() + facet_grid(clone ~ chr, space = "free", scales = "free") + theme(panel.spacing.x = unit(0, "lines"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.background = element_rect(color = "grey92", fill = "white"), panel.grid.major.y = element_line(color = "grey92"), panel.grid.minor.y = element_line(color = "grey92"), strip.background=element_rect(color = "white"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none") + scale_x_continuous(expand = expansion(add = 2)) +
         scale_color_gradient2(low = "blue", mid =  "grey", high = "red", midpoint = 0, limits = c(-max_dif, max_dif)) + ylim(y_lim))
}

get_mse_difs <- function(x, y, sfs, bin_weights, bootstrap_vars){
  mses <- numeric(length(sfs))
  for (i in seq_along(sfs)){
    sf <- sfs[i]
    dif.rounded <- round(sf * x - y)
    precision_weights <- 1 / bootstrap_vars[[i]]
    #precision weights are causing higher scaling factors to be favoured - fixed?
    mses[i] <- weighted.mean((sf * x - y - dif.rounded) ** 2, bin_weights * precision_weights)
  }
  #cn.estimates <- lapply(sfs, `*`, ratio)
  #is.excluded <- sapply(cn.estimates, median) < 1.5
  #mses[is.excluded] <- NA
  mses
}

get_mse <- function(ratio, sfs1, sfs2, bin_weights){
  sfs <- sfs1 / sfs2
  mses <- numeric(length(sfs))
  for (i in seq_along(sfs)){
    sf <- sfs[i]
    ratio.rounded <- round(sf * ratio)
    mses[i] <- weighted.mean((sf * ratio - ratio.rounded) ** 2, bin_weights)
  }
  cn.estimates <- lapply(sfs, `*`, ratio)
  is.excluded <- sapply(cn.estimates, median) < 1.5
  mses[is.excluded] <- NA
  mses
}

bootstrap_cluster_means <- function(x, n_bootstrap){
  bootstrap_means <- list()
  for (i in 1:n_bootstrap){
    bootstrap_means[[i]] <- rowMeans(x[,sample(1:ncol(x), replace = T)])
  }
  do.call(cbind, bootstrap_means)
}

get_mse_x_ref <- function(x, clusters, bootstrap_mean_list, x_idx, ref_cluster, is_ref_female, n_bootstrap, sfs, weight_by_precision, weight_by_cluster){
  x_mean <- rowMeans(x[,clusters == x_idx])
  ref_mean <- rowMeans(x[,clusters == ref_cluster])
  x_ref_bootstrap <- bootstrap_mean_list[[x_idx]] / bootstrap_mean_list[[ref_cluster]]
  if (!is_ref_female){
    x_ref_bootstrap[grepl("chrX|chrY", rownames(x_ref_bootstrap))] <- 0.5 * x_ref_bootstrap[grepl("chrX|chrY", rownames(x_ref_bootstrap))]
  }
  x_ref_bootstrap[is.na(x_ref_bootstrap)] <- Inf
  x_ref_var <- rowVars(x_ref_bootstrap)
  is_excluded_x <- x_ref_var == 0 | rowSums(is.infinite(x_ref_bootstrap) > 0) | ref_mean < 2
  x_ref_ratio <- 2 * (x_mean / ref_mean)
  if (!is_ref_female){
    x_ref_ratio[grepl("chrX|chrY", names(x_ref_ratio))] <- 0.5 * x_ref_ratio[grepl("chrX|chrY", names(x_ref_ratio))]
  }

  x_ref_precision_weights <- 1 / x_ref_var[!is_excluded_x]
  #sex chromosomes in males have higher precision weights due to lower scaling factor than autosomes - let's undo for computing mse
  if (!is_ref_female){
    x_ref_precision_weights[grepl("chrX|chrY", names(x_ref_precision_weights))] <- 0.5 * x_ref_precision_weights[grepl("chrX|chrY", names(x_ref_precision_weights))]
  }
  x_ref_precision_weights <- length(x_ref_precision_weights) * x_ref_precision_weights / sum(x_ref_precision_weights)
  x_ref_precision_weights2 <- 1 / rowVars(log(x_ref_bootstrap[!is_excluded_x,] + 0.01))
  x_ref_precision_weights2 <- length(x_ref_precision_weights2) * x_ref_precision_weights2 / sum(x_ref_precision_weights2)

  if (weight_by_cluster){
    x_clusters <- Ckmeans.1d.dp::Ckmeans.1d.dp(log(x_ref_ratio[!is_excluded_x] + 0.01), y = x_ref_precision_weights2)
    x_cluster_weights <- (1 / table(x_clusters$cluster))[as.character(x_clusters$cluster)]
  } else{
    x_clusters <- list()
    x_clusters$cluster <- rep(1,length(x_ref_ratio[!is_excluded_x]))
    x_cluster_weights <- rep(1,length(x_ref_ratio[!is_excluded_x]))
  }

  if (weight_by_precision){
    x_ref_weights <- x_cluster_weights * x_ref_precision_weights
  } else {
    x_ref_weights <- x_cluster_weights
  }

  x_mse <- get_mse(x_ref_ratio[!is_excluded_x], sfs, 1, x_ref_weights)
  cn_clusters <- rep(NA, length(is_excluded_x))
  cn_clusters[!is_excluded_x] <- x_clusters$cluster
  cluster_weights <- rep(NA, length(is_excluded_x))
  cluster_weights[!is_excluded_x] <- x_cluster_weights
  list(mses = setNames(x_mse, sfs), ratio = x_ref_ratio, bootstrap_mat = x_ref_bootstrap, is_excluded = is_excluded_x, clusters = cn_clusters, cluster_weights = cluster_weights)
}

#' Joint MSE fit of cluster scaling factors for absolute copy number
#'
#' Allocates a `length(sfs) ^ (n_clusters - 1)` MSE array; this scales
#' geometrically and can trivially exceed available memory. The
#' `max_array_size` argument (default 1e8 entries, ~0.8 GB at 8 bytes
#' each) guards the allocation and raises an actionable error before R
#' attempts an out-of-memory allocation. If you hit the guard, either
#' coarsen `sf_step` (increase the step → fewer grid points) or
#' reduce the cluster count fed in.
#'
#' @param max_array_size integer; maximum allowed entry count for the
#'   joint MSE array (`length(sfs) ^ (n_clusters - 1)`). Set to `Inf` to
#'   disable the guard. Default 1e8.
#' @export
fit_cn_clusters <- function(x, clusters, ref_cluster, is_ref_female, is_excluded, n_bootstrap, sf_start, sf_end, sf_step, weight_by_precision, weight_by_cluster, fit_clusters_pairwise = T, weighting_strategy = "min_size", max_array_size = 1e8){
  x <- x[, !is_excluded]
  #hacky way to add support for multiple reference clusters
  if (length(ref_cluster) > 1){
    clusters[clusters %in% ref_cluster & !(clusters == ref_cluster[1])] <- ref_cluster[1]
    ref_cluster <- ref_cluster[1]
  }
  clusters <- factor(clusters[!is_excluded])
  sfs <- seq(sf_start, sf_end, sf_step)
  n_levels <- length(levels(clusters))
  if (n_levels < 2) {
    stop("fit_cn_clusters needs at least 2 cluster levels (a reference + ",
         ">=1 non-reference); got ", n_levels, ".")
  }
  array_entries <- length(sfs) ^ (n_levels - 1)
  if (is.finite(max_array_size) && array_entries > max_array_size) {
    stop(sprintf(
      paste0("fit_cn_clusters MSE array would have %g entries (%d sfs ^ ",
             "(%d clusters - 1)), exceeding max_array_size = %g. ",
             "Coarsen sf_step (e.g. %.2f -> %.2f), reduce the cluster ",
             "count, or pass max_array_size = Inf to override."),
      array_entries, length(sfs), n_levels, max_array_size,
      sf_step,
      max(sf_step, (sf_end - sf_start) / (max_array_size ^ (1 / (n_levels - 1)) - 1))
    ))
  }
  bootstrap_mean_list <- list()
  for (i in seq_along(levels(clusters))){
    bootstrap_mean_list[[i]] <- bootstrap_cluster_means(x[,clusters == levels(clusters)[i]], n_bootstrap)
  }
  names(bootstrap_mean_list) <- levels(clusters)
  mse_array <- array(0,dim = rep(length(sfs), length(levels(clusters)) - 1), dimnames = lapply(1:(length(levels(clusters)) - 1), function(x){sfs}))
  #first, compute mses for each cluster against the reference (distance from rounded values)
  non_ref_clusters <- levels(clusters)[levels(clusters) != ref_cluster]
  ref_cluster <- levels(clusters)[levels(clusters) == ref_cluster]
  cluster_tabs <- table(clusters)
  if (weighting_strategy == "min_size"){
    cluster_weights <- cluster_tabs
  } else{
    cluster_weights <- rep(1, length(non_ref_clusters))
  }
  x_ref_out_list <- list()
  #need to check why this doesn't produce same results as before - fixed?
  print(paste0("Fitting ", length(non_ref_clusters), " clusters individually"))
  for (i in seq_along(non_ref_clusters)){
    print(paste0("Fitting cluster ", non_ref_clusters[i], " (", i, "/", length(non_ref_clusters), ")"))
    new_mse_array <- array(,dim = rep(length(sfs), length(levels(clusters)) - 1), dimnames = lapply(1:(length(levels(clusters)) - 1), function(x){sfs}))
    #also returns bootstrap values to re-use later for pairwise comparisons
    x_ref_out_list[[i]] <- get_mse_x_ref(x, clusters, bootstrap_mean_list, non_ref_clusters[[i]], ref_cluster, is_ref_female, n_bootstrap, sfs, weight_by_precision, weight_by_cluster)
    x_ref_out_list[[i]]$mse <- cluster_weights[i] * x_ref_out_list[[i]]$mse
    mse_idx <- vector(mode = "list", length = length(non_ref_clusters))
    for (j in seq_along(mse_idx)){
      mse_idx[[j]] <- 1:length(sfs)
    }
    #hacky: https://stackoverflow.com/questions/71470426/r-is-there-a-way-to-index-an-array-without-knowing-the-number-of-its-dimensions
    for (j in seq_along(sfs)){
      mse_idx[[i]] <- j
      #mse_array <- do.call(`[<-`, c(list(mse_array), mse_idx, list(do.call(`[`, c(list(mse_array), mse_idx)) + x_ref_out_list[[i]]$mse[j])))
      new_mse <- x_ref_out_list[[i]]$mse[j]
      template <- quote(new_mse_array[] <- new_mse)
      expr <- template
      for (k in 1:length(mse_idx)){
        expr[[c(2,2+k)]] <- mse_idx[[k]]
      }
      eval(expr)
    }
    mse_array <- mse_array + new_mse_array
  }
  gc()
  names(x_ref_out_list) <- non_ref_clusters
  array_dimnames <- setNames(1:length(non_ref_clusters), non_ref_clusters)
  #then, add mses based on distance of pairwise differences (after scaling) from their rounded values (testing alignment of copy number differences between clusters to integers)
  #to improve time efficiency, try making separate arrays for each cluster_pair (of the same dimensions), so that only the assignment operator needs to be used
  #(might reduce copying of data), then add them together at the end
  #more memory-efficient to make a new array for every iteration then add iteratively
  #array_list <- list()
  #for (i in seq_along(cluster_pairs)){
  #  array_list[[i]] <- array(,dim = rep(length(sfs), length(levels(clusters)) - 1), dimnames = lapply(1:(length(levels(clusters)) - 1), function(x){sfs}))
  #}
  if (length(non_ref_clusters) > 1 & fit_clusters_pairwise){
    cluster_pairs <- as.list(as.data.frame(combn(levels(clusters)[levels(clusters) != ref_cluster], 2)))
    #reorder cluster pairs so that largest group is reference (unless it is the normal reference)
    for (i in seq_along(cluster_pairs)){
      cluster_pairs[[i]] <- cluster_pairs[[i]][order(cluster_tabs[cluster_pairs[[i]]], decreasing = F)]
    }
    print(paste0("Fitting ", length(cluster_pairs), " cluster pairs"))
    if (weighting_strategy == "min_size"){
      #first element in each pair is always the smaller cluster
      cluster_weights <- cluster_tabs[sapply(cluster_pairs, `[`, 1)]
    } else{
      cluster_weights <- rep(1, length(cluster_pairs))
    }
    for (i in seq_along(cluster_pairs)){
      print(paste0("Fitting cluster ", cluster_pairs[[i]][1], " to cluster ", cluster_pairs[[i]][2], " (", i, "/", length(cluster_pairs), ")"))
      new_mse_array <- array(,dim = rep(length(sfs), length(levels(clusters)) - 1), dimnames = lapply(1:(length(levels(clusters)) - 1), function(x){sfs}))
      #mse_mat <- matrix(nrow = length(sfs), ncol = length(sfs))
      x_idx <- cluster_pairs[[i]][1]
      y_idx <- cluster_pairs[[i]][2]

      x_mean <- rowMeans(x[,clusters == x_idx])
      y_mean <- rowMeans(x[,clusters == y_idx])
      ref_mean <- rowMeans(x[,clusters == ref_cluster])

      x_y_bootstrap <- (bootstrap_mean_list[[x_idx]] / bootstrap_mean_list[[ref_cluster]]) / (bootstrap_mean_list[[y_idx]] / bootstrap_mean_list[[ref_cluster]])
      x_y_bootstrap[is.na(x_y_bootstrap)] <- Inf

      x_ref_ratio <- 2 * (x_mean / ref_mean)
      y_ref_ratio <- 2 * (y_mean / ref_mean)
      if (!is_ref_female){
        x_ref_ratio[grepl("chrX|chrY", names(x_ref_ratio))] <- 0.5 * x_ref_ratio[grepl("chrX|chrY", names(x_ref_ratio))]
        y_ref_ratio[grepl("chrX|chrY", names(y_ref_ratio))] <- 0.5 * y_ref_ratio[grepl("chrX|chrY", names(y_ref_ratio))]
      }
      x_y_ratio <- x_ref_ratio / y_ref_ratio
      x_y_var <- rowVars(x_y_bootstrap)
      is_excluded_x_y <- x_y_var == 0 | rowSums(is.infinite(x_y_bootstrap) > 0) | y_mean < 2

      x_y_precision_weights2 <- 1 / rowVars(log(x_y_bootstrap[!is_excluded_x_y,] + 0.01))
      x_y_precision_weights2 <- length(x_y_precision_weights2) * x_y_precision_weights2 / sum(x_y_precision_weights2)

      if (weight_by_cluster){
        x_y_clusters <- Ckmeans.1d.dp::Ckmeans.1d.dp(log(x_y_ratio[!is_excluded_x_y] + 0.01), y = x_y_precision_weights2)
        x_y_cluster_weights <- (1 / table(x_y_clusters$cluster))[as.character(x_y_clusters$cluster)]
      } else{
        x_y_cluster_weights <- rep(1, length(x_y_ratio[!is_excluded_x_y]))
      }

      use_array_dims <- c(array_dimnames[x_idx], array_dimnames[y_idx])
      mse_idx <- vector(mode = "list", length = length(non_ref_clusters))
      for (j in seq_along(mse_idx)){
        mse_idx[[j]] <- 1:length(sfs)
      }
      for (j in seq_along(sfs)){
        mse_idx[[use_array_dims[2]]] <- j
        bootstrap_vars <- list()
        if (weight_by_precision){
          for (k in seq_along(sfs)){
            bootstrap_vars[[k]] <- rowVars(sfs[k] * x_ref_out_list[[x_idx]]$bootstrap_mat - sfs[j] * x_ref_out_list[[y_idx]]$bootstrap_mat)[!is_excluded_x_y]
          }
        } else{
          for (k in seq_along(sfs)){
            bootstrap_vars[[k]] <- rep(1, sum(!is_excluded_x_y))
          }
        }
        mse_difs <- cluster_weights[i] * get_mse_difs(x_ref_ratio[!is_excluded_x_y], sfs[j] * y_ref_ratio[!is_excluded_x_y], sfs, x_y_cluster_weights, bootstrap_vars)
        #need to rotate data?
        #mse_array <- do.call(`[<-`, c(list(mse_array), mse_idx, list(do.call(`[`, c(list(mse_array), mse_idx)) + mse_difs)))
        #very slow because of copying of data, nesting of data, and number of updates (n_cluster_pairs * n_sf ** 2 - before was just n_clusters * n_sf)
        for (k in seq_along(mse_difs)){
          mse_idx[[use_array_dims[1]]] <- k
          #mse_array <- do.call(`[<-`, c(list(mse_array), mse_idx, list(do.call(`[`, c(list(mse_array), mse_idx)) + mse_difs[k])))
          #mse_subarray <- do.call(`[`, c(list(mse_array), mse_idx))
          #template <- quote(mse_array[] <- mse_subarray + mse_difs[k])
          #try not assigning mse_subarray
          template <- quote(new_mse_array[] <- mse_difs[k])
          expr <- template
          for (l in 1:length(mse_idx)){
            expr[[c(2,2+l)]] <- mse_idx[[l]]
          }
          eval(expr)
        }
      }
      mse_array <- mse_array + new_mse_array
    }
  }
  #mse_melted <- reshape2::melt(mse_array)
  #colnames(mse_melted) <- c(paste0("cluster_", non_ref_clusters), "mse")
  #mse_melted[which(mse_melted$mse == min(mse_melted$mse, na.rm = T)),]
  #mse_array <- mse_array + Reduce(`+`, array_list)
  list(scale_factors = setNames(as.numeric(rownames(mse_array))[arrayInd(which.min(mse_array), dim(mse_array))], paste0("cluster_", non_ref_clusters)),
       mse_array = mse_array)
}

#' @export
scDblFinder_custom_function <- function(e, dims){
  e_norm <- normalise_counts(e, nb_overdispersion, pseudo_count, lambda)
  e_norm_pca <- get_pca(e_norm, dims)
  #need to use discard_pcs from global environment - scDblFinder can't pass to the custom function and its includePCs argument does not work properly
  e_norm_pca$x[,(1:dims)[!((1:dims) %in% discard_pcs[discard_pcs != 1])]]
}

#' @export
get_doublet_thresholds <- function(x, clusters, is_doublet){
  cluster_levels <- levels(clusters)
  thresholds <- numeric(length(cluster_levels))
  for (i in seq_along(cluster_levels)){
    is_cluster <- clusters == cluster_levels[i]
    x2 <- x[,is_cluster]
    is_doublet2 <- is_doublet[is_cluster]
    if (sum(is_doublet2) == sum(is_cluster)){
      thresholds[i] <- 0
      next
    }
    if (sum(is_doublet2) == 0){
      thresholds[i] <- Inf
      next
    }
    lib_sizes <- colSums(x2)
    quantiles <- quantile(lib_sizes, c(0.1, 0.9))
    try_thresholds <- seq(quantiles[1], quantiles[2], 100)
    f1s <- numeric(length(try_thresholds))
    for (j in seq_along(try_thresholds)){
      threshold <- try_thresholds[j]
      is_less <- lib_sizes < threshold
      tn <- sum(is_less & !is_doublet2)
      fn <- sum(is_less & is_doublet2)
      tp <- sum(!is_less & is_doublet2)
      fp <- sum(!is_less & !is_doublet2)
      precision <- tp / (tp + fp)
      recall <- tp / (tp + fn)
      f1 <- 2 * (precision * recall) / (precision + recall)
      f1s[j] <- f1
    }
    thresholds[i] <- try_thresholds[which.max(f1s)]
  }
  setNames(thresholds, cluster_levels)
}

#computes differences in log odds-ratios between segments
get_log_or_dist <- function(x, clusters, test_idx, ref_idx, pseudo_count){
  #test_idx <- which(clusters %in% test_cluster)
  #ref_idx <- which(clusters %in% ref_cluster)
  #as.matrix(dist(log((rowMeans(x[, test_idx] + pseudo_count) / rowMeans(x[, ref_idx] + pseudo_count))), upper = T))
  #need to add check for is_ref_female
  sf <- ifelse(grepl("chrX|chrY", rownames(x)), 1, 2)
  missing_num <- rowSums(x[, test_idx]) == 0
  missing_denom <- rowSums(x[, ref_idx]) == 0
  log_or <- log(sf * (rowMeans(x[, test_idx] + pseudo_count) / rowMeans(x[, ref_idx] + sf * pseudo_count)))
  #seems to give more consistent results if you don't set to NA
  #log_or[missing_denom | missing_num] <- NA
  outer(log_or, log_or, "-")
}

#bootstraps the differences in log odds-ratios between segments - use to estimate variance of log-odds ratio differences
bootstrap_log_or_dist <- function(x, clusters, test_cluster, ref_cluster, pseudo_count, n_bootstrap = 100){
  dist_list <- list()
  test_idx <- which(clusters %in% test_cluster)
  ref_idx <- which(clusters %in% ref_cluster)
  for (i in 1:n_bootstrap){
    test_idx2 <- sample(test_idx, replace = T)
    ref_idx2 <- sample(ref_idx, replace = T)
    dist_list[[i]] <- get_log_or_dist(x, clusters, test_idx2, ref_idx2, pseudo_count)
  }
  var_list <- list()
  for (i in 1:ncol(dist_list[[1]])){
    #print(i)
    var_list[[i]] <- rowVars(do.call(cbind, lapply(dist_list, function(x, i){x[,i]}, i = i)))
  }
  var_mat <- do.call(cbind, var_list)
}

#smooths copy number estimates by combining information from adjacent bins based both on their distance (in base-pairs) from a bin and the deviation of the log-odds ratio from 1 (normalised by its expected standard deviation)
get_segmented_cn2 <- function(x, clusters, test_cluster, ref_cluster, is_excluded, pseudo_count, n_bootstrap, z_threshold, min_dist_weight = 0, sample_id = NULL){
  x <- x[,!is_excluded]
  clusters <- clusters[!is_excluded]
  sample_id <- sample_id[!is_excluded]
  if (!is.null(sample_id)){
    sample_sfs <- mean(colSums(x[,clusters %in% ref_cluster])) / sapply(split(colSums(x[,clusters %in% ref_cluster]), sample_id[clusters %in% ref_cluster]), mean)[sample_id]
    x <- t(t(x) * sample_sfs)
  }
  x <- x[,]
  chr_names <- sapply(strsplit(rownames(x), "\\."), `[`, 1)
  unique_chr <- sort(unique(chr_names))
  #x_list <- list()
  #for (i in unique_chr){
  #  x_list[[i]] <- x[chr_names == i,]
  #}
  #create full distance matrix and full bootstrapped values - find neigbours later
  log_or_dist <- get_log_or_dist(x, clusters, clusters %in% test_cluster, clusters %in% ref_cluster, pseudo_count)
  log_or_dist_bootstrap_vars <- bootstrap_log_or_dist(x, clusters, test_cluster, ref_cluster, pseudo_count, n_bootstrap = 100)
  log_or_z_initial <- log_or_dist / sqrt(log_or_dist_bootstrap_vars)

  #find neighbours
  is_same_chr <- chr_names[2:length(chr_names)] == chr_names[1:(length(chr_names) - 1)]
  neighbour_idx <- ncol(log_or_z_initial) * (1:ncol(log_or_z_initial) - 1) + 1:ncol(log_or_z_initial) + 1
  neighbour_zs <- log_or_z_initial[neighbour_idx][is_same_chr]

  #estimate standard deviation from neighbours
  #we don't bother centring - we assume each bin is an unbiased estimate of its neigbour on average
  robust_sigma <- 1.4826 * median(abs(neighbour_zs), na.rm = T)

  #re-estimate z
  log_or_z_final <- log_or_z_initial / robust_sigma

  get_log_or_weights <- function(x, threshold = 1){
    x[is.na(x)] <- 0
    x[is.infinite(x)] <- 0
    x[abs(x) <= threshold] <- 1
    x[abs(x) > threshold] <- (threshold / abs(x[abs(x) > threshold])) ** 2
    x
  }
  log_or_weights <- get_log_or_weights(log_or_z_final, threshold = z_threshold)

  get_distance_weights <- function(chr_names, min_weight = 1 / length(chr_names)){
    distance_weights <- outer(chr_names, chr_names, "==")
    distance_weights[!distance_weights] <- 0
    distance_weights[distance_weights] <- 1
    for (i in 1:ncol(distance_weights)){
      #linearly decreasing denominator
      #distance_weights[,i] <- distance_weights[,i] / (1 + abs(i - 1:nrow(distance_weights)))
      #exponentially decreasing denominator
      distance_weights[,i] <- distance_weights[,i] / 2 ** (abs(i - 1:nrow(distance_weights)))
      #0 if not i
      #distance_weights[,i] <- ifelse(1:nrow(distance_weights) == i, 1, 0)
    }
    distance_weights[distance_weights < min_weight] <- min_weight
    distance_weights
  }
  distance_weights <- get_distance_weights(chr_names, min_weight = min_dist_weight)

  smoothed_estimates <- numeric(nrow(x))
  raw_estimates <- rowMeans(x[,clusters == test_cluster]) / rowMeans(x[,clusters %in% ref_cluster])
  for (i in seq_along(smoothed_estimates)){
    combined_weights <- distance_weights[,i] * log_or_weights[,i]
    sf <- ifelse(grepl("chrX|chrY", rownames(x)), 1, 2)
    smoothed_estimates[i] <- mean(sf * x[,clusters == test_cluster] * combined_weights)  / mean(combined_weights * x[,clusters %in% ref_cluster])
  }
  names(smoothed_estimates) <- rownames(x)
  new_feature_names <- get_new_feature_names(rownames(x))
  smoothed_estimates[new_feature_names]
}

#computes absolute copy number estimates after computing scale factors - allows smoothing
#' @export
get_absolute_copy_number2 <- function(x, clusters, ref_cluster, scale_factors, is_excluded, is_ref_female, sample_id = NULL, do_smooth = T, pseudo_count = 1, n_bootstrap = 1000, z_threshold = 1, min_dist_weight = 0){
  new_feature_names <- get_new_feature_names(rownames(x))
  x <- x[,!is_excluded]
  #hacky way to add support for multiple reference clusters
  if (length(ref_cluster) > 1){
    clusters[clusters %in% ref_cluster & !(clusters == ref_cluster[1])] <- ref_cluster[1]
    ref_cluster <- ref_cluster[1]
  }
  clusters <- factor(clusters[!is_excluded])
  sample_id <- sample_id[!is_excluded]
  if (!is.null(sample_id)){
    sample_sfs <- mean(colSums(x[,clusters %in% ref_cluster])) / sapply(split(colSums(x[,clusters %in% ref_cluster]), sample_id[clusters %in% ref_cluster]), mean)[sample_id]
    x <- t(t(x) * sample_sfs)
  }
  cn.list <- list()
  for (i in levels(clusters)[levels(clusters) != ref_cluster]){
    print(paste0("Computing cluster ", i, " copy number"))
    sf <- scale_factors[paste0("cluster_", i)]
    if (!do_smooth){
      cn <- sf * rowMeans(x[new_feature_names, clusters == i]) /
        rowMeans(x[new_feature_names, clusters == ref_cluster])
      chr <- gsub("chr", "", sapply(strsplit(new_feature_names, "\\."), `[`, 1))
      chr <- factor(chr, unique(chr))
      cn <- ifelse(grepl("X|Y", chr) & !is_ref_female, 1, 2) * cn
    } else {
      cn <- scale_factors[paste0("cluster_", i)] * get_segmented_cn2(x, clusters, i, ref_cluster, rep(F, length(clusters)), pseudo_count, n_bootstrap, z_threshold, min_dist_weight, sample_id)
    }
    cn.list[[i]] <- cn
  }
  if (!do_smooth){
    cn <- colMeans(scale_factors[paste0("cluster_", clusters[clusters != ref_cluster])] * t(x[new_feature_names, clusters != ref_cluster])) /
      rowMeans(x[new_feature_names, clusters == ref_cluster])
    chr <- gsub("chr", "", sapply(strsplit(new_feature_names, "\\."), `[`, 1))
    chr <- factor(chr, unique(chr))
    cn <- ifelse(grepl("X|Y", chr) & !is_ref_female, 1, 2) * cn
    cn.list[["all"]] <- cn
  } else {
    clusters <- factor(clusters[!(clusters %in% ref_cluster)])
    cn.list[["all"]] <- sapply(data.table::transpose(cn.list), weighted.mean, as.integer(table(clusters)))
  }
  cn.list
}

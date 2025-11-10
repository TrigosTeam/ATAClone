get_copy_number <- function(count_matrix, clusters, ref.cluster, joint, ploidy.correction.method = "prop", external.ref = NULL){
  require("MatrixGenerics")
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
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- 0.5 * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  } else {
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- 2 * cn.mat[,as.integer(clusters) == i] / as.numeric(cor.list[[ref.cluster]])
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- 0.5 * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  }
  cn.mat[is.infinite(cn.mat)] <- NA
  cn.mat
}

get_new_order <- function(cn.mat, count_matrix_cor, clusters){
  #sort within clusteres based on corrected count matrix
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
  cluster.order <- c(ref.cluster, cluster.order[cluster.order != ref.cluster])
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
get_new_hclust <- function(){
  test_hclust <- hclust(dist(t(joint_cn_estimates[,new_cell_order]), method = "manhattan"), method = "single")
  test_clusters <- cutree(test_hclust, length(levels(iterative.leiden.clusters)))
  idx_split <- split(1:length(iterative.leiden.clusters), test_clusters[new_cell_order])

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

  merge_list <- mapply(get_cluster_merges, idx_split, 1:length(levels(iterative.leiden.clusters)), SIMPLIFY = F)
  #can't remember what this is
  #subgraph_indices <- sapply(merge_list, function(x){x[nrow(x),2] + 1})
  #make graphs of pseudo-bulk
  use.clusters <- 1:length(levels(iterative.leiden.clusters))
  cluster.cn <- list()
  for (i in seq_along(use.clusters)){
    cluster.cn[[i]] <- joint_cn_estimates[,iterative.leiden.clusters == use.clusters[i]][,1]
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

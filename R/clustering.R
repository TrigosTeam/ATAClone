normalise_counts <- function(x, overdispersion){
  transformGamPoi::acosh_transform(x, overdispersion, F)
}

correct_normalised_counts <- function(x, discard_pcs){
  x.pca <- prcomp(t(x))
  use.pcs <- 1:nrow(x)
  use.pcs <- use.pcs[!use.pcs %in% discard_pcs]
  x.cor <- t(x.pca$x[,use.pcs] %*% t(x.pca$rotation[,use.pcs])) + x.pca$center
  x.cor
}

get_pca <- function(x.norm, npcs){
  require("Matrix")
  #prcomp(t(x.norm))
  x.norm.svd <- irlba::irlba(scale(t(x.norm), T, F), nv = npcs)
  diagonal_matrix <- Matrix::Diagonal(npcs, x.norm.svd$d)
  x.norm.pca <- as.matrix(x.norm.svd$u %*% diagonal_matrix)
  x.norm.rotation <- as.matrix(x.norm.svd$v)
  rownames(x.norm.pca) <- colnames(x.norm)
  colnames(x.norm.pca) <- paste0("PC", 1:npcs)
  rownames(x.norm.rotation) <- rownames(x.norm)
  colnames(x.norm.pca) <- paste0("PC", 1:npcs)
  list(x = x.norm.pca, rotation = x.norm.rotation)
}

reorder_clusters <- function(clusters){
  sorted_tabs <- sort(table(clusters), decreasing = T)
  new_cluster_ids <- setNames(1:length(sorted_tabs), names(sorted_tabs))
  #levels shouldn't be here??
  as.integer(unname(new_cluster_ids[as.character(clusters)]), levels = 1:length(sorted_tabs))
}

scan_resolution <- function(knn.graph, start, stop, step, seed = NULL, n.iter = 20, tolerance = 0.01){
  if (!is.null(seed)){
    set.seed(seed)
  }
  for (i in seq(start, stop, step)){
    leiden.resolution <- i / (gorder(knn.graph) - 1)
    #leiden.clusters <- cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership
    #if (sum(leiden.clusters == 2) > 0){
    #  return(i - step)
    #}
    leiden.clusters <- unlist(lapply(1:n.iter, function(j){reorder_clusters(cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership)}))
    if (sum(leiden.clusters != 1) > gorder(knn.graph) * n.iter * tolerance){
      return(i - step)
    }
  }
  i
}

iterative_cluster_sim <- function(x, overdispersion, npcs, discard_pcs, k, start, stop, step, seed, iter.limit = 2, rename.clusters = T){
  set.seed(seed)
  x.norm.pca <- get_pca(normalise_counts(x, overdispersion), npcs)
  x.sim.norm.pca <- get_pca(normalise_counts(simulate_counts(x, overdispersion), overdispersion), npcs)
  leiden.clusters <- rep("start", ncol(x))
  use_pcs <- 1:npcs
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  use_pcs_sim <- 1:npcs
  discard_pcs_sim <- ifelse(1 %in% discard_pcs, 1, NULL)
  use_pcs_sim <- use_pcs_sim[!use_pcs_sim %in% discard_pcs_sim][1:length(use_pcs)]
  for (i in 1:iter.limit){
    leiden.clusters2 <- leiden.clusters
    for (j in unique(leiden.clusters)){
      knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j, use_pcs]), k = k, directed = F, d = NA)
      knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,use_pcs_sim]), k = k, directed = F, d = NA)
      leiden.resolution <- scan_resolution(knn.graph.sim, 0, k, 1, seed) / (gorder(knn.graph) - 1)
      leiden.clusters2[leiden.clusters == j] <- paste0(leiden.clusters[leiden.clusters == j], "_", reorder_clusters(cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership))
    }
    leiden.clusters <- leiden.clusters2
    if (length(unique(leiden.clusters)) > 1 & i != iter.limit){
      x.list <- list()
      for (j in unique(leiden.clusters)){
        x.list[[j]] <- x[,leiden.clusters == j]
      }
      x.sim.norm.pca <- get_pca(normalise_counts(do.call(cbind, lapply(x.list, simulate_counts, overdispersion))[rownames(x),colnames(x)], overdispersion), npcs)
    } else {
      break
    }
  }
  leiden.clusters <- gsub("start_", "", leiden.clusters)
  if (rename.clusters){
    sorted_tabs <- sort(table(leiden.clusters), decreasing = T)
    new_cluster_ids <- setNames(1:length(sorted_tabs), names(sorted_tabs))
    leiden.clusters <- factor(unname(new_cluster_ids[as.character(leiden.clusters)]), levels = 1:length(sorted_tabs))
    names(leiden.clusters) <- colnames(peak_counts_filtered)
  }
  leiden.clusters
}

get_cluster_outliers <- function(knn.graph, clusters, min_neighbour_quantile){
  outlier.list <- list()
  for (i in seq_along(levels(clusters))){
    cluster.id <- levels(clusters)[i]
    knn.subgraph <- subgraph(knn.graph, which(clusters == cluster.id))
    neighbour.prop <- degree(knn.subgraph) / degree(knn.graph)[which(clusters == cluster.id)]
    outlier.list[[cluster.id]] <- setNames(neighbour.prop < quantile(neighbour.prop, min_neighbour_quantile), names(clusters)[clusters == cluster.id])
  }
  is.outlier <- unlist(outlier.list)[paste0(clusters, ".", names(clusters))]
  names(is.outlier) <- sapply(strsplit(names(is.outlier), "\\."), `[`, 2)
  is.outlier
}

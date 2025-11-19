#' @export
normalise_counts <- function(x, overdispersion){
  transformGamPoi::acosh_transform(x, overdispersion, F)
}

#' @export
correct_normalised_counts <- function(x, discard_pcs){
  x.pca <- prcomp(t(x))
  use.pcs <- 1:nrow(x)
  use.pcs <- use.pcs[!use.pcs %in% discard_pcs]
  x.cor <- t(x.pca$x[,use.pcs] %*% t(x.pca$rotation[,use.pcs])) + x.pca$center
  x.cor
}

#' @export
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
    leiden.resolution <- i / (igraph::gorder(knn.graph) - 1)
    #leiden.clusters <- cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership
    #if (sum(leiden.clusters == 2) > 0){
    #  return(i - step)
    #}
    leiden.clusters <- unlist(lapply(1:n.iter, function(j){reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership)}))
    if (sum(leiden.clusters != 1) > igraph::gorder(knn.graph) * n.iter * tolerance){
      return(i - step)
    }
  }
  i
}

#' @export
iterative_cluster_sim <- function(x, overdispersion, npcs, discard_pcs, k, start = 0, stop = 2 * k, step = 1, seed = 100, iter.limit = 2, rename.clusters = T){
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
      leiden.resolution <- scan_resolution(knn.graph.sim, 0, k, 1, seed) / (igraph::gorder(knn.graph) - 1)
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
    names(leiden.clusters) <- colnames(x)
  }
  leiden.clusters
}

#' @export
get_umap <- function(pca_obj, npcs, discard_pcs){
  use_pcs <- 1:npcs
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  umap_embeddings <- uwot::umap(pca_obj[,use_pcs])
  colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")
  umap_embeddings
}

#' @export
get_cluster_outliers <- function (pca_obj, clusters, npcs, discard_pcs, knn_k, min_neighbour_quantile){
  use_pcs <- 1:npcs
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  knn.graph <- scran::buildKNNGraph(t(pca_obj[,use_pcs]), k = knn_k, directed = F, d = NA)
  outlier.list <- list()
  for (i in seq_along(levels(clusters))) {
    cluster.id <- levels(clusters)[i]
    knn.subgraph <- igraph::subgraph(knn.graph, which(clusters ==
                                                        cluster.id))
    neighbour.prop <- igraph::degree(knn.subgraph)/igraph::degree(knn.graph)[which(clusters ==
                                                                                     cluster.id)]
    outlier.list[[cluster.id]] <- setNames(neighbour.prop <
                                             quantile(neighbour.prop, min_neighbour_quantile),
                                           names(clusters)[clusters == cluster.id])
  }
  is.outlier <- unlist(outlier.list)[paste0(clusters, ".",
                                            names(clusters))]
  names(is.outlier) <- sapply(strsplit(names(is.outlier), "\\."),
                              `[`, 2)
  is.outlier
}


#' @export
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

get_internal_edges <- function(knn_graph, leiden_clusters){
  cluster.ids <- sort(unique(leiden_clusters))
  subgraph.list <- list()
  for (i in cluster.ids){
    subgraph.list[[i]] <- subgraph(knn_graph, which(leiden_clusters == i))
  }
  sum(sapply(subgraph.list, ecount))
}

#this was meant to provide a visual aid to picking best Leiden partition
#deprecated - changed to simulation approach
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

#this was meant to pick a best Leiden partition from many runs
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

#switched to pre-computed barcode probability
get_barcode_probability <- function(barcodes, weights = rep(1, length(barcodes))){
  barcode_1_9 <- substr(colnames(peak_counts), 1, 9)
  barcode_1_9_tabs <- table(rep(barcode_1_9, weights))
  barcode_1_9_probs <- (barcode_1_9_tabs / sum(barcode_1_9_tabs))[barcode_1_9]
  barcode_8_16 <- substr(colnames(peak_counts), 8, 16)
  barcode_8_16_tabs <- table(rep(barcode_8_16, weights))
  barcode_8_16_probs <- (barcode_8_16_tabs / sum(barcode_8_16_tabs))[barcode_8_16]
  setNames(barcode_1_9_probs * barcode_8_16_probs, barcodes)
}

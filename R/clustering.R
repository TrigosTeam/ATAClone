#' @export
find_lambda <- function(x, overdispersion, pseudo_count){
  mu_mat <- do.call(cbind, lapply(1:ncol(x), function(x2){rowMeans(x)}))
  sim_counts <- t(simulate_counts(mu_mat, overdispersion))
  #don't want to fit on bins with all zeros
  sim_counts2 <- sim_counts[,colMeans(sim_counts) > 0] + pseudo_count
  lm_out <- lm(sim_counts2~1, y = T)
  fit_box_cox <- boxcox(lm_out, plotit = F)
  lambda <- fit_box_cox$x[order(fit_box_cox$y, decreasing = T)][1]
  lambda
}

#' @export
normalise_counts <- function(x, overdispersion, pseudo_count = 0, lambda = NULL){
  if (is.null(lambda)){
    if (overdispersion == 0){
      2 * sqrt(x + pseudo_count)
    } else {
      #equivalent to acosh transformation in transformGamPoi when pseudo_count = 0
      2 * sqrt(1 / overdispersion - 0.5) * asinh(sqrt((x + pseudo_count) / (1 / overdispersion - 2 * pseudo_count)))
    }
  } else {
    #power transformation
    (x + pseudo_count) ** lambda
  }
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

scan_resolution2 <- function(knn.graph, start, stop, step, seed = NULL, n.iter = ceiling(500000 / gorder(knn.graph)), tolerance = 0.01){
  print("Beginning binary search")
  start_p <- 0
  stop_p <- 1
  if (!is.null(seed)){
    set.seed(seed)
  }
  leiden.resolution <- (start + stop) / 2
  ps <- numeric(100)
  resolutions <- numeric(100)
  #count number of times "start" of binary search range is 0 - if we get stuck alternating between 0 and a non-zero value, we probably won't reach
  #required tolerance level - so cancel search when we hit max_zeros and return largest "start" value seen so far
  zero_count <- 0
  max_zeros <- 5
  for (i in 1:100){
    #print(leiden.resolution)
    leiden.clusters <- unlist(lapply(1:n.iter, function(j){ATAClone:::reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution / (igraph::gorder(knn.graph) - 1))$membership)}))
    p <- (sum(leiden.clusters != 1)) / (igraph::gorder(knn.graph) * n.iter)
    #print(p)
    print(paste0("Tried resolution ", leiden.resolution, " (p=", round(p,4), ")"))
    resolutions[i] <- leiden.resolution
    ps[i] <- p
    if (p > stop_p | p < start_p){
      print("unstable!")
      #plot(resolutions[1:i], ps[1:i])
      #abline(h = tolerance)
      return(max(resolutions[1:i][ps[1:i] < tolerance]))
    }
    if (p > tolerance){
      stop <- leiden.resolution
      stop_p <- p
      new.leiden.resolution <- (leiden.resolution + start) / 2
      #new.leiden.resolution <- exp((log(leiden.resolution) + log(start)) / 2)
    } else {
      start <- leiden.resolution
      start_p <- p
      if (start_p == 0){
        zero_count <- zero_count + 1
        if (zero_count == max_zeros){
          return(leiden.resolution)
        }
      } else {
        zero_count <- 0
      }
      new.leiden.resolution <- (leiden.resolution + stop) / 2
      #new.leiden.resolution <- exp((log(leiden.resolution) + log(stop)) / 2)
    }
    #if ((stop - start) < 1){
    #  return(leiden.resolution)
    #}
    if (abs(tolerance - p) < (tolerance / 10)){
      #plot(resolutions[1:i], ps[1:i])
      #abline(h = tolerance)
      return(leiden.resolution)
    } else{
      leiden.resolution <- new.leiden.resolution
    }
  }
  print("binary search did not finish")
  leiden.resolution
}

scan_resolution <- function(knn.graph, start, stop, step, seed = NULL, n.iter = ceiling(500000 / gorder(knn.graph)), tolerance = 0.01){
  if (!is.null(seed)){
    set.seed(seed)
  }
  ps <- numeric(length(seq(start, stop, step)))
  resolutions <- numeric(length(seq(start, stop, step)))
  for (j in seq_along(seq(start, stop, step))){
    i <- seq(start, stop, step)[j]
    leiden.resolution <- i / (igraph::gorder(knn.graph) - 1)
    #leiden.clusters <- cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership
    #if (sum(leiden.clusters == 2) > 0){
    #  return(i - step)
    #}
    leiden.clusters <- unlist(lapply(1:n.iter, function(j){ATAClone:::reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership)}))
    p <- (sum(leiden.clusters != 1)) / (igraph::gorder(knn.graph) * n.iter)
    ps[j] <- p
    resolutions[j] <- i
    #print(i)
    #print(p)
    if (sum(leiden.clusters != 1) > igraph::gorder(knn.graph) * n.iter * tolerance){
      #plot(resolutions[1:j], ps[1:j])
      #abline(h = tolerance)
      return(i - step)
    }
  }
  #plot(resolutions[1:j], ps[1:j])
  #abline(h = tolerance)
  i
}

#' @export
iterative_cluster_sim <- function(x, overdispersion, pseudo_count, lambda = NULL, npcs, discard_pcs, k, start = 0, stop = 2 * k, step = 1, seed = 100, iter.limit = 2, rename.clusters = T, tolerance = 0.01){
  set.seed(seed)
  x.norm.pca <- get_pca(normalise_counts(x, overdispersion, pseudo_count, lambda), npcs)
  x.sim.norm.pca <- get_pca(normalise_counts(simulate_counts(x, 0.015), overdispersion, pseudo_count, lambda), npcs)
  leiden.clusters <- rep("start", ncol(x))
  use_pcs <- 1:npcs
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  use_pcs_sim <- 1:npcs
  discard_pcs_sim <- ifelse(1 %in% discard_pcs, 1, NULL)
  use_pcs_sim <- use_pcs_sim[!use_pcs_sim %in% discard_pcs_sim][1:length(use_pcs)]
  for (i in 1:iter.limit){
    print(paste0("Beginning clustering iteration ", i))
    leiden.clusters2 <- leiden.clusters
    for (j in unique(leiden.clusters)){
      print(paste0("Finding resolution for simulated cluster ", j))
      knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j, use_pcs]), k = k, directed = F, d = NA)
      knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,use_pcs_sim]), k = k, directed = F, d = NA)
      leiden.resolution <- scan_resolution2(knn.graph.sim, start, stop, step, seed, tolerance = tolerance) / (igraph::gorder(knn.graph) - 1)
      leiden.clusters2[leiden.clusters == j] <- paste0(leiden.clusters[leiden.clusters == j], "_", ATAClone:::reorder_clusters(cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership))
    }
    leiden.clusters <- leiden.clusters2
    print(paste0("Found ", length(unique(leiden.clusters)), " clusters"))
    if (length(unique(leiden.clusters)) > 1 & i != iter.limit){
      x.list <- list()
      for (j in unique(leiden.clusters)){
        x.list[[j]] <- x[,leiden.clusters == j]
      }
      x.sim.norm.pca <- get_pca(normalise_counts(do.call(cbind, lapply(x.list, simulate_counts, 0.015))[rownames(x),colnames(x)], overdispersion, pseudo_count, lambda), npcs)
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

#new normalisation function for two-step normalisation - input into PCA on a log-scale, project without PC1, then exponentiate and re-input to PCA to stabilise variance
#while removing library size effect from both "stable" and highly-variable genes (projecting just variance-stabilised data without PC1 only removes library size from stable features)
#' @export
normalise_counts2 <- function(x, overdispersion, pseudo_count, lambda){
  log((x + pseudo_count) ** lambda)
}

#simulation-based clustering using the above-described two-step normalisation strategy
#' @export
iterative_cluster_sim2 <- function(x, overdispersion, pseudo_count, lambda = NULL, npcs, discard_pcs, k, start = 0, stop = 2 * k, step = 1, seed = 100, iter.limit = 2, rename.clusters = T, tolerance = 0.01){
  set.seed(seed)
  use_pcs <- 1:npcs
  #use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  use_pcs_sim <- 1:npcs
  #discard_pcs_sim <- ifelse(1 %in% discard_pcs, 1, NULL)
  discard_pcs_sim <- 1
  use_pcs_sim <- use_pcs_sim[!use_pcs_sim %in% discard_pcs_sim][1:length(use_pcs)]
  #x.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(x, overdispersion, pseudo_count, lambda), discard_pcs)), npcs)
  x.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(x, overdispersion, pseudo_count, lambda), 1)), npcs)
  x.sim.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(simulate_counts(x, overdispersion), overdispersion, pseudo_count, lambda), discard_pcs_sim)), npcs)
  leiden.clusters <- rep("start", ncol(x))
  for (i in 1:iter.limit){
    print(paste0("Beginning clustering iteration ", i))
    leiden.clusters2 <- leiden.clusters
    for (j in unique(leiden.clusters)){
      print(paste0("Finding resolution for simulated cluster ", j))
      #knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j,use_pcs]), k = k, directed = F, d = NA)
      #knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,use_pcs_sim]), k = k, directed = F, d = NA)
      leiden.resolution <- ATAClone:::scan_resolution2(knn.graph.sim, start, stop, step, seed, tolerance = tolerance) / (igraph::gorder(knn.graph) - 1)
      leiden.clusters2[leiden.clusters == j] <- paste0(leiden.clusters[leiden.clusters == j], "_", ATAClone:::reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership))
    }
    leiden.clusters <- leiden.clusters2
    print(paste0("Found ", length(unique(leiden.clusters)), " clusters"))
    if (length(unique(leiden.clusters)) > 1 & i != iter.limit){
      x.list <- list()
      for (j in unique(leiden.clusters)){
        x.list[[j]] <- x[,leiden.clusters == j]
      }
      x.sim.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(do.call(cbind, lapply(x.list, simulate_counts, overdispersion))[rownames(x),colnames(x)], overdispersion, pseudo_count, lambda), discard_pcs_sim)), npcs)
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
regress_bins <- function(x, regress_var){
  x.t <- t(x)
  fit.lm <- lm(x.t ~ regress_var)
  t(x.t - fit.lm$fitted.values) + rowMeans(x)
}

#testing allowing regression of unwanted factors
#' @export
iterative_cluster_sim3 <- function(x, overdispersion, pseudo_count, lambda = NULL, npcs, discard_pcs, k, start = 0, stop = 2 * k, step = 1, seed = 100, iter.limit = 2, rename.clusters = T, tolerance = 0.01, regress.var = NULL){
  set.seed(seed)
  use_pcs <- 1:npcs
  #use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  use_pcs <- use_pcs[!use_pcs %in% discard_pcs]
  #need to use more PCs if you don't discard any from real data - we will always discard PC1 from the simulated data
  npcs_sim <- npcs + 1 - length(discard_pcs)
  use_pcs_sim <- 1:npcs_sim
  #discard_pcs_sim <- ifelse(1 %in% discard_pcs, 1, NULL)
  discard_pcs_sim <- 1
  use_pcs_sim <- use_pcs_sim[!use_pcs_sim %in% discard_pcs_sim][1:(length(use_pcs))]
  #x.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(x, overdispersion, pseudo_count, lambda), discard_pcs)), npcs)
  x.cor <- exp(correct_normalised_counts(normalise_counts2(x, overdispersion, pseudo_count, lambda), 1))
  if (!is.null(regress.var)){
    x.cor <- regress_bins(x.cor, regress.var)
  }
  x.norm.pca <- get_pca(x.cor, npcs)
  x.sim.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(simulate_counts(x, overdispersion), overdispersion, pseudo_count, lambda), discard_pcs_sim)), npcs_sim)
  leiden.clusters <- rep("start", ncol(x))
  for (i in 1:iter.limit){
    print(paste0("Beginning clustering iteration ", i))
    leiden.clusters2 <- leiden.clusters
    for (j in unique(leiden.clusters)){
      print(paste0("Finding resolution for simulated cluster ", j))
      #knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j,use_pcs]), k = k, directed = F, d = NA)
      #knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,use_pcs_sim]), k = k, directed = F, d = NA)
      leiden.resolution <- ATAClone:::scan_resolution2(knn.graph.sim, start, stop, step, seed, tolerance = tolerance) / (igraph::gorder(knn.graph) - 1)
      leiden.clusters2[leiden.clusters == j] <- paste0(leiden.clusters[leiden.clusters == j], "_", ATAClone:::reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership))
    }
    leiden.clusters <- leiden.clusters2
    print(paste0("Found ", length(unique(leiden.clusters)), " clusters"))
    if (length(unique(leiden.clusters)) > 1 & i != iter.limit){
      x.list <- list()
      for (j in unique(leiden.clusters)){
        x.list[[j]] <- x[,leiden.clusters == j]
      }
      x.sim.norm.pca <- get_pca(exp(correct_normalised_counts(normalise_counts2(do.call(cbind, lapply(x.list, simulate_counts, overdispersion))[rownames(x),colnames(x)], overdispersion, pseudo_count, lambda), discard_pcs_sim)), npcs_sim)
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

estimate_lib_size_pca <- function(x){
  x_pca <- prcomp(sqrt(t(x)))
  pc1_est <- t(x_pca$x[,1] %*% t(x_pca$rotation[,1])) + rowMeans(sqrt(x))
  colSums(pc1_est ** 2)
}

simulate_counts2 <- function(mu.mat, overdispersion){
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

#new, simpler glm approach - more robust PCA approaches
#' @export
iterative_cluster_sim_glm <- function(x, overdispersion = NULL, regress.var = NULL, pseudo_count = 0.25, lambda = 0.4, npcs, k, start = 0, stop = 2 * k, step = 1, seed = 100, iter.limit = 2, rename.clusters = T, tolerance = 0.01){
  set.seed(seed)
  lib_size <- estimate_lib_size_pca(x)
  nb.glm <- glmGamPoi::glm_gp(data = x, design = ~ log(regress.var), size_factors = lib_size, overdispersion = overdispersion)
  x.cor <- exp(log((x + pseudo_count) ** lambda) - log((nb.glm$Mu + pseudo_count) ** lambda) + log(rowMeans((x + pseudo_count) ** lambda)))
  #could simulate directly from mu?
  x.sim <- simulate_counts2(nb.glm$Mu, overdispersion)
  x.sim.cor <- exp(log((x.sim + pseudo_count) ** lambda) - log((nb.glm$Mu + pseudo_count) ** lambda) + log(rowMeans((x.sim + pseudo_count) ** lambda)))
  x.norm.pca <- get_pca(x.cor, npcs)
  x.sim.norm.pca <- get_pca(x.sim.cor, npcs)
  leiden.clusters <- rep("start", ncol(x))
  for (i in 1:iter.limit){
    print(paste0("Beginning clustering iteration ", i))
    leiden.clusters2 <- leiden.clusters
    for (j in unique(leiden.clusters)){
      print(paste0("Finding resolution for simulated cluster ", j))
      knn.graph <- scran::buildKNNGraph(t(x.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      knn.graph.sim <- scran::buildKNNGraph(t(x.sim.norm.pca$x[leiden.clusters == j,]), k = k, directed = F, d = NA)
      leiden.resolution <- ATAClone:::scan_resolution2(knn.graph.sim, start, stop, step, seed, tolerance = tolerance) / (igraph::gorder(knn.graph) - 1)
      leiden.clusters2[leiden.clusters == j] <- paste0(leiden.clusters[leiden.clusters == j], "_", ATAClone:::reorder_clusters(igraph::cluster_leiden(knn.graph, "CPM", resolution = leiden.resolution)$membership))
    }
    leiden.clusters <- leiden.clusters2
    print(paste0("Found ", length(unique(leiden.clusters)), " clusters"))
    if (length(unique(leiden.clusters)) > 1 & i != iter.limit){
      #resimulate after re-estimating means
      #testing allowing interaction (switch * to + to make additive)
      #col_data <- data.frame(cluster = factor(leiden.clusters), beta = log(regress.var))
      #design_mat <- model.matrix(~ cluster * beta, data = col_data)
      col_data <- data.frame(cluster = factor(leiden.clusters))
      design_mat <- model.matrix(~ cluster, data = col_data)
      design_mat <- cbind(design_mat, beta = log(regress.var))
      nb.glm2 <- glmGamPoi::glm_gp(data = x, design = design_mat, size_factors = lib_size, overdispersion = overdispersion)
      x.sim <- simulate_counts2(nb.glm2$Mu, overdispersion)
      x.sim.cor <- exp(log((x.sim + pseudo_count) ** lambda) - log((nb.glm$Mu + pseudo_count) ** lambda) + log(rowMeans((x.sim + pseudo_count) ** lambda)))
      x.sim.norm.pca <- get_pca(x.sim.cor, npcs)
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

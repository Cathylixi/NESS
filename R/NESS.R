#' NESS: Neighbor Embedding Stability Scoring
#'
#' Performs dimensionality reduction (t-SNE, UMAP, or PHATE) multiple times to evaluate
#' local neighbor stability across repeated embeddings. This function helps assess
#' the robustness of low-dimensional embeddings for high-dimensional data.
#'
#' @param data A numeric matrix or data frame with rows as observations and columns as features.
#' @param ... Additional arguments passed to the dimensionality reduction methods
#'            (`Rtsne`, `uwot::umap`, or `phateR::phate`), such as `theta` for t-SNE,
#'            `min_dist` for UMAP, or `decay` for PHATE.
#' @param data.name Character string used in plot titles to label the dataset.
#' @param GCP Optional numeric vector of neighborhood sizes (e.g., perplexity for t-SNE or
#'            number of neighbors for UMAP/PHATE). If `NULL`, a default sequence is generated.
#' @param cluster Optional vector of cluster labels for coloring the embedding plots.
#' @param rareness Logical; if `TRUE`, computes rareness metrics based on neighbor consistency across embeddings.
#' @param method Dimensionality reduction method to use. One of `"tsne"`, `"umap"`, or `"phateR"`.
#' @param initialization Initialization method: `1` for random, `2` for PCA-based initialization.
#' @param stability_threshold Quantile threshold (default = 0.75) for determining local neighbor stability.
#' @param early_stop Logical; if `TRUE`, stops early if global stability saturates.
#' @param seed_base Base random seed used for reproducibility.
#' @param N Number of repeated embedding runs.
#' @param k Number of nearest neighbors to use when computing stability metrics (default = 50).
#' @param svd_cutoff_ratio Threshold ratio used to estimate intrinsic dimensionality via SVD (default = 1.1).
#' @param svd_max_k Maximum number of SVD components to check when estimating dimensionality (default = 30).
#' @param stop_global_stability_threshold Early stopping threshold for global stability (default = 0.9).
#' @param stop_relative_change_threshold Early stopping threshold for relative improvement in global stability (default = 0.05).
#'
#' @return A list containing:
#' \describe{
#'   \item{GCP}{Vector of neighborhood sizes used for evaluation.}
#'   \item{GCP.optim}{The selected GCP value corresponding to the median rareness score.}
#'   \item{rare.mean}{(optional) Vector of rareness mean scores across GCP values (if `rareness = TRUE`).}
#'   \item{rare.var}{(optional) Vector of rareness variance scores across GCP values (if `rareness = TRUE`).}
#'   \item{embedding}{Embedding coordinates for the optimal GCP value.}
#'   \item{local_stability}{Vector of local kNN stability scores (no names).}
#'   \item{global_stability}{Vector of global stability scores across GCP values.}
#'   \item{embedding_stability_colored}{A ggplot2 plot of the embedding colored by local stability score.}
#'   \item{global_stability_plot}{A ggplot2 line plot showing global stability across GCP values.}
#'   \item{embedding_cluster_colored}{(optional) Embedding plot colored by cluster labels, if `cluster` is provided.}
#'   \item{rareness_mean}{(optional) ggplot2 plot of rareness score mean, if `rareness = TRUE`.},
#'   \item{par}{A list of input parameters used to run the function for reproducibility.}
#' }
#'
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom phateR phate
#' @importFrom RSpectra svds
#' @importFrom BiocNeighbors findKNN
#' @importFrom stats median
#' @export


NESS <- function(
    data,
    ...,
    data.name = "",
    GCP = NULL,
    cluster = NULL,
    method = "tsne",
    initialization = 1,
    stability_threshold = 0.75,
    early_stop = FALSE,
    seed_base = 1,
    N = 20,
    k = 50,
    svd_cutoff_ratio = 1.1,
    svd_max_k = 30,
    stop_global_stability_threshold = 0.98,
    stop_relative_change_threshold = 0.02
) {
  # Setup parameters
  set.seed(seed_base)
  eg.v <- svds((data), k = svd_max_k)$d
  pc <- max(which(eg.v[1:(svd_max_k - 1)] / eg.v[2:svd_max_k] > svd_cutoff_ratio))
  svd.data <- svds(data, k = pc)
  data.denoise <- svd.data$v %*% diag(svd.data$d)

  if (is.null(GCP)) {
    GCP <- round(seq(10, 10 * log(nrow(data)), length.out = 5))
  }

  global_stability <- c()
  rare.median <- c()
  rare.var <- c()
  rare_violin <- list()
  best_embedding <- list()
  embedding_plot_list <- list()
  plot_list_stability <- list()
  local_stability <- list()
  prev_stability <- NULL

  # Running NESS main analysis
  for (j in seq_along(GCP)) {
    cat(paste0("Running GCP value: ", GCP[j], " (", j, "/", length(GCP), ")\n"))
    gcp <- GCP[j]
    set.seed(seed_base + j)
    knn.mat <- list()
    knn.info <- NULL
    all_embeddings <- list()

    for (i in seq_len(N)) {
      if (i %% 10 == 0 || i == N) cat(paste0("  ", toupper(method), " iteration: ", i, "/", N, "\n"))
      set.seed(seed_base + i)

      embedding <- switch(method,
                          tsne = Rtsne(data.denoise, perplexity = gcp, check_duplicates = FALSE, pca = (initialization == 2), seed = seed_base + i, ...)$Y,
                          umap = {
                            if (is.null(knn.info)) {
                              out <- uwot::umap(data.denoise, n_neighbors = gcp, init = ifelse(initialization == 2, "pca", "lvrandom"), n_threads = 1, ret_nn = TRUE, ...)
                              knn.info <- list(idx = out$nn$euclidean$idx, dist = out$nn$euclidean$dist)
                              out$embedding
                            } else {
                              uwot::umap(data.denoise, n_neighbors = gcp, init = ifelse(initialization == 2, "pca", "lvrandom"), n_threads = 1, nn_method = knn.info, ret_nn = FALSE, ...)
                            }
                          },
                          phateR = {
                            if (i == 1) phate.init <- NULL
                            out <- phate(data.denoise, ndim = 2, knn = gcp, seed = seed_base + i, verbose = TRUE, init = if (exists("phate.init")) phate.init else NULL, ...)
                            if (i == 1) phate.init <- out
                            out$embedding
                          },
                          stop("Unsupported method")
      )
      all_embeddings[[i]] <- embedding
      knn.mat[[i]] <- findKNN(embedding, k = k)$index
    }

    # Compute rareness score
    temp.out <- matrix(ncol = N, nrow = N)
    for (p in 1:(N - 1)) {
      for (q in (p + 1):N) {
        temp <- sapply(1:nrow(knn.mat[[1]]), function(mn) {
          length(intersect(knn.mat[[p]][mn, 1:k], knn.mat[[q]][mn, 1:k])) / k
        })
        temp.out[p, q] <- median(temp)
      }
    }

    rareness_mean <- sapply(1:N, function(aa) {
      if (aa == 1) {
        mean(temp.out[1, 2:N])
      } else if (aa == N) {
        mean(temp.out[1:(N - 1), N])
      } else {
        mean(c(temp.out[1:(aa - 1), aa], temp.out[aa, (aa + 1):N]))
      }
    })
    rareness_variance <- sapply(1:N, function(aa) {
      temp_list <- if (aa == 1) temp.out[1, 2:N]
      else if (aa == N) temp.out[1:(N - 1), N]
      else c(temp.out[1:(aa - 1), aa], temp.out[aa, (aa + 1):N])
      sum((temp_list - rareness_mean[aa])^2) / length(temp_list)
    })

    rare.median[j] <- median(temp.out[upper.tri(temp.out)], na.rm = TRUE)
    rare.var[j] <- var(temp.out[upper.tri(temp.out)], na.rm = TRUE)
    rare_violin[[j]] <- plot_rareness_violin(temp.out, gcp)

    # Find the best embedding
    best_embedding_index <- which.min(rareness_mean)
    best_embedding[[j]] <- all_embeddings[[best_embedding_index]]
    embedding_plot_list[[j]] <- plot_embedding(best_embedding[[j]], cluster, data.name, method, gcp, pc)

    # Calculate stability score
    knn.graph <- matrix(0, nrow = nrow(knn.mat[[1]]), ncol = nrow(knn.mat[[1]]))
    for (i in seq_len(N)) {
      for (l in seq_len(nrow(knn.mat[[1]]))) {
        knn.graph[l, knn.mat[[i]][l, 1:k]] <- knn.graph[l, knn.mat[[i]][l, 1:k]] + 1
      }
    }
    knn.graph <- pmax(knn.graph, t(knn.graph))

    knn.score <- sapply(seq_len(nrow(knn.mat[[1]])), function(i) {
      quantile(knn.graph[i, knn.graph[i, ] != 0] / N, stability_threshold)
    })
    names(knn.score) <- NULL
    local_stability[[j]] <- knn.score
    plot_list_stability[[j]] <- plot_embedding_stability(best_embedding[[j]], knn.score, data.name, method, gcp)

    global_stability[j] <- mean(knn.score, na.rm = TRUE)
    print(paste("GCP =", GCP[j], "→ Global Stability =", round(global_stability[j], 4)))

    stability_plot <- plot_line(GCP[1:length(global_stability)], global_stability, "Global Stability", paste0(data.name, method, " - Global Stability"))


    # Find the best GCP score
    if (any(global_stability > stop_global_stability_threshold)) {
      best_gcp_idx <- which(global_stability > stop_global_stability_threshold)[1]
      if (early_stop) break
    } else if (any(diff(global_stability) < stop_relative_change_threshold)) {
      best_gcp_idx <- which(diff(global_stability) < stop_relative_change_threshold)[1] + 1
      if (early_stop) break
    } else {
      best_gcp_idx <- length(global_stability)
    }
  }

  result_list <- list(
    GCP = GCP,
    GCP.optim = GCP[best_gcp_idx],
    embedding = best_embedding[[best_gcp_idx]],
    embedding_plot = embedding_plot_list[[best_gcp_idx]],
    local_stability = local_stability[[best_gcp_idx]],
    global_stability = global_stability,
    embedding_stability_colored = plot_list_stability[[best_gcp_idx]],
    global_stability_plot = stability_plot,
    rare.median = rare.median,
    rare.var = rare.var,
    rare_violin = rare_violin[[best_gcp_idx]],
    par = list(
      data.name = data.name,
      GCP = GCP,
      cluster = cluster,
      rareness = TRUE,
      method = method,
      initialization = initialization,
      stability_threshold = stability_threshold,
      early_stop = early_stop,
      seed_base = seed_base,
      N = N,
      k = k,
      svd_cutoff_ratio = svd_cutoff_ratio,
      svd_max_k = svd_max_k,
      stop_global_stability_threshold = stop_global_stability_threshold,
      stop_relative_change_threshold = stop_relative_change_threshold
    )
  )
  return(result_list)
}


# Plot helper for stability embedding
plot_embedding_stability <- function(embedding, knn_score, data.name, method, gcp) {
  plot.data <- data.frame(dim1 = embedding[, 1], dim2 = embedding[, 2], knn_score = knn_score)
  ggplot(plot.data, aes(dim1, dim2, colour = knn_score)) +
    geom_point(size = 1) +
    ggtitle(paste0(data.name, "_", method, "_p", gcp, "_colored by local stability score")) +
    theme_minimal()
}

# Plot helper for cluster-colored embedding
plot_embedding <- function(embedding, cluster, data.name, method, gcp, pc) {
  plot.data <- data.frame(dim1 = embedding[, 1], dim2 = embedding[, 2])

  if (!is.null(cluster)) {
    plot.data$cluster <- cluster
    ggplot(plot.data, aes(dim1, dim2, colour = cluster)) +
      geom_point(size = 1) +
      ggtitle(paste0(data.name, "_", method, "_p", gcp, "_pc", pc, "_colored by cluster")) +
      theme_minimal()
  } else {
    ggplot(plot.data, aes(dim1, dim2)) +
      geom_point(size = 1, color = "black") +
      ggtitle(paste0(data.name, "_", method, "_p", gcp, "_pc", pc)) +
      theme_minimal()
  }
}


# Generic plot helper for line plots
plot_line <- function(x, y, y_label, title) {
  df <- data.frame(x = x, y = y)
  ggplot(df, aes(x = x, y = y)) +
    geom_line(color = "blue") +
    geom_point(size = 2) +
    xlab("GCP") +
    ylab(y_label) +
    ggtitle(title) +
    theme_minimal()
}

# Rareness boxplot helper
plot_rareness_violin <- function(temp.out, gcp) {
  df <- data.frame(value = temp.out[upper.tri(temp.out)],
                   gcp = paste0("GCP=", gcp))

  ggplot(df, aes(x = gcp, y = value)) +
    geom_violin(fill = "steelblue", color = "black", alpha = 0.8) +
    labs(
      title = paste0("Violin Plot of Rareness Score", "_p", gcp),
      x = "",  # You could put "GCP" if you want to label the axis generally
      y = "Pairwise Rareness Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
}





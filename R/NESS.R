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
#' @export


NESS <- function(
    data,
    ...,
    data.name = "",
    GCP = NULL,
    cluster = NULL,
    rareness = FALSE,
    method = "tsne",
    initialization = 1,
    stability_threshold = 0.75,
    early_stop = TRUE,
    seed_base = 1,
    N = 20,
    k = 50,
    svd_cutoff_ratio = 1.1,
    svd_max_k = 30,
    stop_global_stability_threshold = 0.9,
    stop_relative_change_threshold = 0.05
) {
  set.seed(seed_base)
  eg.v <- svds((data), k = svd_max_k)$d
  pc <- max(which(eg.v[1:(svd_max_k - 1)] / eg.v[2:svd_max_k] > svd_cutoff_ratio))
  svd.data <- svds(data, k = pc)
  data.denoise <- svd.data$v %*% diag(svd.data$d)

  if (is.null(GCP)) {
    GCP <- round(seq(10, 10 * log(nrow(data)), length.out = 5))
  }

  global_stability <- c()
  rare.mean <- c()
  rare.var <- c()
  all_embedding <- list()
  plot_list_cluster <- list()
  plot_list_stability <- list()
  local_stability <- list()
  init_type <- ifelse(initialization == 2, "PCA Initialization", "Random Initialization")
  prev_stability <- NULL
  tracking_GCP_length <- 0

  for (j in seq_along(GCP)) {
    cat(paste0("Running GCP value: ", GCP[j], " (", j, "/", length(GCP), ")\n"))
    tracking_GCP_length <- j
    gcp <- GCP[j]
    set.seed(seed_base + j)
    knn.mat <- list()
    knn.info <- NULL

    if (method == "tsne") {
      for (i in seq_len(N)) {
        if (i %% 10 == 0 || i == N) {
          cat(paste0("  TSNE iteration: ", i, "/", N, "\n"))
        }
        out <- Rtsne(data.denoise,
                          perplexity = gcp,
                          check_duplicates = FALSE,
                          pca = (initialization == 2),
                          seed = seed_base + i,
                          ...)
        embedding <- out$Y
        knn.mat[[i]] <- findKNN(embedding, k = k)$index
      }
    }
   else if (method == "umap") {
      for (i in seq_len(N)) {
        if (i %% 10 == 0 || i == N) {
          cat(paste0("  UMAP iteration: ", i, "/", N, "\n"))
        }
        set.seed(seed_base + i)
        if (is.null(knn.info)) {
          out <- uwot::umap(data.denoise,
                            n_neighbors = gcp,
                            init = ifelse(initialization == 2, "pca", "lvrandom"),
                            n_threads = 1,
                            ret_nn = TRUE,
                            ...)
          knn.info <- list(idx = out$nn$euclidean$idx, dist = out$nn$euclidean$dist)
          embedding <- out$embedding
        } else {
          embedding <- uwot::umap(data.denoise,
                                  n_neighbors = gcp,
                                  init = ifelse(initialization == 2, "pca", "lvrandom"),
                                  n_threads = 1,
                                  nn_method = knn.info,
                                  ret_nn = FALSE,
                                  ...)
        }
        knn.mat[[i]] <- findKNN(embedding, k = k)$index
      }
    }
    else if (method == "phateR") {
      phate.init <- NULL
      for (i in seq_len(N)) {
        if (i %% 10 == 0 || i == N) {
          cat(paste0("  PHATE iteration: ", i, "/", N, "\n"))
        }
        set.seed(seed_base + i)
        out <- phate(data.denoise,
                     ndim = 2,
                     knn = gcp,
                     seed = seed_base + i,
                     verbose = TRUE,
                     init = phate.init,
                     ...)
        if (i == 1) phate.init <- out
        embedding <- out$embedding
        knn.mat[[i]] <- findKNN(embedding, k = k)$index
      }
    }

    Y <- embedding

    knn.graph <- matrix(0, nrow = nrow(knn.mat[[1]]), ncol = nrow(knn.mat[[1]]))
    for (i in seq_len(N)) {
      for (l in seq_len(nrow(knn.mat[[1]]))) {
        knn.graph[l, knn.mat[[i]][l, 1:k]] <- knn.graph[l, knn.mat[[i]][l, 1:k]] + 1
      }
      knn.graph <- pmax(knn.graph, t(knn.graph))
    }

    knn.score <- sapply(seq_len(nrow(knn.mat[[1]])), function(i) {
      quantile(knn.graph[i, knn.graph[i, ] != 0] / N, stability_threshold)
    })
    names(knn.score) <- NULL
    local_stability[[j]] <- knn.score

    all_embedding[[j]] <- Y

    if (rareness) {
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

      rare.mean[j] <- median(temp.out[upper.tri(temp.out)], na.rm = TRUE)
      rare.var[j] <- var(temp.out[upper.tri(temp.out)], na.rm = TRUE)
    }

    if (!is.null(cluster)) {
      plot.data <- data.frame(dim1 = Y[, 1], dim2 = Y[, 2], cluster = cluster)
      plot11 <- ggplot(plot.data, aes(dim1, dim2, colour = cluster)) +
        geom_point(size = 1) +
        ggtitle(paste0(data.name, "_", method, "(", init_type, ")" , "_p", gcp, "_pc", pc, "_colored by cluster"))
      plot_list_cluster[[j]] <- plot11
    }

    plot.data <- data.frame(dim1 = Y[, 1], dim2 = Y[, 2], knn_score = knn.score)
    plot22 <- ggplot(plot.data, aes(dim1, dim2, colour = knn_score)) +
      geom_point(size = 1) +
      ggtitle(paste0(data.name, "_", method, "(", init_type, ")" , "_p", gcp, "_pc", pc, "_colored by local stability score"))
    plot_list_stability[[j]] <- plot22

    global_stability[j] <- mean(knn.score, na.rm = TRUE)

    print(paste("GCP =", GCP[j], "→ Global Stability =", round(global_stability[j], 4)))
    if (!is.null(prev_stability)) {
      relative_change <- (global_stability[j] - prev_stability) / abs(prev_stability)
      prev_stability <- global_stability[j]


      if (early_stop && (global_stability[j] >= stop_global_stability_threshold || relative_change <= stop_relative_change_threshold)) {
        break
      }
    }
  }

  names(global_stability) <- GCP
  names(rare.mean) <- GCP
  names(rare.var) <- GCP
  rare.mean.vec <- unlist(rare.mean)[1:tracking_GCP_length]
  med_value <- median(rare.mean.vec, na.rm = TRUE)
  med_index <- which.min(abs(rare.mean.vec - med_value))
  print(paste("Index of MEDIAN rareness score:", med_index))
  print(paste("Corresponding GCP value for median rareness score:", GCP[med_index]))

  plot_improved <- function(x, y, y_label, title) {
    df <- data.frame(x = x, y = y)
    ggplot(df, aes(x = x, y = y)) +
      geom_line(color = "blue") +
      geom_point(size = 2) +
      xlab("GCP") +
      ylab(y_label) +
      ggtitle(title) +
      theme_minimal()
  }

  stability_plot <- plot_improved(GCP[1:tracking_GCP_length], global_stability, "Global Stability", paste0(data.name, method, " - Global Stability (", init_type, ")"))

  if (early_stop) {
    message <- paste("The final GCP is", GCP[med_index], "with median rareness score =", rare.mean.vec[med_index])
    print(message)
  }

  result_list <- list(
    GCP = GCP,
    GCP.optim = GCP[med_index],
    rare.mean = rare.mean,
    rare.var = rare.var,
    embedding = all_embedding[[med_index]],
    local_stability = local_stability[[med_index]],
    global_stability = global_stability,
    embedding_stability_colored = plot_list_stability[[med_index]],
    global_stability_plot = stability_plot,
    par = list(
      data.name = data.name,
      GCP = GCP,
      cluster = cluster,
      rareness = rareness,
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
  if (!is.null(cluster)) {
    result_list$embedding_cluster_colored <- plot_list_cluster[[med_index]]
  }
  if (rareness) {
    result_list$rareness_mean <- plot_improved(GCP[1:tracking_GCP_length], rare.mean, "Rareness Score(Mean)", paste0(data.name, method, " - Rareness Score(Mean) (", init_type, ")"))
  }

  return(result_list)
}

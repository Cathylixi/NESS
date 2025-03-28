#' NESS: Neighbor Embedding Stability Scoring
#'
#' Performs dimensionality reduction (t-SNE, UMAP, or PHATE) multiple times to evaluate
#' local neighbor stability across repeated embeddings. This function is useful for
#' evaluating the robustness of low-dimensional embeddings for high-dimensional data.
#'
#' @param GCP Optional numeric vector of neighborhood sizes (e.g., perplexity for t-SNE).
#'            If `NULL`, the function generates a default sequence.
#' @param data A numeric matrix or data frame with rows as observations and columns as features.
#' @param cell_type Optional vector of cell type labels for coloring the embedding plots.
#' @param rareness Logical; if `TRUE`, computes rareness metrics based on neighbor stability.
#' @param data.name Character string used in plot titles to label the dataset.
#' @param method Dimensionality reduction method to use. One of `"tsne"`, `"umap"`, or `"phateR"`.
#' @param initialization Initialization method: `1` for random, `2` for PCA.
#' @param stability_threshold Quantile threshold (default = 0.75) for determining neighbor stability.
#' @param early_stop Logical; if `TRUE`, stops early if global stability saturates.
#'
#' @return A list containing:
#' \describe{
#'   \item{local_stability}{A list of local kNN stability scores across GCP values.}
#'   \item{plot_list_stability}{A list of ggplot2 objects showing embeddings colored by stability.}
#'   \item{global_stability_plot}{A ggplot2 plot showing global stability vs. GCP.}
#'   \item{plot_list_cell_type}{(optional) Embedding plots colored by cell type.}
#'   \item{rareness_mean}{(optional) A plot of the rareness score (mean) if `rareness = TRUE`.}
#' }
#'
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom phateR phate
#' @importFrom RSpectra svds
#' @importFrom BiocNeighbors findKNN
#' @export


NESS <- function(GCP = NULL, data, cell_type = NULL, rareness = FALSE, data.name,
                 method = "tsne", initialization = 1, stability_threshold = 0.75, early_stop = TRUE) {
  set.seed(123)
  eg.v <- svds((data), k = 30)$d
  pc <- max(which(eg.v[1:29] / eg.v[2:30] > 1.1))
  svd.data <- svds(data, k = pc)
  data.denoise <- svd.data$v %*% diag(svd.data$d)

  if (is.null(GCP)) {
    GCP <- round(seq(10, 10 * log(nrow(data)), length.out = 5)) #equally divide into 5 GCPs
  }

  N <- 30
  k <- 50
  global_stability <- c()
  rare.mean <- c() # To store rareness score (mean)
  rare.var <- c() # To store rareness score (variance)

  plot_list_cell_type <- list()  # To store cell type plots
  plot_list_stability <- list()  # To store stability plots
  local_stability <- list()      # To store knn.score for each GCP iteration

  # Determine initialization type for naming
  init_type <- ifelse(initialization == 2, "PCA Initialization", "Random Initialization")

  prev_stability <- 0

  tracking_GCP_length <- 0 #track how many GCP is being processed to avoid

  distance_matrix <- as.matrix(dist(data.denoise, method = "euclidean"))
  distance_matrix_phateR <- distance_matrix
  diag(distance_matrix_phateR) <- 0





  for (j in seq_along(GCP)) {
    tracking_GCP_length <- j #track current GCP index

    perp <- GCP[j]
    set.seed(123 + j)
    knn.mat <- list()
    knn.info <- NULL

    if (method == "tsne") {
      for (i in seq_len(N)) {
        set.seed(123 + i)
        if (i == 1) {
          out <- Rtsne(distance_matrix,
                       is_distance = TRUE,
                       perplexity = perp,
                       check_duplicates = FALSE,
                       pca = (initialization == 2),
                       seed = 123)
          embedding <- out$Y
        } else {
          init_embed <- embedding
          out <- Rtsne(distance_matrix,
                       is_distance = TRUE,
                       perplexity = perp,
                       check_duplicates = FALSE,
                       pca = (initialization == 2),
                       Y_init = init_embed,
                       seed = 123)
          embedding <- out$Y
        }
        knn.mat[[i]] <- findKNN(embedding, k = k)$index
      }
    } else if (method == "umap") {
      for (i in seq_len(N)) {
        if (is.null(knn.info)) {
          set.seed(123 + i)
          out <- uwot::umap(data.denoise,
                            n_neighbors = perp,
                            init = ifelse(initialization == 2, "pca", "lvrandom"),
                            n_threads = 1,
                            seed = 123 + i,
                            ret_nn = TRUE) # Include nearest neighbor information in the output
          # Ensure NN information is present
          knn.info <- list(idx = out$nn$euclidean$idx, dist = out$nn$euclidean$dist)
          embedding <- out$embedding
        } else {
          embedding <- uwot::umap(data.denoise,
                                  n_neighbors = perp,
                                  init = ifelse(initialization == 2, "pca", "lvrandom"),
                                  n_threads = 1,
                                  seed = 123 + i,
                                  nn_method = knn.info,  # Reuse NN info  with proper format
                                  ret_nn = FALSE)
        }
        knn.mat[[i]] <- findKNN(embedding, k = k)$index
      }
    } else if (method == "phateR") {
      phate.init <- NULL
      for (i in seq_len(N)) {
        set.seed(123 + i)
        out <- phate(distance_matrix_phateR,
                     knn.dist.method = "precomputed",
                     ndim = 2,
                     knn = perp,
                     seed = 123 + i,
                     verbose = TRUE,
                     init = phate.init) # reuse PHATE object if available
        if (i == 1) phate.init <- out
        embedding <- out$embedding
        knn.mat[[i]] <- findKNN(embedding, k = 50)$index
      }
    }

    Y <- embedding






    # Accumulate neighbor counts (construct knn graph)
    knn.graph <- matrix(0, nrow = nrow(knn.mat[[1]]), ncol = nrow(knn.mat[[1]]))
    for (i in seq_len(N)) {
      for (l in seq_len(nrow(knn.mat[[1]]))) {
        knn.graph[l, knn.mat[[i]][l, 1:k]] <- knn.graph[l, knn.mat[[i]][l, 1:k]] + 1
      }
      knn.graph <- pmax(knn.graph, t(knn.graph))
    }

    # Compute knn.score using user-defined stability_threshold
    knn.score <- sapply(seq_len(nrow(knn.mat[[1]])), function(i) {
      quantile(knn.graph[i, knn.graph[i, ] != 0] / N, stability_threshold)
    })
    local_stability[[j]] <- knn.score









    #rareness score:
    if (rareness) {
      temp.out <- matrix(ncol = N, nrow = N)

      # Use p and q as inner loop variables
      for (p in 1:(N - 1)) {
        for (q in (p + 1):N) {
          temp <- c()
          for (mn in 1:dim(knn.mat[[1]])[1]) {
            temp[mn] <- length(intersect(knn.mat[[p]][mn, 1:k], knn.mat[[q]][mn, 1:k])) / k
          }
          temp.out[p, q] <- median(temp)
        }
      }

      rareness_mean <- c()
      # Calculate mean rareness for indices 2 to (N-1)
      for (aa in 2:(N - 1)) {
        rareness_mean[aa] <- mean(c(temp.out[1:(aa - 1), aa], temp.out[aa, (aa + 1):N]))
      }
      rareness_mean[1] <- mean(temp.out[1, 2:N])
      rareness_mean[N] <- mean(temp.out[1:(N - 1), N])

      rareness_variance <- c()
      for (aa in 2:(N - 1)) {
        temp_list <- c(temp.out[1:(aa - 1), aa], temp.out[aa, (aa + 1):N])
        rareness_variance[aa] <- sum((temp_list - rareness_mean[aa])^2) / length(temp_list)
      }
      rareness_variance[1] <- sum((temp.out[1, 2:N] - rareness_mean[1])^2) / length(temp.out[1, 2:N])
      rareness_variance[N] <- sum((temp.out[1:(N - 1), N] - rareness_mean[N])^2) / length(temp.out[1:(N - 1), N])

      rare.mean[j] <- rareness_mean[N]
      rare.var[j] <- rareness_variance[N]
    }

    boxplot_rareness_mean <- ggplot(data.frame(GCP = GCP[1:tracking_GCP_length]), aes(y = rare.mean)) +
      geom_boxplot() +
      ggtitle(paste("Rareness Score Boxplot "))








    # Metrics depending on cell_type
    if (!is.null(cell_type)) {
      # Plot by cell type
      plot.data <- data.frame(dim1 = Y[, 1], dim2 = Y[, 2], cell_type = cell_type)
      plot11 <- ggplot(plot.data, aes(dim1, dim2, colour = cell_type)) +
        geom_point(size = 1) +
        ggtitle(paste0(data.name, "_", method, "(", init_type, ")" , "_p", perp, "_pc", pc, "_colored by cell type"))
      plot_list_cell_type[[j]] <- plot11

    }

    # Plot by stability
    plot.data <- data.frame(dim1 = Y[, 1], dim2 = Y[, 2], knn_score = knn.score)
    plot22 <- ggplot(plot.data, aes(dim1, dim2, colour = knn_score)) +
      geom_point(size = 1) +
      ggtitle(paste0(data.name, "_", method, "(", init_type, ")" , "_p", perp, "_pc", pc, "_colored by local stability score"))
    plot_list_stability[[j]] <- plot22

    # Stability metrics
    global_stability[j] <- mean(knn.score, na.rm = TRUE)

    if (early_stop && (global_stability[j] >= 0.9 ||  (j > 1 && (global_stability[j] - prev_stability) / prev_stability <= 0.05))) {
      break  # early stop indication
    }
    prev_stability <- global_stability[j]

  } # End of the for-loop over GCP






  #print(rare.mean)
  #print(rare.var)
  rare.mean.vec <- unlist(rare.mean)
  med_value <- median(rare.mean.vec, na.rm = TRUE)
  med_index <- which(rare.mean.vec == med_value)
  # If nothing matched exactly (e.g., due to float rounding), find closest
  if (length(med_index) == 0) {
    med_index <- which.min(abs(rare.mean.vec - med_value))
  }
  # Print for debugging
  print(paste("Index of median rareness score:", med_index))
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

  message <- paste("The final GCP is", GCP[med_index], "with a global stability of", global_stability[med_index])
  if (early_stop) print(message)

  # Return results
  result_list <- list(
    local_stability = local_stability[[med_index]],
    plot_list_stability = plot_list_stability[[med_index]], #stability score for plot with median rareness
    global_stability_plot = stability_plot
  )

  # Add plots if applicable
  if (!is.null(cell_type)) {
    result_list$plot_list_cell_type <- plot_list_cell_type[[med_index]]
  }
  if (rareness) {
    result_list$rareness_mean <- plot_improved(GCP[1:tracking_GCP_length], rare.mean, "Rareness Score(Mean)", paste0(data.name, method, " - Rareness Score(Mean) (", init_type, ")"))
  }

  # Return the result
  return(result_list)
}

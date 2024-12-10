# Install required packages if not already installed
required_packages <- c("MASS", "mclust", "ggplot2", "tidyverse", "Matrix", 
                       "expm", "dbscan", "cluster", "RSpectra", "Rtsne")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set random seed for reproducibility
set.seed(123)

#' Generate samples from mixture of low-rank Gaussians
#' @param n Number of samples
#' @param ambient_dim Dimension of ambient space
#' @param K Number of components
#' @param ranks Vector of ranks for each component
#' @param mixing_props Vector of mixing proportions
#' @param noise_sd Standard deviation of additive noise (0 for no noise)
#' @return List containing data matrix X and true labels
generate_low_rank_mixture <- function(n, ambient_dim, K, ranks, mixing_props, noise_sd = 0) {
  # Input validation
  stopifnot(
    length(ranks) == K,
    length(mixing_props) == K,
    abs(sum(mixing_props) - 1) < 1e-10,
    all(mixing_props >= 0),
    all(ranks <= ambient_dim)
  )
  
  # Initialize outputs
  X <- matrix(0, n, ambient_dim)
  true_labels <- numeric(n)
  
  # Generate orthonormal basis for each component
  U_stars <- lapply(ranks, function(r) {
    qr.Q(qr(matrix(rnorm(ambient_dim * r), ambient_dim, r)))
  })
  
  # Generate samples according to mixing proportions
  sample_counts <- rmultinom(1, n, mixing_props)[,1]
  current_idx <- 1
  
  for(k in 1:K) {
    if(sample_counts[k] > 0) {
      idx <- current_idx:(current_idx + sample_counts[k] - 1)
      Z <- matrix(rnorm(sample_counts[k] * ranks[k]), sample_counts[k], ranks[k])
      X[idx,] <- Z %*% t(U_stars[[k]])
      
      if(noise_sd > 0) {
        X[idx,] <- X[idx,] + matrix(rnorm(sample_counts[k] * ambient_dim, 0, noise_sd), 
                                    sample_counts[k], ambient_dim)
      }
      
      true_labels[idx] <- k
      current_idx <- current_idx + sample_counts[k]
    }
  }
  
  return(list(
    X = X,                 # Data matrix
    labels = true_labels,  # True component labels
    U = U_stars,          # Generated basis matrices
    counts = sample_counts # Samples per component
  ))
}

#' Perform spectral clustering
#' @param similarity_matrix Similarity matrix
#' @param k Number of clusters
#' @return Vector of cluster assignments
spectral_clustering <- function(similarity_matrix, k) {
  D <- diag(rowSums(similarity_matrix))
  L <- D - similarity_matrix
  eig <- RSpectra::eigs_sym(L, k)
  embedding <- eig$vectors
  kmeans(embedding, centers = k)$cluster
}

#' Evaluate clustering methods and compute ARI
#' @param X Data matrix
#' @param true_labels True cluster labels
#' @param K Number of clusters
#' @return Named vector of ARI values
evaluate_clustering <- function(X, true_labels, K) {
  # K-means clustering
  km <- kmeans(X, centers = K, nstart = 10)
  
  # Hierarchical clustering
  hc <- cutree(hclust(dist(X), method = "ward.D2"), k = K)
  
  # Gaussian Mixture Model
  gmm <- Mclust(X, G = K)
  
  # Spectral clustering
  dist_mat <- as.matrix(dist(X))
  sigma <- mean(dist_mat)
  similarity <- exp(-dist_mat^2 / (2 * sigma^2))
  spec <- spectral_clustering(similarity, K)
  
  # Calculate Adjusted Rand Index for each method
  ari_values <- c(
    kmeans = adjustedRandIndex(true_labels, km$cluster),
    hierarchical = adjustedRandIndex(true_labels, hc),
    gmm = adjustedRandIndex(true_labels, gmm$classification),
    spectral = adjustedRandIndex(true_labels, spec)
  )
  
  return(ari_values)
}

#' Create t-SNE visualization
#' @param X Data matrix
#' @param labels Cluster labels
#' @param title Plot title
#' @return ggplot object
create_tsne_plot <- function(X, labels, title) {
  tsne_result <- Rtsne(X, perplexity = 30, verbose = FALSE)
  
  plot_data <- data.frame(
    x = tsne_result$Y[,1],
    y = tsne_result$Y[,2],
    label = factor(labels)
  )
  
  ggplot(plot_data, aes(x = x, y = y, color = label)) +
    geom_point(alpha = 0.6, size = 1) +
    labs(title = title,
         x = "t-SNE 1",
         y = "t-SNE 2",
         color = "Component") +
    theme_minimal() +
    scale_color_discrete(name = "Component")
}

#' Run experiment with multiple conditions and repetitions
#' @param n_reps Number of repetitions
#' @return Dataframe with results
run_experiment <- function(n_reps = 30) {
  # Fixed parameters
  n_samples <- 1000
  ambient_dim <- 100
  
  # Define experimental conditions
  conditions <- list(
    list(K = 5, ranks = rep(5, 5), noise = 0, condition = "default"),
    list(K = 10, ranks = rep(5, 10), noise = 0, condition = "more_components"),
    list(K = 15, ranks = rep(5, 15), noise = 0, condition = "most_components"),
    list(K = 5, ranks = c(5,5,10,15,20), noise = 0, condition = "mixed_ranks"),
    list(K = 5, ranks = rep(20, 5), noise = 0, condition = "high_ranks"),
    list(K = 5, ranks = rep(5, 5), noise = 0.1, condition = "low_noise"),
    list(K = 5, ranks = rep(5, 5), noise = 0.3, condition = "high_noise")
  )
  
  results <- data.frame()
  
  for(cond in conditions) {
    for(rep in 1:n_reps) {
      cat(sprintf("Running condition: %s, repetition: %d\n", cond$condition, rep))
      
      data <- generate_low_rank_mixture(
        n = n_samples,
        ambient_dim = ambient_dim,
        K = cond$K,
        ranks = cond$ranks,
        mixing_props = rep(1/cond$K, cond$K),
        noise_sd = cond$noise
      )
      
      ari_values <- evaluate_clustering(data$X, data$labels, cond$K)
      
      results <- rbind(results, 
                       data.frame(
                         condition = cond$condition,
                         method = names(ari_values),
                         ari = ari_values,
                         rep = rep
                       ))
    }
  }
  
  return(results)
}

# Run the main experiment
results <- run_experiment(n_reps = 30)

# Calculate summary statistics
summary_stats <- results %>%
  group_by(condition, method) %>%
  summarise(
    mean_ari = mean(ari),
    sd_ari = sd(ari),
    ci_lower = mean_ari - 1.96 * sd_ari/sqrt(n()),
    ci_upper = mean_ari + 1.96 * sd_ari/sqrt(n()),
    .groups = "drop"
  )

# Create visualization plots
p1 <- ggplot(results, aes(x = method, y = ari, fill = method)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  facet_wrap(~condition, scales = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Distribution of ARI Scores by Method and Condition",
    x = "Clustering Method",
    y = "Adjusted Rand Index"
  )

p2 <- ggplot(summary_stats, aes(x = method, y = mean_ari, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~condition, scales = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Mean ARI Scores with 95% Confidence Intervals",
    x = "Clustering Method",
    y = "Mean Adjusted Rand Index"
  )

# Save results and plots
dir.create("results", showWarnings = FALSE)
saveRDS(results, "results/experiment_results.rds")
saveRDS(summary_stats, "results/summary_stats.rds")
ggsave("results/violin_plot.pdf", p1, width = 12, height = 8)
ggsave("results/bar_plot.pdf", p2, width = 12, height = 8)

# Print summary statistics
print(summary_stats)
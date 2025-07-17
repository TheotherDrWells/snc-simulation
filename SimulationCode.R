# ========== USER SETTINGS ==========
out_dir <- "E:/SNC/"  # <--- set this once
n_reps <- 1000        # Set to 1000 for full simulation
set.seed(123)
# ===================================

library(psych)
library(polycor)
library(ggplot2)

# ========== SIMULATION DESIGN ==========
design_grid <- expand.grid(
  factors = c(1, 2, 3),
  items_per_factor = c(6, 12, 18),
  reversals = c(0, 2, 3),
  n = c(250, 500, 1000)
)

design_grid$condition <- 1:nrow(design_grid)

# ========== STORAGE ==========
results <- data.frame()

# ========== SNC FUNCTIONS ==========
snc <- function(R, k = 2, factors = NULL, digits = 3) {
  # Enforce symmetry
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  p <- nrow(R)
  if (p < 3) stop("Correlation matrix must be at least 3×3.")
  if (k < 1 || k >= p) stop(paste0("k must be between 1 and ", p - 1))
  
  snc_values <- numeric(p)
  for (i in 1:p) {
    neighbors <- order(abs(R[i, -i]), decreasing = TRUE)[1:k]
    adjusted <- ifelse(neighbors >= i, neighbors + 1, neighbors)
    snc_values[i] <- mean(abs(R[i, adjusted]))
  }
  
  snc_values <- round(snc_values, digits)
  overall_val <- round(mean(snc_values), digits)
  item_names <- rownames(R)
  if (is.null(item_names)) item_names <- paste0("Item", 1:p)
  
  if (is.null(factors)) {
    item_df <- data.frame(item = item_names, snc = snc_values, stringsAsFactors = FALSE)
    result <- list(overall = overall_val, items = item_df)
  } else {
    if (length(factors) != p) stop("Length of 'factors' must match number of items.")
    item_df <- data.frame(item = item_names, snc = snc_values, factor = as.character(factors), stringsAsFactors = FALSE)
    factor_df <- aggregate(snc ~ factor, data = item_df, FUN = mean)
    factor_df$snc <- round(factor_df$snc, digits)
    result <- list(overall = overall_val, items = item_df, factors = factor_df)
  }
  
  class(result) <- "snc"
  return(result)
}

print.snc <- function(x, ...) {
  cat("Strongest Neighbor Coherence (SNC)\n")
  cat("----------------------------------\n")
  cat("Overall SNC:", x$overall, "\n\n")
  
  if (!is.null(x$factors)) {
    cat("Item-level SNC:\n")
    print(x$items, row.names = FALSE)
    cat("\nFactor-level SNC:\n")
    print(x$factors, row.names = FALSE)
  } else {
    cat("Item-level SNC:\n")
    print(x$items, row.names = FALSE)
  }
}

# ========== SIM FUNCTION ==========
simulate_condition <- function(factors, items_per_factor, reversals, n, rep_id) {
  total_items <- items_per_factor
  items_per_dim <- total_items / factors
  if (items_per_dim %% 1 != 0) return(NULL)  # skip illegal configs
  
  # 1. Create loadings matrix
  Lambda <- matrix(0, nrow = total_items, ncol = factors)
  for (f in 1:factors) {
    start <- (f - 1) * items_per_dim + 1
    end <- f * items_per_dim
    Lambda[start:end, f] <- 0.7
  }
  
  # 2. Simulate data
  true_scores <- MASS::mvrnorm(n, mu = rep(0, factors), Sigma = diag(factors))
  error <- matrix(rnorm(n * total_items, 0, sqrt(1 - 0.49)), nrow = n)
  data <- true_scores %*% t(Lambda) + error
  
  # 3. Convert to ordinal (simulate Likert-style)
  ordinal_data <- apply(data, 2, function(x) cut(x,
                                                 breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                                 labels = FALSE, include.lowest = TRUE
  ))
  
  # 4. Reverse items
  if (reversals > 0) {
    rev_idx <- sample(1:total_items, reversals)
    ordinal_data[, rev_idx] <- 6 - ordinal_data[, rev_idx]  # 1–5 reversed
  }
  
  # 5. Factor assignment vector
  factor_vector <- rep(paste0("F", 1:factors), each = items_per_dim)
  
  # 6. Polychoric correlation matrix
  R <- polycor::hetcor(as.data.frame(ordinal_data))$correlations
  if (is.null(R)) return(NULL)
  
  # 7. Alpha per factor
  alpha_vals <- sapply(unique(factor_vector), function(f) {
    idx <- which(factor_vector == f)
    subR <- R[idx, idx]
    psych::alpha(subR)$total$raw_alpha
  })
  
  # 8. SNC per factor
  snc_obj <- snc(R, factors = factor_vector, digits = 3)
  if (is.null(snc_obj)) return(NULL)
  snc_df <- snc_obj$factors
  
  # 9. Combine
  out <- data.frame(
    rep = rep_id,
    factors = factors,
    items = items_per_factor,
    reversals = reversals,
    n = n,
    factor = snc_df$factor,
    alpha = as.numeric(alpha_vals[snc_df$factor]),
    snc = snc_df$snc
  )
  return(out)
}

# ========== RUN SIMULATION ==========
for (i in 1:nrow(design_grid)) {
  params <- design_grid[i, ]
  for (r in 1:n_reps) {
    sim <- simulate_condition(
      factors = params$factors,
      items_per_factor = params$items_per_factor,
      reversals = params$reversals,
      n = params$n,
      rep_id = r
    )
    if (!is.null(sim)) {
      sim$condition <- params$condition
      results <- rbind(results, sim)
    }
  }
}

# ========== SAVE ==========
write.csv(results, file.path(out_dir, "snc_sim_results1000.csv"), row.names = FALSE)

# ========== AGGREGATE ==========
summary_table <- aggregate(cbind(alpha, snc) ~ factors + items + reversals + n + factor, data = results, FUN = mean, na.rm = TRUE)
write.csv(summary_table, file.path(out_dir, "snc_sim_summary1000.csv"), row.names = FALSE)


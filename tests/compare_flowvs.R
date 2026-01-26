# =============================================================================
# Comparison script: Native flowVS implementation vs original flowVS package
#
# This script compares the cofactor estimation between:
# 1. Native implementation (from main.R)
# 2. Original flowVS package (requires Bioconductor)
#
# Run this script in an R environment with flowVS installed to compare results.
# =============================================================================

# Generate synthetic flow cytometry-like data with multiple populations
generate_synthetic_fcs_data <- function(n_cells = 10000, seed = 42) {
  set.seed(seed)

  # Simulate 5 channels (like CD45, CD3, CD4, CD8, CD19)
  channels <- c("CD45", "CD3", "CD4", "CD8", "CD19")

  # Create data with multiple populations per channel
  # Population structure mimics real flow cytometry data
  data <- list()

  for (ch in channels) {
    # Generate 2-4 populations with different means and variances
    n_pops <- sample(2:4, 1)
    pop_sizes <- rmultinom(1, n_cells, rep(1/n_pops, n_pops))[,1]

    values <- numeric(0)
    for (i in seq_len(n_pops)) {
      # Populations on different scales (mimics negative, dim, bright populations)
      pop_mean <- 10^runif(1, 1, 4)  # means between 10 and 10000
      pop_sd <- pop_mean * runif(1, 0.2, 0.5)  # CV between 20-50%
      pop_values <- rnorm(pop_sizes[i], mean = pop_mean, sd = pop_sd)
      # Add some negative values (background)
      pop_values <- pop_values + rnorm(pop_sizes[i], mean = 0, sd = pop_mean * 0.1)
      values <- c(values, pop_values)
    }

    data[[ch]] <- values
  }

  as.data.frame(data)
}

# =============================================================================
# Native implementation (copied from main.R for standalone testing)
# =============================================================================

find_density_peaks <- function(x, n = 512) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)

  dens <- density(x, n = n, na.rm = TRUE)

  dy <- diff(dens$y)
  peaks_idx <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1
  valleys_idx <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1

  if (length(peaks_idx) == 0) return(NULL)

  peak_heights <- dens$y[peaks_idx]
  min_height <- max(dens$y) * 0.05
  significant_peaks <- peaks_idx[peak_heights >= min_height]

  if (length(significant_peaks) == 0) return(NULL)

  peak_locations <- dens$x[significant_peaks]
  boundaries <- c(min(x), dens$x[valleys_idx], max(x))
  boundaries <- sort(unique(boundaries))

  list(
    peaks = peak_locations,
    boundaries = boundaries,
    density = dens
  )
}

assign_to_peaks <- function(x, peak_info) {
  if (is.null(peak_info)) return(rep(1L, length(x)))

  peaks <- peak_info$peaks
  boundaries <- peak_info$boundaries
  assignments <- rep(NA_integer_, length(x))

  for (i in seq_along(peaks)) {
    lower_bound <- max(boundaries[boundaries <= peaks[i]])
    upper_bound <- min(boundaries[boundaries >= peaks[i]])
    if (i == 1) lower_bound <- min(boundaries)
    if (i == length(peaks)) upper_bound <- max(boundaries)
    in_range <- x >= lower_bound & x < upper_bound
    assignments[in_range] <- i
  }

  remaining <- is.na(assignments)
  if (any(remaining)) {
    for (j in which(remaining)) {
      distances <- abs(x[j] - peaks)
      assignments[j] <- which.min(distances)
    }
  }

  assignments
}

bartlett_statistic <- function(x, groups) {
  valid <- is.finite(x) & !is.na(groups)
  x <- x[valid]
  groups <- groups[valid]

  unique_groups <- unique(groups)
  m <- length(unique_groups)

  if (m < 2) return(Inf)

  n_i <- numeric(m)
  var_i <- numeric(m)

  for (i in seq_along(unique_groups)) {
    group_data <- x[groups == unique_groups[i]]
    n_i[i] <- length(group_data)
    if (n_i[i] < 2) return(Inf)
    var_i[i] <- var(group_data)
  }

  if (any(var_i <= 0) || any(!is.finite(var_i))) return(Inf)

  n <- sum(n_i)
  df_i <- n_i - 1
  var_pooled <- sum(df_i * var_i) / sum(df_i)

  if (var_pooled <= 0 || !is.finite(var_pooled)) return(Inf)

  numerator <- (n - m) * log(var_pooled) - sum(df_i * log(var_i))
  correction <- 1 + (1 / (3 * (m - 1))) * (sum(1 / df_i) - 1 / (n - m))
  B <- numerator / correction

  if (!is.finite(B)) return(Inf)
  B
}

flowvs_objective <- function(cofactor, x) {
  if (cofactor <= 0) return(Inf)

  transformed <- asinh(x / cofactor)
  peak_info <- find_density_peaks(transformed)

  if (is.null(peak_info) || length(peak_info$peaks) < 2) {
    return(1e10)
  }

  groups <- assign_to_peaks(transformed, peak_info)
  B <- bartlett_statistic(transformed, groups)
  B
}

estimate_cofactor_native <- function(values, cf_low = -1, cf_high = 10) {
  values <- values[is.finite(values)]

  if (length(values) < 100) {
    return(5)
  }

  intervals <- seq(cf_low, cf_high, by = 1)
  best_cofactor <- 5
  best_stat <- Inf

  for (i in seq_len(length(intervals) - 1)) {
    low <- 10^intervals[i]
    high <- 10^intervals[i + 1]

    result <- tryCatch({
      optimize(
        f = flowvs_objective,
        interval = c(low, high),
        x = values,
        tol = 0.01 * (high - low)
      )
    }, error = function(e) NULL)

    if (!is.null(result) && is.finite(result$objective)) {
      if (result$objective < best_stat) {
        best_stat <- result$objective
        best_cofactor <- result$minimum
      }
    }
  }

  if (!is.finite(best_stat) || best_stat > 1e9) {
    positive_values <- values[values > 0]
    if (length(positive_values) >= 10) {
      best_cofactor <- as.numeric(quantile(positive_values, 0.05, na.rm = TRUE))
      best_cofactor <- max(0.1, min(1e10, best_cofactor))
    }
  }

  best_cofactor
}

# =============================================================================
# Run comparison
# =============================================================================

cat("=== flowVS Algorithm Comparison ===\n\n")

# Generate test data
cat("Generating synthetic flow cytometry data...\n")
test_data <- generate_synthetic_fcs_data(n_cells = 10000, seed = 42)
cat(sprintf("Generated %d cells with %d channels\n\n", nrow(test_data), ncol(test_data)))

# Run native implementation
cat("Running native implementation...\n")
native_results <- sapply(test_data, estimate_cofactor_native)
cat("Native cofactors:\n")
print(round(native_results, 2))
cat("\n")

# Try to run original flowVS if available
cat("Checking for flowVS package...\n")
flowvs_available <- requireNamespace("flowVS", quietly = TRUE) &&
                    requireNamespace("flowCore", quietly = TRUE)

if (flowvs_available) {
  cat("flowVS package found! Running original implementation...\n")

  library(flowCore)
  library(flowVS)

  # Convert to flowFrame
  expr_matrix <- as.matrix(test_data)
  ff <- flowFrame(exprs = expr_matrix)
  fs <- flowSet(list(sample1 = ff))

  # Run flowVS
  channels <- colnames(test_data)
  flowvs_results <- estParamFlowVS(fs, channels)

  cat("flowVS cofactors:\n")
  print(round(flowvs_results, 2))
  cat("\n")

  # Compare results
  cat("=== Comparison ===\n")
  comparison <- data.frame(
    Channel = channels,
    Native = round(native_results, 2),
    flowVS = round(flowvs_results, 2),
    Ratio = round(native_results / flowvs_results, 3),
    LogDiff = round(log10(native_results) - log10(flowvs_results), 3)
  )
  print(comparison)

  cat("\nSummary statistics:\n")
  cat(sprintf("Mean absolute log10 difference: %.3f\n", mean(abs(comparison$LogDiff))))
  cat(sprintf("Max absolute log10 difference: %.3f\n", max(abs(comparison$LogDiff))))
  cat(sprintf("Correlation (log scale): %.3f\n", cor(log10(native_results), log10(flowvs_results))))

} else {
  cat("flowVS package not available.\n")
  cat("To install flowVS, run:\n")
  cat("  BiocManager::install('flowVS')\n\n")

  cat("Native implementation results only:\n")
  print(data.frame(
    Channel = names(native_results),
    Cofactor = round(native_results, 2),
    Log10_Cofactor = round(log10(native_results), 2)
  ))
}

# =============================================================================
# Test with HD dataset if flowVS is available
# =============================================================================

if (flowvs_available) {
  cat("\n\n=== Testing with flowVS HD Dataset ===\n")

  tryCatch({
    data(HD, package = "flowVS")

    # Get first sample for testing
    ff_hd <- HD[[1]]
    channels_hd <- c("CD45", "CD3", "CD4", "CD8", "CD19")

    # Check which channels are available
    available_channels <- intersect(channels_hd, colnames(exprs(ff_hd)))
    cat(sprintf("Available channels: %s\n", paste(available_channels, collapse = ", ")))

    if (length(available_channels) > 0) {
      # Extract data
      hd_data <- as.data.frame(exprs(ff_hd)[, available_channels])

      # Run native
      cat("\nRunning native implementation on HD data...\n")
      native_hd <- sapply(hd_data, estimate_cofactor_native)

      # Run flowVS
      cat("Running flowVS on HD data...\n")
      flowvs_hd <- estParamFlowVS(HD, available_channels)

      # Compare
      cat("\n=== HD Dataset Comparison ===\n")
      hd_comparison <- data.frame(
        Channel = available_channels,
        Native = round(native_hd, 2),
        flowVS = round(flowvs_hd, 2),
        Ratio = round(native_hd / flowvs_hd, 3),
        LogDiff = round(log10(native_hd) - log10(flowvs_hd), 3)
      )
      print(hd_comparison)

      cat("\nHD Summary statistics:\n")
      cat(sprintf("Mean absolute log10 difference: %.3f\n", mean(abs(hd_comparison$LogDiff))))
      cat(sprintf("Correlation (log scale): %.3f\n", cor(log10(native_hd), log10(flowvs_hd))))
    }

  }, error = function(e) {
    cat(sprintf("Error loading HD dataset: %s\n", e$message))
  })
}

cat("\n=== Comparison Complete ===\n")

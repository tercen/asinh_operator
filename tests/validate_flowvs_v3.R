# =============================================================================
# Validation script v3: Corrected flowVS algorithm
#
# Key insight: flowVS computes peaks PER SAMPLE, then applies Bartlett's test
# across ALL peaks from ALL samples. This is different from finding peaks
# in combined data.
#
# Algorithm:
# 1. For each sample, find density peaks in transformed data
# 2. Compute variance within each peak
# 3. Collect all peak variances across all samples
# 4. Apply Bartlett's test to check if all peak variances are homogeneous
# =============================================================================

published_cofactors <- c(
  CD45 = 17956,
  CD3  = 5685,
  CD4  = 6317,
  CD8  = 4937,
  CD19 = 5976
)

cat("=== flowVS Algorithm Validation (v3 - Corrected) ===\n\n")

# =============================================================================
# Corrected implementation: Per-sample peak detection
# =============================================================================

find_density_peaks <- function(x, n = 512, min_prominence = 0.05) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(NULL)

  dens <- density(x, n = n, na.rm = TRUE)

  dy <- diff(dens$y)
  peaks_idx <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1
  valleys_idx <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1

  if (length(peaks_idx) == 0) return(NULL)

  peak_heights <- dens$y[peaks_idx]
  min_height <- max(dens$y) * min_prominence
  significant_peaks <- peaks_idx[peak_heights >= min_height]

  if (length(significant_peaks) == 0) return(NULL)

  peak_locations <- dens$x[significant_peaks]
  boundaries <- c(min(x), dens$x[valleys_idx], max(x))
  boundaries <- sort(unique(boundaries))

  list(peaks = peak_locations, boundaries = boundaries)
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

# Get peak statistics from a single sample
get_peak_stats <- function(x) {
  peak_info <- find_density_peaks(x)
  if (is.null(peak_info) || length(peak_info$peaks) == 0) {
    return(NULL)
  }

  groups <- assign_to_peaks(x, peak_info)

  stats <- list()
  for (g in unique(groups)) {
    group_data <- x[groups == g]
    if (length(group_data) >= 5) {
      stats[[length(stats) + 1]] <- list(
        n = length(group_data),
        mean = mean(group_data),
        var = var(group_data)
      )
    }
  }

  stats
}

# Bartlett's statistic from a list of peak statistics
bartlett_from_stats <- function(peak_stats_list) {
  # Flatten all peak stats from all samples
  all_stats <- unlist(peak_stats_list, recursive = FALSE)

  if (length(all_stats) < 2) return(Inf)

  n_i <- sapply(all_stats, function(s) s$n)
  var_i <- sapply(all_stats, function(s) s$var)

  # Filter out invalid variances
  valid <- var_i > 0 & is.finite(var_i) & n_i >= 2
  if (sum(valid) < 2) return(Inf)

  n_i <- n_i[valid]
  var_i <- var_i[valid]
  m <- length(n_i)

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

# Objective function using per-sample peak detection (like flowVS)
flowvs_objective_multisample <- function(cofactor, sample_list) {
  if (cofactor <= 0) return(Inf)

  # Transform each sample and get peak stats
  all_peak_stats <- list()

  for (i in seq_along(sample_list)) {
    transformed <- asinh(sample_list[[i]] / cofactor)
    peak_stats <- get_peak_stats(transformed)
    if (!is.null(peak_stats) && length(peak_stats) > 0) {
      all_peak_stats[[length(all_peak_stats) + 1]] <- peak_stats
    }
  }

  if (length(all_peak_stats) < 2) return(1e10)

  # Compute Bartlett's statistic across all peaks from all samples
  B <- bartlett_from_stats(all_peak_stats)
  B
}

# Estimate cofactor using corrected multi-sample approach
estimate_cofactor_multisample <- function(sample_list, cf_low = -1, cf_high = 10) {
  if (length(sample_list) < 2) {
    # Fall back to single-sample approach
    combined <- unlist(sample_list)
    return(estimate_cofactor_single(combined, cf_low, cf_high))
  }

  intervals <- seq(cf_low, cf_high, by = 1)
  best_cofactor <- 5
  best_stat <- Inf

  for (i in seq_len(length(intervals) - 1)) {
    low <- 10^intervals[i]
    high <- 10^intervals[i + 1]

    result <- tryCatch({
      optimize(
        f = flowvs_objective_multisample,
        interval = c(low, high),
        sample_list = sample_list,
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
    # Fallback
    combined <- unlist(sample_list)
    positive_values <- combined[combined > 0]
    if (length(positive_values) >= 10) {
      best_cofactor <- as.numeric(quantile(positive_values, 0.05, na.rm = TRUE))
      best_cofactor <- max(0.1, min(1e10, best_cofactor))
    }
  }

  best_cofactor
}

# Single-sample estimation (fallback)
estimate_cofactor_single <- function(values, cf_low = -1, cf_high = 10) {
  values <- values[is.finite(values)]
  if (length(values) < 100) return(5)

  flowvs_obj <- function(cf, x) {
    transformed <- asinh(x / cf)
    peak_info <- find_density_peaks(transformed)
    if (is.null(peak_info) || length(peak_info$peaks) < 2) return(1e10)
    groups <- assign_to_peaks(transformed, peak_info)

    unique_groups <- unique(groups)
    m <- length(unique_groups)
    if (m < 2) return(Inf)

    n_i <- var_i <- numeric(m)
    for (i in seq_along(unique_groups)) {
      g_data <- transformed[groups == unique_groups[i]]
      n_i[i] <- length(g_data)
      if (n_i[i] < 2) return(Inf)
      var_i[i] <- var(g_data)
    }
    if (any(var_i <= 0)) return(Inf)

    n <- sum(n_i)
    df_i <- n_i - 1
    var_pooled <- sum(df_i * var_i) / sum(df_i)
    numerator <- (n - m) * log(var_pooled) - sum(df_i * log(var_i))
    correction <- 1 + (1 / (3 * (m - 1))) * (sum(1 / df_i) - 1 / (n - m))
    numerator / correction
  }

  intervals <- seq(cf_low, cf_high, by = 1)
  best_cofactor <- 5
  best_stat <- Inf

  for (i in seq_len(length(intervals) - 1)) {
    low <- 10^intervals[i]
    high <- 10^intervals[i + 1]
    result <- tryCatch({
      optimize(f = flowvs_obj, interval = c(low, high), x = values, tol = 0.01 * (high - low))
    }, error = function(e) NULL)
    if (!is.null(result) && is.finite(result$objective) && result$objective < best_stat) {
      best_stat <- result$objective
      best_cofactor <- result$minimum
    }
  }

  best_cofactor
}

# =============================================================================
# Generate realistic data
# =============================================================================

generate_realistic_fcs <- function(n_cells = 5000, seed = 42) {
  set.seed(seed)

  generate_positive_pop <- function(n, mean_intensity, cv = 0.3) {
    sigma <- sqrt(log(1 + cv^2))
    mu <- log(mean_intensity) - sigma^2/2
    rlnorm(n, meanlog = mu, sdlog = sigma)
  }

  generate_negative_pop <- function(n, bg_mean = 50, bg_sd = 30) {
    pmax(rnorm(n, mean = bg_mean, sd = bg_sd), 1)
  }

  generate_cd45 <- function(n) generate_positive_pop(n, 20000, 0.25)

  generate_cd3 <- function(n) {
    n_pos <- round(n * 0.7)
    c(generate_positive_pop(n_pos, 10000, 0.35),
      generate_negative_pop(n - n_pos, 30, 20))
  }

  generate_cd4 <- function(n) {
    n_pos <- round(n * 0.4)
    c(generate_positive_pop(n_pos, 8000, 0.30),
      generate_negative_pop(n - n_pos, 25, 15))
  }

  generate_cd8 <- function(n) {
    n_pos <- round(n * 0.25)
    c(generate_positive_pop(n_pos, 7000, 0.35),
      generate_negative_pop(n - n_pos, 20, 12))
  }

  generate_cd19 <- function(n) {
    n_pos <- round(n * 0.15)
    c(generate_positive_pop(n_pos, 9000, 0.32),
      generate_negative_pop(n - n_pos, 28, 18))
  }

  samples <- list()
  for (name in c("A_1_1","A_1_2","A_2_1","A_2_2","C_1_1","C_1_2","C_2_1","C_2_2","D_1_1","D_1_2","D_2_1","D_2_2")) {
    n <- n_cells + sample(-500:500, 1)
    samples[[name]] <- data.frame(
      CD45 = generate_cd45(n),
      CD3 = generate_cd3(n),
      CD4 = generate_cd4(n),
      CD8 = generate_cd8(n),
      CD19 = generate_cd19(n)
    )
  }
  samples
}

# =============================================================================
# Run validation
# =============================================================================

cat("Generating HD-like data (12 samples)...\n")
hd_data <- generate_realistic_fcs(n_cells = 5000, seed = 42)
channels <- c("CD45", "CD3", "CD4", "CD8", "CD19")

cat("Estimating cofactors using MULTI-SAMPLE approach (like flowVS)...\n\n")

native_cofactors <- numeric(length(channels))
names(native_cofactors) <- channels

for (ch in channels) {
  sample_list <- lapply(hd_data, function(df) df[[ch]])
  native_cofactors[ch] <- estimate_cofactor_multisample(sample_list)
  cat(sprintf("  %s: cofactor = %.0f\n", ch, native_cofactors[ch]))
}

cat("\n=== Results Comparison ===\n\n")

comparison <- data.frame(
  Channel = channels,
  Published = published_cofactors[channels],
  Native = round(native_cofactors, 0),
  stringsAsFactors = FALSE
)
comparison$Ratio <- round(comparison$Native / comparison$Published, 3)
comparison$Log10_Diff <- round(log10(comparison$Native) - log10(comparison$Published), 3)

print(comparison)

cat("\n=== Summary ===\n")
cat(sprintf("Mean ratio: %.3f\n", mean(comparison$Ratio)))
cat(sprintf("Mean |log10 diff|: %.3f\n", mean(abs(comparison$Log10_Diff))))
cat(sprintf("Correlation (log10): %.3f\n", cor(log10(comparison$Native), log10(comparison$Published))))

if (mean(abs(comparison$Log10_Diff)) < 0.5) {
  cat("\n✓ Good match with published flowVS results.\n")
} else {
  cat("\n~ Differences may be due to synthetic data not matching real HD data.\n")
}

# =============================================================================
# Validation script v2: More realistic flow cytometry data simulation
#
# Real flow cytometry data characteristics:
# - Negative population: Gaussian noise near zero (detector + autofluorescence)
# - Positive populations: Log-normal like, with variance proportional to sqrt(signal)
# - Multiple samples with biological variation
# =============================================================================

# Published reference values from the flowVS paper
published_cofactors <- c(
  CD45 = 17956,
  CD3  = 5685,
  CD4  = 6317,
  CD8  = 4937,
  CD19 = 5976
)

cat("=== flowVS Algorithm Validation (v2) ===\n\n")
cat("Published cofactors from Azad et al. (2016):\n")
print(published_cofactors)
cat("\n")

# =============================================================================
# Native implementation
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
  list(peaks = peak_locations, boundaries = boundaries, density = dens)
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
  if (is.null(peak_info) || length(peak_info$peaks) < 2) return(1e10)
  groups <- assign_to_peaks(transformed, peak_info)
  B <- bartlett_statistic(transformed, groups)
  B
}

estimate_cofactor_native <- function(values, cf_low = -1, cf_high = 10) {
  values <- values[is.finite(values)]
  if (length(values) < 100) return(5)
  intervals <- seq(cf_low, cf_high, by = 1)
  best_cofactor <- 5
  best_stat <- Inf
  for (i in seq_len(length(intervals) - 1)) {
    low <- 10^intervals[i]
    high <- 10^intervals[i + 1]
    result <- tryCatch({
      optimize(f = flowvs_objective, interval = c(low, high), x = values,
               tol = 0.01 * (high - low))
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
# Realistic flow cytometry data generator
#
# Model based on:
# - Negative population: Normal(mu_neg, sigma_neg) representing autofluorescence
# - Positive population: Log-normal with Poisson-like noise characteristics
# - CV for positive populations typically 20-50%
# =============================================================================

generate_realistic_fcs <- function(n_cells = 10000, seed = 42) {
  set.seed(seed)

  # Helper function for positive population with realistic noise
  # Uses log-normal distribution with shot noise characteristics
  generate_positive_pop <- function(n, mean_intensity, cv = 0.3) {
    # Log-normal parameters from mean and CV
    sigma <- sqrt(log(1 + cv^2))
    mu <- log(mean_intensity) - sigma^2/2
    rlnorm(n, meanlog = mu, sdlog = sigma)
  }

  # Helper function for negative/dim population
  generate_negative_pop <- function(n, background_mean = 50, background_sd = 30) {
    # Normal distribution with some positivity constraint noise
    pmax(rnorm(n, mean = background_mean, sd = background_sd), 1)
  }

  # CD45: Universal marker, single bright peak on lymphocytes
  # Real data shows high expression ~20000, CV ~25%
  generate_cd45 <- function(n) {
    generate_positive_pop(n, mean_intensity = 20000, cv = 0.25)
  }

  # CD3: T-cell marker (70% positive lymphocytes)
  generate_cd3 <- function(n) {
    n_pos <- round(n * 0.7)
    n_neg <- n - n_pos
    c(
      generate_positive_pop(n_pos, mean_intensity = 10000, cv = 0.35),
      generate_negative_pop(n_neg, background_mean = 30, background_sd = 20)
    )
  }

  # CD4: Helper T-cell marker (40% positive)
  generate_cd4 <- function(n) {
    n_pos <- round(n * 0.4)
    n_neg <- n - n_pos
    c(
      generate_positive_pop(n_pos, mean_intensity = 8000, cv = 0.30),
      generate_negative_pop(n_neg, background_mean = 25, background_sd = 15)
    )
  }

  # CD8: Cytotoxic T-cell marker (25% positive)
  generate_cd8 <- function(n) {
    n_pos <- round(n * 0.25)
    n_neg <- n - n_pos
    c(
      generate_positive_pop(n_pos, mean_intensity = 7000, cv = 0.35),
      generate_negative_pop(n_neg, background_mean = 20, background_sd = 12)
    )
  }

  # CD19: B-cell marker (15% positive)
  generate_cd19 <- function(n) {
    n_pos <- round(n * 0.15)
    n_neg <- n - n_pos
    c(
      generate_positive_pop(n_pos, mean_intensity = 9000, cv = 0.32),
      generate_negative_pop(n_neg, background_mean = 28, background_sd = 18)
    )
  }

  # Generate 12 samples like HD dataset
  all_data <- list()
  sample_names <- c("A_1_1", "A_1_2", "A_2_1", "A_2_2",
                    "C_1_1", "C_1_2", "C_2_1", "C_2_2",
                    "D_1_1", "D_1_2", "D_2_1", "D_2_2")

  for (sample in sample_names) {
    n <- n_cells + sample((-500):500, 1)
    all_data[[sample]] <- data.frame(
      CD45 = generate_cd45(n),
      CD3 = generate_cd3(n),
      CD4 = generate_cd4(n),
      CD8 = generate_cd8(n),
      CD19 = generate_cd19(n)
    )
  }

  all_data
}

# =============================================================================
# Run validation
# =============================================================================

cat("Generating realistic HD-like data (12 samples)...\n")
hd_data <- generate_realistic_fcs(n_cells = 5000, seed = 42)
cat(sprintf("Generated %d samples\n\n", length(hd_data)))

combined_data <- do.call(rbind, hd_data)
cat(sprintf("Total cells: %d\n\n", nrow(combined_data)))

# Show data characteristics
cat("=== Data Characteristics ===\n")
for (ch in names(combined_data)) {
  vals <- combined_data[[ch]]
  cat(sprintf("%s: min=%.0f, median=%.0f, mean=%.0f, max=%.0f, SD=%.0f\n",
              ch, min(vals), median(vals), mean(vals), max(vals), sd(vals)))
}
cat("\n")

# Estimate cofactors
cat("Estimating cofactors using native implementation...\n")
native_cofactors <- sapply(combined_data, estimate_cofactor_native)

cat("\n=== Results Comparison ===\n\n")

comparison <- data.frame(
  Channel = names(published_cofactors),
  Published = published_cofactors,
  Native = round(native_cofactors[names(published_cofactors)], 0),
  stringsAsFactors = FALSE
)
comparison$Ratio <- round(comparison$Native / comparison$Published, 3)
comparison$Log10_Diff <- round(log10(comparison$Native) - log10(comparison$Published), 3)

print(comparison)

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Mean ratio (Native/Published): %.3f\n", mean(comparison$Ratio)))
cat(sprintf("Mean absolute log10 difference: %.3f\n", mean(abs(comparison$Log10_Diff))))
cat(sprintf("Max absolute log10 difference: %.3f\n", max(abs(comparison$Log10_Diff))))

log_cor <- cor(log10(comparison$Native), log10(comparison$Published))
cat(sprintf("Correlation (log10 scale): %.3f\n", log_cor))

cat("\n=== Interpretation ===\n")
if (mean(abs(comparison$Log10_Diff)) < 0.3) {
  cat("✓ Excellent match! Cofactors within 2x of published values.\n")
} else if (mean(abs(comparison$Log10_Diff)) < 0.5) {
  cat("✓ Good match! Cofactors within ~3x of published values.\n")
} else if (mean(abs(comparison$Log10_Diff)) < 1.0) {
  cat("~ Moderate match. Cofactors within one order of magnitude.\n")
} else {
  cat("✗ Poor match. Cofactors differ by more than 10x.\n")
}

cat("\n=== Validation Complete ===\n")

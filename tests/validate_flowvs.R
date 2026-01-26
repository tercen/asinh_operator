# =============================================================================
# Validation script: Compare native implementation against published flowVS results
#
# Published cofactors from Azad et al. (2016) BMC Bioinformatics:
# "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry"
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4964071/
#
# HD dataset cofactors from paper (Figure 4):
# - CD45: 17,956
# - CD3:  5,685
# - CD4:  6,317
# - CD8:  4,937
# - CD19: 5,976
# =============================================================================

# Published reference values from the flowVS paper
published_cofactors <- c(
  CD45 = 17956,
  CD3  = 5685,
  CD4  = 6317,
  CD8  = 4937,
  CD19 = 5976
)

cat("=== flowVS Algorithm Validation ===\n\n")
cat("Published cofactors from Azad et al. (2016):\n")
print(published_cofactors)
cat("\n")

# =============================================================================
# Native implementation (copied from main.R)
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
# Generate realistic HD-like data based on paper description
#
# The HD data consists of:
# - 12 samples from 3 healthy individuals (A, C, D)
# - 2 days per individual, 2 replicates per day
# - Lymphocytes pre-gated
# - Channels: CD45, CD3, CD4, CD8, CD19
# =============================================================================

generate_hd_like_data <- function(n_cells_per_sample = 5000, seed = 42) {
  set.seed(seed)

  # CD45: Universal leukocyte marker - single bright population
  # Very high expression on all lymphocytes
  generate_cd45 <- function(n) {
    # Single bright population with high expression
    rnorm(n, mean = 18000, sd = 3000)
  }

  # CD3: T-cell marker - bimodal (T-cells vs non-T-cells)
  generate_cd3 <- function(n) {
    # ~70% T-cells (bright), ~30% non-T-cells (dim/negative)
    n_tcells <- round(n * 0.7)
    n_other <- n - n_tcells
    c(
      rnorm(n_tcells, mean = 8000, sd = 2000),  # T-cells
      rnorm(n_other, mean = 200, sd = 100)       # Non-T-cells
    )
  }

  # CD4: Helper T-cell marker - bimodal
  generate_cd4 <- function(n) {
    # ~40% CD4+ T-cells, ~60% CD4-
    n_cd4pos <- round(n * 0.4)
    n_cd4neg <- n - n_cd4pos
    c(
      rnorm(n_cd4pos, mean = 7000, sd = 1500),
      rnorm(n_cd4neg, mean = 150, sd = 80)
    )
  }

  # CD8: Cytotoxic T-cell marker - bimodal
  generate_cd8 <- function(n) {
    # ~25% CD8+ T-cells, ~75% CD8-
    n_cd8pos <- round(n * 0.25)
    n_cd8neg <- n - n_cd8pos
    c(
      rnorm(n_cd8pos, mean = 6000, sd = 1200),
      rnorm(n_cd8neg, mean = 100, sd = 60)
    )
  }

  # CD19: B-cell marker - bimodal
  generate_cd19 <- function(n) {
    # ~15% B-cells, ~85% non-B-cells
    n_bcells <- round(n * 0.15)
    n_other <- n - n_bcells
    c(
      rnorm(n_bcells, mean = 7500, sd = 1800),
      rnorm(n_other, mean = 120, sd = 70)
    )
  }

  # Generate 12 samples
  all_data <- list()
  sample_names <- c("A_1_1", "A_1_2", "A_2_1", "A_2_2",
                    "C_1_1", "C_1_2", "C_2_1", "C_2_2",
                    "D_1_1", "D_1_2", "D_2_1", "D_2_2")

  for (sample in sample_names) {
    # Add some variation between samples
    n <- n_cells_per_sample + sample((-500):500, 1)

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

cat("Generating HD-like synthetic data (12 samples)...\n")
hd_data <- generate_hd_like_data(n_cells_per_sample = 5000, seed = 42)
cat(sprintf("Generated %d samples\n\n", length(hd_data)))

# Combine all samples for cofactor estimation (like flowVS does)
combined_data <- do.call(rbind, hd_data)
cat(sprintf("Total cells: %d\n\n", nrow(combined_data)))

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

# Correlation on log scale
log_cor <- cor(log10(comparison$Native), log10(comparison$Published))
cat(sprintf("Correlation (log10 scale): %.3f\n", log_cor))

cat("\n=== Interpretation ===\n")
cat("A log10 difference of 0.3 means ~2x difference in cofactor.\n")
cat("A log10 difference of 1.0 means 10x difference in cofactor.\n")
cat("Cofactors within 1 order of magnitude are generally acceptable\n")
cat("since the exact value depends on data characteristics.\n")

if (mean(abs(comparison$Log10_Diff)) < 0.5) {
  cat("\n✓ Native implementation produces similar cofactors to published flowVS.\n")
} else if (mean(abs(comparison$Log10_Diff)) < 1.0) {
  cat("\n~ Native implementation produces moderately different cofactors.\n")
  cat("  This may be due to differences in peak detection or data simulation.\n")
} else {
  cat("\n✗ Native implementation produces significantly different cofactors.\n")
  cat("  Review the algorithm implementation.\n")
}

cat("\n=== Validation Complete ===\n")

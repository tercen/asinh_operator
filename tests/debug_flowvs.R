# =============================================================================
# Debug script: Investigate differences in cofactor estimation
# =============================================================================

# Native implementation functions
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
    density = dens,
    all_peaks_idx = peaks_idx,
    significant_peaks_idx = significant_peaks
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

  if (m < 2) return(list(stat = Inf, n_groups = m, group_vars = NULL))

  n_i <- numeric(m)
  var_i <- numeric(m)

  for (i in seq_along(unique_groups)) {
    group_data <- x[groups == unique_groups[i]]
    n_i[i] <- length(group_data)
    if (n_i[i] < 2) return(list(stat = Inf, n_groups = m, group_vars = NULL))
    var_i[i] <- var(group_data)
  }

  if (any(var_i <= 0) || any(!is.finite(var_i))) return(list(stat = Inf, n_groups = m, group_vars = var_i))

  n <- sum(n_i)
  df_i <- n_i - 1
  var_pooled <- sum(df_i * var_i) / sum(df_i)

  if (var_pooled <= 0 || !is.finite(var_pooled)) return(list(stat = Inf, n_groups = m, group_vars = var_i))

  numerator <- (n - m) * log(var_pooled) - sum(df_i * log(var_i))
  correction <- 1 + (1 / (3 * (m - 1))) * (sum(1 / df_i) - 1 / (n - m))
  B <- numerator / correction

  list(
    stat = if (is.finite(B)) B else Inf,
    n_groups = m,
    group_sizes = n_i,
    group_vars = var_i,
    var_pooled = var_pooled
  )
}

flowvs_objective_debug <- function(cofactor, x) {
  if (cofactor <= 0) return(list(stat = Inf, n_peaks = 0))

  transformed <- asinh(x / cofactor)
  peak_info <- find_density_peaks(transformed)

  if (is.null(peak_info) || length(peak_info$peaks) < 2) {
    return(list(stat = 1e10, n_peaks = if(is.null(peak_info)) 0 else length(peak_info$peaks)))
  }

  groups <- assign_to_peaks(transformed, peak_info)
  B_info <- bartlett_statistic(transformed, groups)

  list(
    stat = B_info$stat,
    n_peaks = length(peak_info$peaks),
    peak_locations = peak_info$peaks,
    group_vars = B_info$group_vars
  )
}

# Generate realistic data
set.seed(42)

# Generate CD3 data (bimodal - T-cells vs non-T-cells)
n <- 10000
n_tcells <- round(n * 0.7)
n_other <- n - n_tcells
cd3_data <- c(
  rnorm(n_tcells, mean = 8000, sd = 2000),  # T-cells (bright)
  rnorm(n_other, mean = 200, sd = 100)       # Non-T-cells (dim)
)

cat("=== Debugging CD3 Channel ===\n\n")
cat(sprintf("Data range: %.1f to %.1f\n", min(cd3_data), max(cd3_data)))
cat(sprintf("Data mean: %.1f, SD: %.1f\n\n", mean(cd3_data), sd(cd3_data)))

# Test different cofactors
test_cofactors <- c(100, 200, 500, 1000, 2000, 5000, 5685, 10000, 20000)

cat("Cofactor | Bartlett Stat | N Peaks | Peak Locations\n")
cat("---------|---------------|---------|---------------\n")

for (cf in test_cofactors) {
  result <- flowvs_objective_debug(cf, cd3_data)
  peaks_str <- if (result$n_peaks > 0) {
    paste(round(result$peak_locations, 2), collapse = ", ")
  } else {
    "none"
  }
  cat(sprintf("%8d | %13.1f | %7d | %s\n",
              cf, result$stat, result$n_peaks, peaks_str))
}

# Find optimal cofactor
cat("\n=== Searching for optimal cofactor ===\n")

intervals <- seq(-1, 5, by = 0.5)
best_cofactor <- 5
best_stat <- Inf

for (i in seq_len(length(intervals) - 1)) {
  low <- 10^intervals[i]
  high <- 10^intervals[i + 1]

  result <- tryCatch({
    optimize(
      f = function(cf) flowvs_objective_debug(cf, cd3_data)$stat,
      interval = c(low, high),
      tol = 0.01 * (high - low)
    )
  }, error = function(e) NULL)

  if (!is.null(result) && is.finite(result$objective)) {
    cat(sprintf("Interval [%.1f, %.1f]: optimal = %.1f, stat = %.1f\n",
                low, high, result$minimum, result$objective))
    if (result$objective < best_stat) {
      best_stat <- result$objective
      best_cofactor <- result$minimum
    }
  }
}

cat(sprintf("\nBest cofactor found: %.1f (Bartlett stat: %.1f)\n", best_cofactor, best_stat))
cat(sprintf("Published flowVS cofactor for CD3: 5685\n"))

# Detailed analysis at published cofactor
cat("\n=== Analysis at published cofactor (5685) ===\n")
result_published <- flowvs_objective_debug(5685, cd3_data)
cat(sprintf("Bartlett statistic: %.1f\n", result_published$stat))
cat(sprintf("Number of peaks: %d\n", result_published$n_peaks))
if (result_published$n_peaks > 0) {
  cat(sprintf("Peak locations (asinh scale): %s\n",
              paste(round(result_published$peak_locations, 3), collapse = ", ")))
}

# Detailed analysis at our optimal cofactor
cat(sprintf("\n=== Analysis at native optimal cofactor (%.1f) ===\n", best_cofactor))
result_native <- flowvs_objective_debug(best_cofactor, cd3_data)
cat(sprintf("Bartlett statistic: %.1f\n", result_native$stat))
cat(sprintf("Number of peaks: %d\n", result_native$n_peaks))
if (result_native$n_peaks > 0) {
  cat(sprintf("Peak locations (asinh scale): %s\n",
              paste(round(result_native$peak_locations, 3), collapse = ", ")))
}

# The key insight: Show variance at different cofactors
cat("\n=== Variance analysis ===\n")
cat("At published cofactor 5685, checking if variances are actually more homogeneous:\n")

for (cf in c(best_cofactor, 5685)) {
  transformed <- asinh(cd3_data / cf)
  peak_info <- find_density_peaks(transformed)

  if (!is.null(peak_info) && length(peak_info$peaks) >= 2) {
    groups <- assign_to_peaks(transformed, peak_info)
    B_info <- bartlett_statistic(transformed, groups)

    cat(sprintf("\nCofactor %.1f:\n", cf))
    cat(sprintf("  Group variances: %s\n", paste(round(B_info$group_vars, 4), collapse = ", ")))
    cat(sprintf("  Variance ratio (max/min): %.2f\n", max(B_info$group_vars) / min(B_info$group_vars)))
    cat(sprintf("  Bartlett stat: %.1f\n", B_info$stat))
  }
}

library(tercenApi)
library(tercen)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)

ctx <- tercenCtx()
method <- ctx$op.value("method", type = as.character, default = "fixed")
scale <- ctx$op.value("scale", type = as.integer, default = 5)

# =============================================================================
# Native implementation of flowVS algorithm for variance stabilization
# Based on: Azad et al. (2016) "flowVS: Channel-Specific Variance Stabilization
# in Flow Cytometry", BMC Bioinformatics
# =============================================================================

#' Find density peaks using kernel density estimation
#' @param x numeric vector of values
#' @param n number of points for density estimation
#' @return list with peak locations and boundaries
find_density_peaks <- function(x, n = 512) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)

  # Compute kernel density estimate
  dens <- density(x, n = n, na.rm = TRUE)

  # Find local maxima (peaks) and minima (valleys)
  # A peak is where the derivative changes from positive to negative
  dy <- diff(dens$y)
  peaks_idx <- which(dy[-length(dy)] > 0 & dy[-1] < 0) + 1
  valleys_idx <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1

  if (length(peaks_idx) == 0) return(NULL)

  # Filter peaks by prominence (at least 5% of max density)
  peak_heights <- dens$y[peaks_idx]
  min_height <- max(dens$y) * 0.05
  significant_peaks <- peaks_idx[peak_heights >= min_height]

  if (length(significant_peaks) == 0) return(NULL)

  # Get peak locations
  peak_locations <- dens$x[significant_peaks]

  # Determine boundaries between peaks using valleys
  boundaries <- c(min(x), dens$x[valleys_idx], max(x))
  boundaries <- sort(unique(boundaries))

  list(
    peaks = peak_locations,
    boundaries = boundaries,
    density = dens
  )
}

#' Assign data points to peaks based on boundaries
#' @param x numeric vector of values
#' @param peak_info output from find_density_peaks
#' @return integer vector of peak assignments (1, 2, 3, ...)
assign_to_peaks <- function(x, peak_info) {
  if (is.null(peak_info)) return(rep(1L, length(x)))

  peaks <- peak_info$peaks
  boundaries <- peak_info$boundaries

  # Assign each point to a peak region
  assignments <- rep(NA_integer_, length(x))

  for (i in seq_along(peaks)) {
    # Find boundary indices that surround this peak
    lower_bound <- max(boundaries[boundaries <= peaks[i]])
    upper_bound <- min(boundaries[boundaries >= peaks[i]])

    # Handle edge cases
    if (i == 1) lower_bound <- min(boundaries)
    if (i == length(peaks)) upper_bound <- max(boundaries)

    # Assign points within this range
    in_range <- x >= lower_bound & x < upper_bound
    assignments[in_range] <- i
  }

  # Assign remaining points to nearest peak
  remaining <- is.na(assignments)
  if (any(remaining)) {
    for (j in which(remaining)) {
      distances <- abs(x[j] - peaks)
      assignments[j] <- which.min(distances)
    }
  }

  assignments
}

#' Calculate Bartlett's test statistic for homogeneity of variances
#' @param x numeric vector of transformed values
#' @param groups integer vector of group assignments
#' @return Bartlett's statistic (lower = more homogeneous variances)
bartlett_statistic <- function(x, groups) {
  # Remove NA and infinite values
  valid <- is.finite(x) & !is.na(groups)
  x <- x[valid]
  groups <- groups[valid]

  unique_groups <- unique(groups)
  m <- length(unique_groups)  # number of groups

  if (m < 2) return(Inf)

  # Calculate per-group statistics
  n_i <- numeric(m)      # group sizes
  var_i <- numeric(m)    # group variances

  for (i in seq_along(unique_groups)) {
    group_data <- x[groups == unique_groups[i]]
    n_i[i] <- length(group_data)
    if (n_i[i] < 2) return(Inf)
    var_i[i] <- var(group_data)
  }

  # Check for zero or invalid variances
  if (any(var_i <= 0) || any(!is.finite(var_i))) return(Inf)

  # Total sample size
  n <- sum(n_i)

  # Pooled variance estimate
  df_i <- n_i - 1
  var_pooled <- sum(df_i * var_i) / sum(df_i)

  if (var_pooled <= 0 || !is.finite(var_pooled)) return(Inf)

  # Bartlett's statistic numerator
  numerator <- (n - m) * log(var_pooled) - sum(df_i * log(var_i))

  # Bartlett's correction factor
  correction <- 1 + (1 / (3 * (m - 1))) * (sum(1 / df_i) - 1 / (n - m))

  # Bartlett's statistic
  B <- numerator / correction

  if (!is.finite(B)) return(Inf)

  B
}

#' Objective function for flowVS optimization
#' @param cofactor the cofactor to evaluate
#' @param x raw data values
#' @return Bartlett's statistic after asinh transformation
flowvs_objective <- function(cofactor, x) {
  if (cofactor <= 0) return(Inf)

  # Transform data
  transformed <- asinh(x / cofactor)

  # Find peaks in transformed data
  peak_info <- find_density_peaks(transformed)

  if (is.null(peak_info) || length(peak_info$peaks) < 2) {
    # If we can't find multiple peaks, return a penalty
    # but not Inf to allow optimization to continue
    return(1e10)
  }

  # Assign points to peaks
  groups <- assign_to_peaks(transformed, peak_info)

  # Calculate Bartlett's statistic
  B <- bartlett_statistic(transformed, groups)

  B
}

#' Estimate optimal cofactor using flowVS algorithm
#' Performs golden section search over logarithmic cofactor space
#' @param values numeric vector of channel values
#' @param cf_low lower bound for cofactor search (log10 scale), default -1
#' @param cf_high upper bound for cofactor search (log10 scale), default 10
#' @return optimal cofactor value
estimate_cofactor_flowvs <- function(values, cf_low = -1, cf_high = 10) {
  values <- values[is.finite(values)]

  if (length(values) < 100) {
    # Not enough data for reliable peak detection, use fallback
    return(5)
  }

  # Grid search over logarithmic intervals
  # Divide [cf_low, cf_high] into unit intervals and search each
  intervals <- seq(cf_low, cf_high, by = 1)

  best_cofactor <- 5
  best_stat <- Inf

  for (i in seq_len(length(intervals) - 1)) {
    low <- 10^intervals[i]
    high <- 10^intervals[i + 1]

    # Golden section search within this interval
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

  # If optimization failed completely, use percentile fallback
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
# Main operator logic
# =============================================================================

if (method == "auto") {
  # Automatic cofactor estimation per channel using flowVS algorithm
  if (length(ctx$rnames) < 1)
    stop("auto method requires channel names in row projection")

  channel_colname <- ctx$rnames[[1]]

  # Get data with channel names
  data_df <- ctx %>%
    select(.ri, .ci, .y) %>%
    as_tibble()

  row_df <- ctx$rselect(channel_colname) %>%
    mutate(.ri = as.integer(row_number() - 1L))

  data_with_channels <- data_df %>%
    left_join(row_df, by = ".ri")

  # Estimate cofactors per channel using flowVS algorithm
  cofactor_df <- data_with_channels %>%
    group_by(!!sym(channel_colname)) %>%
    summarise(cofactor = estimate_cofactor_flowvs(.y), .groups = "drop")

  # Join cofactors back and apply transformation
  result <- data_with_channels %>%
    left_join(cofactor_df, by = channel_colname) %>%
    mutate(asinh = asinh(.y / cofactor)) %>%
    mutate(.ri = as.integer(.ri), .ci = as.integer(.ci)) %>%
    select(.ri, .ci, asinh)

  result %>%
    ctx$addNamespace() %>%
    ctx$save()

} else if (method == "manual") {
  # Manual per-channel scale values from row factor
  if (length(ctx$rnames) < 2)
    stop("manual method requires a scaling value after channel name in projection")

  row_df <- ctx$rselect(ctx$rnames[[2]]) %>%
    mutate(.ri = as.integer(row_number() - 1L)) %>%
    lazy_dt()

  scale_colname <- sym(ctx$rnames[[2]])

  ctx %>%
    select(.ri, .ci, .y) %>%
    lazy_dt() %>%
    left_join(row_df, by = ".ri") %>%
    mutate(asinh = asinh(.y / !!scale_colname)) %>%
    select(.ri, .ci, asinh) %>%
    as_tibble() %>%
    ctx$addNamespace() %>%
    ctx$save()

} else {
  # Fixed global scale (default)
  ctx %>%
    select(.ri, .ci, .y) %>%
    lazy_dt() %>%
    mutate(asinh = asinh(.y / scale)) %>%
    select(.ri, .ci, asinh) %>%
    as_tibble() %>%
    ctx$addNamespace() %>%
    ctx$save()
}

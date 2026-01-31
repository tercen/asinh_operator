library(tercenApi)
library(tercen)
library(data.table)
library(dplyr, warn.conflicts = FALSE)

ctx <- tercenCtx()
method <- ctx$op.value("method", type = as.character, default = "fixed")
scale <- ctx$op.value("scale", type = as.integer, default = 5)

# =============================================================================
# Standalone flowVS implementation for automatic cofactor estimation
# Based on: Azad et al. (2016) "flowVS: Channel-Specific Variance Stabilization
# in Flow Cytometry", BMC Bioinformatics
# This implementation does not require the flowVS/flowCore/flowStats packages.
# =============================================================================

# Source the standalone flowVS implementation
source("flowvs_standalone.R")

#' Estimate cofactors using standalone flowVS implementation
#' @param data_df data frame with columns: sample_id, and channel columns
#' @param channels character vector of channel names to estimate
#' @return data frame with columns: channel, cofactor
estimate_cofactors_flowvs <- function(data_df, channels) {
  # Ensure sample_id column exists
  if (!"sample_id" %in% names(data_df)) {
    stop("data_df must have a 'sample_id' column")
  }

  cat("Estimating cofactors using flowVS algorithm...\n")
  cat("Channels:", paste(channels, collapse = ", "), "\n")
  cat("Samples:", length(unique(data_df$sample_id)), "\n")

  # Use standalone implementation
  cofactors <- est_param_flowvs(data_df, channels, verbose = TRUE)

  result <- data.frame(channel = channels, cofactor = cofactors, stringsAsFactors = FALSE)

  cat("\nEstimated cofactors:\n")
  print(result)

  result
}

# =============================================================================
# Main operator logic
# =============================================================================

if (method == "auto") {
  # Automatic cofactor estimation per channel using flowVS
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

  # Get column names (the columns in ctx)
  col_df <- ctx$cselect() %>%
    mutate(.ci = as.integer(row_number() - 1L))

  # Check if there's a sample identifier in the column factors
  # Default to using .ci as sample_id if no sample column exists
  if (ncol(col_df) > 1) {
    sample_col <- names(col_df)[1]  # First column factor as sample_id
    data_with_samples <- data_with_channels %>%
      left_join(col_df, by = ".ci") %>%
      rename(sample_id = !!sym(sample_col))
  } else {
    # Use .ci as sample_id
    data_with_samples <- data_with_channels %>%
      mutate(sample_id = paste0("sample_", .ci))
  }

  # Get unique channels
  channels <- unique(data_with_samples[[channel_colname]])

  # Pivot data to wide format for flowVS (each channel as a column)
  # Using data.table::dcast instead of tidyr::pivot_wider
  pivot_data <- data_with_samples %>%
    select(sample_id, .ci, !!sym(channel_colname), .y) %>%
    mutate(cell_id = row_number()) %>%
    as.data.table()

  wide_data <- dcast(pivot_data, sample_id + cell_id ~ get(channel_colname),
                     value.var = ".y") %>%
    select(-cell_id)

  # Remove rows with any NA values (incomplete cases)
  wide_data <- wide_data[complete.cases(wide_data), ]

  # Estimate cofactors using flowVS
  cofactor_df <- estimate_cofactors_flowvs(as.data.frame(wide_data), channels)

  # Convert to named vector for lookup
  cofactors <- setNames(cofactor_df$cofactor, cofactor_df$channel)

  # Join cofactors back and apply transformation
  result <- data_with_channels %>%
    mutate(
      cofactor = cofactors[as.character(!!sym(channel_colname))],
      asinh = asinh(.y / cofactor)
    ) %>%
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
    mutate(.ri = as.integer(row_number() - 1L))

  scale_colname <- sym(ctx$rnames[[2]])

  ctx %>%
    select(.ri, .ci, .y) %>%
    as_tibble() %>%
    left_join(row_df, by = ".ri") %>%
    mutate(asinh = asinh(.y / !!scale_colname)) %>%
    mutate(.ri = as.integer(.ri), .ci = as.integer(.ci)) %>%
    select(.ri, .ci, asinh) %>%
    ctx$addNamespace() %>%
    ctx$save()

} else {
  # Fixed global scale (default)
  ctx %>%
    select(.ri, .ci, .y) %>%
    as_tibble() %>%
    mutate(asinh = asinh(.y / scale)) %>%
    mutate(.ri = as.integer(.ri), .ci = as.integer(.ci)) %>%
    select(.ri, .ci, asinh) %>%
    ctx$addNamespace() %>%
    ctx$save()
}

library(tercenApi)
library(tercen)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)

ctx <- tercenCtx()
method <- ctx$op.value("method", type = as.character, default = "fixed")
scale <- ctx$op.value("scale", type = as.integer, default = 5)

# Function to estimate optimal asinh cofactor for a channel
# Uses percentile-based estimation (5th percentile of positive values)
# This approach is commonly used in cytometry for variance stabilization
estimate_cofactor <- function(values, percentile = 0.05) {
  positive_values <- values[values > 0]
  if (length(positive_values) < 10) {
    # Not enough positive values, use default
    return(5)
  }
  cofactor <- quantile(positive_values, probs = percentile, na.rm = TRUE)
  # Ensure cofactor is reasonable (between 1 and 1000)
  cofactor <- max(1, min(1000, cofactor))
  return(as.numeric(cofactor))
}

if (method == "auto") {
  # Automatic cofactor estimation per channel
  if (length(ctx$rnames) < 1)
    stop("auto method requires channel names in row projection")

  channel_colname <- ctx$rnames[[1]]

  # Get data with channel names
  data_df <- ctx %>%
    select(.ri, .ci, .y) %>%
    as_tibble()

  row_df <- ctx$rselect(channel_colname) %>%
    mutate(.ri = row_number() - 1)

  data_with_channels <- data_df %>%
    left_join(row_df, by = ".ri")

  # Estimate cofactors per channel
  cofactor_df <- data_with_channels %>%
    group_by(!!sym(channel_colname)) %>%
    summarise(cofactor = estimate_cofactor(.y), .groups = "drop")

  # Join cofactors back and apply transformation
  result <- data_with_channels %>%
    left_join(cofactor_df, by = channel_colname) %>%
    mutate(asinh = asinh(.y / cofactor)) %>%
    select(.ri, .ci, asinh)

  result %>%
    ctx$addNamespace() %>%
    ctx$save()

} else if (method == "manual") {
  # Manual per-channel scale values from row factor
  if (length(ctx$rnames) < 2)
    stop("manual method requires a scaling value after channel name in projection")

  row_df <- ctx$rselect(ctx$rnames[[2]]) %>%
    mutate(.ri = row_number() - 1) %>%
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

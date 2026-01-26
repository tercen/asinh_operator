library(tercenApi)
library(tercen)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(flowCore)
library(flowVS)

ctx <- tercenCtx()
method <- ctx$op.value("method", type = as.character, default = "fixed")
scale <- ctx$op.value("scale", type = as.integer, default = 5)

if (method == "flowvs") {
  # flowVS automatic cofactor estimation
  if (length(ctx$rnames) < 1)
    stop("flowvs method requires channel names in row projection")

  channel_colname <- ctx$rnames[[1]]

  # Get data with channel names
  data_df <- ctx %>%
    select(.ri, .ci, .y) %>%
    as_tibble()

  row_df <- ctx$rselect(channel_colname) %>%
    mutate(.ri = row_number() - 1)

  data_with_channels <- data_df %>%
    left_join(row_df, by = ".ri")

  # Pivot to wide format (rows=events/cells, cols=channels)
  wide_data <- data_with_channels %>%
    pivot_wider(
      id_cols = .ci,
      names_from = !!sym(channel_colname),
      values_from = .y
    )

  # Get channel names (all columns except .ci)
  channels <- setdiff(colnames(wide_data), ".ci")

  # Create expression matrix
  expr_matrix <- as.matrix(wide_data[, channels])
  colnames(expr_matrix) <- channels

  # Create flowFrame
  ff <- flowFrame(exprs = expr_matrix)

  # Create flowSet (required by estParamFlowVS)
  fs <- flowSet(list(sample1 = ff))

  # Estimate cofactors per channel
  cofactors <- estParamFlowVS(fs, channels)

  # Create cofactor lookup table
  cofactor_df <- tibble(
    channel = channels,
    cofactor = cofactors
  )
  names(cofactor_df)[1] <- channel_colname

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

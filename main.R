library(tercenApi)
library(tercen)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)

ctx <- tercenCtx()
method <- ctx$op.value("method", type = as.character, default = "fixed")
scale <- ctx$op.value("scale", type = as.integer, default = 5)
docker_image <- ctx$op.value("docker.image", type = as.character, default = "ghcr.io/tercen/flowvs:latest")

# =============================================================================
# flowVS cofactor estimation using Docker container
# The container runs the actual flowVS package from Bioconductor
# Image: ghcr.io/tercen/flowvs (configurable via docker.image property)
# =============================================================================

#' Estimate cofactors using flowVS Docker container
#' @param data_df data frame with columns: sample_id, and channel columns
#' @param channels character vector of channel names to estimate
#' @param docker_image Docker image to use
#' @return data frame with columns: channel, cofactor
estimate_cofactors_docker <- function(data_df, channels, docker_image) {
  # Create temporary directory for input/output
  tmp_dir <- tempdir()
  input_file <- file.path(tmp_dir, "flowvs_input.csv")
  output_file <- file.path(tmp_dir, "flowvs_output.csv")
  script_file <- file.path(tmp_dir, "run_flowvs.R")

  # Ensure cleanup happens even on error
  on.exit({
    if (file.exists(input_file)) unlink(input_file)
    if (file.exists(output_file)) unlink(output_file)
    if (file.exists(script_file)) unlink(script_file)
  }, add = TRUE)

  # Ensure sample_id column exists
  if (!"sample_id" %in% names(data_df)) {
    stop("data_df must have a 'sample_id' column")
  }

  # Write input data
  write.csv(data_df, input_file, row.names = FALSE)

  # Build channels argument
  channels_arg <- paste(channels, collapse = ",")

  # Write R script to file (avoids escaping issues)
  r_script <- sprintf('
.libPaths(c("/app/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
pdf(NULL)
suppressPackageStartupMessages({
  library(flowVS)
  library(flowCore)
})

data <- read.csv("/data/flowvs_input.csv", stringsAsFactors = FALSE)
channels <- strsplit("%s", ",")[[1]]

sample_ids <- unique(data$sample_id)
frames <- lapply(sample_ids, function(sid) {
  subset_data <- data[data$sample_id == sid, channels, drop = FALSE]
  flowFrame(as.matrix(subset_data))
})
names(frames) <- sample_ids
fs <- flowSet(frames)

cofactors <- estParamFlowVS(fs, channels)

result <- data.frame(channel = channels, cofactor = cofactors)
write.csv(result, "/data/flowvs_output.csv", row.names = FALSE)
', channels_arg)

  writeLines(r_script, script_file)

  # Build Docker command using script file
  docker_cmd <- sprintf(
    'docker run --rm -w /tmp -v "%s:/data" --entrypoint /usr/local/bin/Rscript %s /data/run_flowvs.R',
    tmp_dir,
    docker_image
  )

  # Run Docker container with error handling
  cat("Running flowVS via Docker container...\n")
  cat("Docker image:", docker_image, "\n")
  cat("Channels:", channels_arg, "\n")

  result <- tryCatch({
    output <- system(docker_cmd, intern = TRUE, ignore.stderr = FALSE)
    attr(output, "status")
  }, error = function(e) {
    stop("Failed to run Docker container: ", e$message, "\n",
         "Ensure the Docker image '", docker_image, "' is available. ",
         "You can pull it with: docker pull ", docker_image)
  })

  exit_code <- if (is.null(result)) 0 else result

  if (exit_code != 0) {
    stop("Docker flowVS container failed with exit code: ", exit_code, "\n",
         "Ensure the Docker image '", docker_image, "' is available and working.")
  }

  # Read output
  if (!file.exists(output_file)) {
    stop("flowVS output file not found. The Docker container may have failed silently. ",
         "Check that the image '", docker_image, "' contains the flowVS package.")
  }

  cofactor_df <- read.csv(output_file, stringsAsFactors = FALSE)

  cat("Estimated cofactors:\n")
  print(cofactor_df)

  cofactor_df
}

# =============================================================================
# Main operator logic
# =============================================================================

if (method == "auto") {
  # Automatic cofactor estimation per channel using flowVS Docker container
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
  wide_data <- data_with_samples %>%
    select(sample_id, .ci, !!sym(channel_colname), .y) %>%
    mutate(cell_id = row_number()) %>%
    pivot_wider(
      id_cols = c(sample_id, cell_id),
      names_from = !!sym(channel_colname),
      values_from = .y
    ) %>%
    select(-cell_id)

  # Remove rows with any NA values (incomplete cases)
  wide_data <- wide_data[complete.cases(wide_data), ]

  # Estimate cofactors using Docker
  cofactor_df <- estimate_cofactors_docker(wide_data, channels, docker_image)

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

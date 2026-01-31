# =============================================================================
# Standalone implementation of flowVS estParamFlowVS algorithm
# Based on: Azad et al. (2016) "flowVS: Channel-Specific Variance Stabilization
# in Flow Cytometry", BMC Bioinformatics
#
# This implementation replicates the exact behavior of flowVS::estParamFlowVS
# without requiring the full flowVS/flowCore/flowStats package stack.
# =============================================================================

#' Compute bandwidth using flowStats formula
#' @param y numeric vector
#' @param bwFac bandwidth factor (default 2)
#' @return bandwidth value
compute_bandwidth <- function(y, bwFac = 2) {
  n <- length(y)
  st.dev <- sqrt(var(y, na.rm = TRUE))
  Q1.val <- quantile(y, 0.25, na.rm = TRUE)
  Q3.val <- quantile(y, 0.75, na.rm = TRUE)
  IQR.val <- (Q3.val - Q1.val) / (qnorm(0.75) - qnorm(0.25))
  bwNS <- min(st.dev, IQR.val) * (4 / (7 * n))^(1/9)
  bwFac * bwNS
}

#' Find peaks and valleys in density estimate
#' This replicates the peak detection behavior of flowVS
#'
#' Uses an adaptive bandwidth approach: starts with a smaller bandwidth to find
#' candidate peaks, then uses R's default density() for finding valleys.
#' This better approximates the behavior of flowVS which uses curv1Filter for
#' peak detection and default density for valley/separator detection.
#'
#' @param y numeric vector of (transformed) values
#' @param bwFac bandwidth factor (default 2, as in flowVS)
#' @param gridsize grid size for density estimation (default 512)
#' @param min_prominence minimum prominence for a peak to be considered (relative to max density)
#' @return list with:
#'   - peaks: matrix with columns x, y for peak locations
#'   - valleys: vector of valley x locations between peaks
#'   - density: the default density estimate (used for valley detection)
find_density_peaks_valleys <- function(y, bwFac = 2, gridsize = 512, min_prominence = 0.05) {
  n <- length(y)
  if (n < 10) return(list(peaks = matrix(nrow = 0, ncol = 2), valleys = numeric(0), density = NULL))

  # Compute base bandwidth using flowStats formula
  h_base <- compute_bandwidth(y, bwFac)

  # Use a smaller bandwidth for peak detection to avoid over-smoothing
  # This better approximates featureSignif's behavior which detects more subtle peaks
  h_detect <- h_base * 0.7

  # Compute kernel density estimate with detection bandwidth for peak finding
  dens_detect <- density(y, bw = h_detect, n = gridsize)

  # Find local maxima (peaks)
  d1 <- diff(dens_detect$y)
  peaks_idx <- which(d1[-length(d1)] > 0 & d1[-1] <= 0) + 1

  if (length(peaks_idx) == 0) {
    peaks_idx <- which.max(dens_detect$y)
  }

  # Filter peaks by prominence (at least min_prominence of max density)
  max_dens <- max(dens_detect$y)
  prominent <- dens_detect$y[peaks_idx] >= max_dens * min_prominence
  peaks_idx <- peaks_idx[prominent]

  if (length(peaks_idx) == 0) {
    peaks_idx <- which.max(dens_detect$y)
  }

  peaks <- cbind(x = dens_detect$x[peaks_idx], y = dens_detect$y[peaks_idx])

  # Use default density for valley detection (as flowVS does)
  # flowVS's densityPeaks uses: dens <- density(y) without bandwidth specification
  dens_default <- density(y)

  # Find valleys (local minima) between consecutive peaks using default density
  valleys <- numeric(0)
  if (nrow(peaks) > 1) {
    for (i in seq_len(nrow(peaks) - 1)) {
      # Find minimum between peak i and peak i+1
      sel <- dens_default$x > peaks[i, 1] & dens_default$x < peaks[i + 1, 1]
      if (any(sel)) {
        valley_idx <- which(sel)[which.min(dens_default$y[sel])]
        valleys <- c(valleys, dens_default$x[valley_idx])
      }
    }
  }

  list(peaks = peaks, valleys = valleys, density = dens_default)
}

#' Find peaks in density and compute population statistics
#' This replicates flowVS::densityPeaks
#'
#' @param y numeric vector of (transformed) values
#' @param bwFac bandwidth factor (default 2)
#' @param populationQuant quantile for trimming populations (default 0.01)
#' @param borderQuant quantile for trimming borders (default 0.001)
#' @param peak.distance.thr threshold for merging close peaks (default 0.05)
#' @return matrix with columns: n, mean, median, variance
density_peaks <- function(y, bwFac = 2, populationQuant = 0.01,
                          borderQuant = 0.001, peak.distance.thr = 0.05) {
  # Trim border values (same as flowVS)
  y <- y[y > min(y) & y < max(y)]
  min.range <- quantile(y, borderQuant)
  max.range <- quantile(y, 1 - borderQuant)
  y <- y[y > min.range & y < max.range]

  if (length(y) < 10) {
    return(matrix(nrow = 0, ncol = 4,
                  dimnames = list(NULL, c("n", "mean", "median", "variance"))))
  }

  # Find peaks and valleys using density estimation
  pv <- find_density_peaks_valleys(y, bwFac = bwFac)
  peaks <- pv$peaks
  valleys <- pv$valleys
  dens <- pv$density

  if (nrow(peaks) == 0) {
    return(matrix(nrow = 0, ncol = 4,
                  dimnames = list(NULL, c("n", "mean", "median", "variance"))))
  }

  # Merge peaks that are too close (same logic as flowVS)
  mrange <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
  npeaks <- nrow(peaks)
  peak.rm.idx <- c()

  if (npeaks >= 2) {
    p1 <- 1
    p2 <- 2
    while (p2 <= npeaks) {
      peak1 <- peaks[p1, ]
      peak2 <- peaks[p2, ]
      if (abs(peak1[1] - peak2[1]) < mrange * peak.distance.thr) {
        # Remove the smaller peak
        if (peak1[2] < peak2[2]) {
          peak.rm.idx <- c(peak.rm.idx, p1)
          p1 <- p2
        } else {
          peak.rm.idx <- c(peak.rm.idx, p2)
        }
        p2 <- p2 + 1
      } else {
        p1 <- p2
        p2 <- p2 + 1
      }
    }
  }

  if (length(peak.rm.idx) > 0) {
    peaks <- peaks[-peak.rm.idx, , drop = FALSE]
    # Also remove corresponding valleys
    if (length(valleys) > 0) {
      # Recalculate valleys after removing peaks
      if (nrow(peaks) > 1) {
        new_valleys <- numeric(0)
        for (i in seq_len(nrow(peaks) - 1)) {
          sel <- dens$x > peaks[i, 1] & dens$x < peaks[i + 1, 1]
          if (any(sel)) {
            valley_idx <- which(sel)[which.min(dens$y[sel])]
            new_valleys <- c(new_valleys, dens$x[valley_idx])
          }
        }
        valleys <- new_valleys
      } else {
        valleys <- numeric(0)
      }
    }
  }

  if (nrow(peaks) == 0) {
    return(matrix(nrow = 0, ncol = 4,
                  dimnames = list(NULL, c("n", "mean", "median", "variance"))))
  }

  # Compute population statistics for each peak
  # Use valleys between peaks as separators (same as flowVS)
  peaksStats <- matrix(nrow = 0, ncol = 4)
  colnames(peaksStats) <- c("n", "mean", "median", "variance")

  peaks <- peaks[order(peaks[, 1]), , drop = FALSE]  # Sort by x position
  npeaks <- nrow(peaks)

  left <- min(y, na.rm = TRUE)
  for (i in seq_len(npeaks)) {
    p <- peaks[i, ]

    if (i == npeaks) {
      right <- max(y, na.rm = TRUE)
    } else if (length(valleys) >= i) {
      right <- valleys[i]
    } else {
      # Find minimum density between this peak and next
      sel <- dens$x > peaks[i, 1] & dens$x < peaks[i + 1, 1]
      if (any(sel)) {
        right <- dens$x[sel][which.min(dens$y[sel])]
      } else {
        right <- (peaks[i, 1] + peaks[i + 1, 1]) / 2
      }
    }

    # Select population
    y.sel <- y[y >= left & y <= right]

    if (length(y.sel) > 0) {
      # Trim population by quantiles
      min.range.pop <- quantile(y.sel, populationQuant)
      max.range.pop <- quantile(y.sel, 1 - populationQuant)

      # Handle asymmetric peaks (same logic as flowVS)
      if (!is.na(p[1])) {
        r <- (max.range.pop - p[1]) / (p[1] - min.range.pop + 1e-10)
        if (r >= 3 && i == npeaks) {
          right <- min(max.range.pop, p[1] + 3 * (p[1] - min.range.pop))
          y.sel <- y[y >= left & y <= right]
          if (length(y.sel) > 0) {
            min.range.pop <- quantile(y.sel, populationQuant)
            max.range.pop <- quantile(y.sel, 1 - populationQuant)
          }
        } else if (r <= 1/3 && i == 1) {
          left <- max(min.range.pop, p[1] - 3 * (max.range.pop - p[1]))
          y.sel <- y[y >= left & y <= right]
          if (length(y.sel) > 0) {
            min.range.pop <- quantile(y.sel, populationQuant)
            max.range.pop <- quantile(y.sel, 1 - populationQuant)
          }
        }
      }

      # Final population
      population <- y.sel[y.sel >= min.range.pop & y.sel <= max.range.pop]

      if (length(population) > 1) {
        peaksStats <- rbind(peaksStats, c(
          length(population),
          mean(population),
          median(population),
          var(population)
        ))
      }
    }

    left <- right
  }

  peaksStats
}

#' Bartlett's test for homogeneity of variances
#' This is an exact copy of flowVS::bartlettTest
#'
#' @param peakStats matrix with columns n, mean, median, variance, [sample]
#' @return Bartlett's test statistic
bartlett_test <- function(peakStats) {
  # Remove peaks with only 1 observation
  peakStats <- peakStats[peakStats[, 1] > 1, , drop = FALSE]

  if (nrow(peakStats) <= 1) return(1e9)

  size <- peakStats[, 1]
  variance <- peakStats[, 4]

  # Check for invalid variances
  if (any(variance <= 0) || any(!is.finite(variance))) return(1e9)

  sum1 <- sum((size - 1) * log(variance))
  sum2 <- sum((size - 1) * variance)
  N <- sum(size)
  sum3 <- sum(1 / (size - 1))
  k <- nrow(peakStats)

  # Bartlett's statistic
  bt <- ((N - k) * log(sum2 / (N - k)) - sum1) /
    (1 + (sum3 - 1 / (N - k)) / (3 * (k - 1)))

  if (!is.finite(bt)) return(1e9)

  bt
}

#' Compute flowVS objective (Bartlett statistic) for a given cofactor
#' This replicates flowVS::flowVS1D
#'
#' @param cofactor the cofactor to evaluate
#' @param data_list list of numeric vectors, one per sample
#' @param bwFac bandwidth factor (default 2)
#' @param populationQuant quantile for trimming populations (default 0.01)
#' @param borderQuant quantile for trimming borders (default 0.001)
#' @return Bartlett's statistic (lower = better variance stabilization)
flowvs_objective <- function(cofactor, data_list, bwFac = 2,
                             populationQuant = 0.01, borderQuant = 0.001) {
  MAX_BT <- 1e9

  peakStats <- matrix(nrow = 0, ncol = 5)
  colnames(peakStats) <- c("n", "mean", "median", "variance", "sample")

  for (i in seq_along(data_list)) {
    y <- data_list[[i]]
    y <- asinh(y / cofactor)

    peakStat <- density_peaks(y, bwFac = bwFac,
                              populationQuant = populationQuant,
                              borderQuant = borderQuant)

    if (nrow(peakStat) > 0) {
      peakStat <- cbind(peakStat, rep(i, nrow(peakStat)))
      peakStats <- rbind(peakStats, peakStat)
    }
  }

  if (nrow(peakStats) <= 1) return(MAX_BT)

  # Remove outlier samples (10%) - same as flowVS
  sampleIdx <- peakStats[, 5]
  nsamples <- length(unique(sampleIdx))
  npeaks_per_sample <- table(sampleIdx)
  sampleOrder <- order(abs(npeaks_per_sample - median(npeaks_per_sample)),
                       decreasing = TRUE)
  rm <- round(nsamples / 10)

  if (rm > 0 && nsamples > (rm + 1)) {
    rmIdx <- sampleOrder[1:rm]
    rmSamples <- as.integer(names(npeaks_per_sample[rmIdx]))
    selected <- !(sampleIdx %in% rmSamples)
    peakStats <- peakStats[selected, , drop = FALSE]
  }

  if (nrow(peakStats) <= 1) return(MAX_BT)

  # Remove outlier peaks (10%) - same as flowVS
  pkVar <- peakStats[, 4]
  npeaks <- length(pkVar)
  varOrder <- order(abs(pkVar - median(pkVar)), decreasing = TRUE)
  rm <- round(npeaks / 10)

  if (rm > 0 && npeaks > (rm + 1)) {
    rmIdx <- varOrder[1:rm]
    peakStats <- peakStats[-rmIdx, , drop = FALSE]
  }

  if (nrow(peakStats) <= 1) return(MAX_BT)

  bartlett_test(peakStats)
}

#' Find optimal cofactor for a single channel
#' This replicates flowVS::optimStat
#'
#' @param data_list list of numeric vectors, one per sample
#' @param cfLow lower bound for cofactor search (log10 scale), default -1
#' @param cfHigh upper bound for cofactor search (log10 scale), default 10
#' @param verbose print progress (default TRUE)
#' @return optimal cofactor value
optim_cofactor <- function(data_list, cfLow = -1, cfHigh = 10, verbose = TRUE) {
  MAX_BT <- 1e9

  if (cfLow >= cfHigh) {
    cfLow <- -1
    cfHigh <- 10
  }

  cf <- cfLow:cfHigh
  ncf <- length(cf)
  cfopt <- rep(0, ncf - 1)
  btopt <- rep(0, ncf - 1)

  if (verbose) {
    cat(sprintf("%18s       %10s    %15s\n", "cf range", "opt cf", "Bartlett's stat"))
    cat("==========================================================\n")
  }

  # Grid search over logarithmic intervals
  for (i in seq_len(ncf - 1)) {
    low <- exp(cf[i])
    high <- exp(cf[i + 1])
    tol <- (high - low) / 10

    # Optimize within this interval
    opt <- suppressWarnings(
      optimize(f = flowvs_objective, interval = c(low, high),
               data_list = data_list, tol = tol)
    )

    btopt[i] <- opt$objective
    cfopt[i] <- opt$minimum

    if (verbose) {
      cat(sprintf("[%9.2f, %9.2f ] %10.2f ", low, high, opt$minimum))
      if (opt$objective == MAX_BT) {
        cat(sprintf("%10s\n", "MAX (10^9)"))
      } else {
        cat(sprintf("%15.2f\n", opt$objective))
      }
    }
  }

  # Local refinement around the best interval
  minIdx <- which.min(btopt)
  del <- cfopt[minIdx] / 10
  btLocal <- rep(0, 11)
  btLocal[6] <- btopt[minIdx]
  cfLocal <- c(5:1, cfopt[minIdx], 1:5)
  cfLocal[1:5] <- cfopt[minIdx] - 5:1 * del
  cfLocal[7:11] <- cfopt[minIdx] + 1:5 * del

  for (i in c(1:5, 7:11)) {
    btLocal[i] <- flowvs_objective(cfLocal[i], data_list)
  }

  minIdx <- which.min(btLocal)
  cfLocal[minIdx]
}

#' Estimate optimal cofactors for all channels
#' This replicates flowVS::estParamFlowVS
#'
#' @param data data.frame with sample_id column and channel columns
#' @param channels character vector of channel names to process
#' @param verbose print progress (default TRUE)
#' @return named numeric vector of optimal cofactors
est_param_flowvs <- function(data, channels, verbose = TRUE) {
  # Validate inputs
  if (!"sample_id" %in% names(data)) {
    stop("data must contain a 'sample_id' column")
  }

  missing_channels <- setdiff(channels, names(data))
  if (length(missing_channels) > 0) {
    stop("Channels not found in data: ", paste(missing_channels, collapse = ", "))
  }

  samples <- unique(data$sample_id)
  cofactors <- numeric(length(channels))
  names(cofactors) <- channels

  for (col in channels) {
    if (verbose) {
      cat("====================================================================\n")
      cat("Channel ", col, " : Finding optimum cofactor for asinh transformation\n")
      cat("====================================================================\n")
    }

    # Create list of sample data
    data_list <- lapply(samples, function(s) {
      data[data$sample_id == s, col]
    })

    cf <- optim_cofactor(data_list, verbose = verbose)
    cofactors[col] <- cf

    if (verbose) {
      cat("\nOptimum cofactor for ", col, " : ", cf, "\n")
      cat("====================================================================\n\n")
    }
  }

  cofactors
}

# Export for use in main.R
if (FALSE) {
  # Test code - run this to validate against flowVS output
  data <- read.csv("flowvs_input.csv", stringsAsFactors = FALSE)
  channels <- c("CD4", "CD8", "CD3")
  cofactors <- est_param_flowvs(data, channels)
  print(cofactors)

  # Expected output:
  # CD4: 6450.34368035551
  # CD8: 4766.69839674341
  # CD3: 6317.34353589174
}

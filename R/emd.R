#' Earth Mover's Distance between spectral histograms
#'
#' @description pcv.emd can be used to calculate Earth Mover's Distance between pairwise histograms
#' in a wide dataframe of multi value traits. The is expected to be used with output from \code{mv_ag}.
#' See also \link{pcv.euc} for euclidean distance between histograms.
#'
#' @param df Data frame to use with multi value traits in wide format or long format
#' @param cols Columns to use. Defaults to NULL in which case all columns are used.
#' Single strings will be used to regex a pattern in column names (see examples).
#'  A vector of names, positions, or booleans will also work.
#'  For long data this is taken as a regex pattern (or full name)
#'  to use in filtering the trait column.
#' @param reorder Should data be reordered to put similar rows together in the resulting plot?
#'  This takes a vector of column names of length 1 or more (see examples).
#' @param include if a long dataframe is returned then these columns will be added to the dataframe,
#'  labelled for i and j (the row positions for compared histograms).
#'  If a matrix is returned then this information is stored in the row names.
#'  This defaults to reorder.
#' @param mat Logical, should data be returned as an nrow x nrow matrix or as a long dataframe?
#'  By Default this is FALSE and a long dataframe is returned.
#'  Both options are comparable in terms of speed,
#'  although for large datasets the matrix version may be slightly faster.
#' @param plot Logical, should a plot be returned? For a matrix this is made with heatmap(),
#'  for a dataframe this uses ggplot.
#' @param parallel Number of cores to use. Defaults to 1 unless the "mc.cores" option is set.
#' @param trait Column name for long data to identify traits. This defaults to "trait". If this and
#' value are in the column names of the data then it is assumed to be in long format,
#' otherwise it is assumed to be in wide format.
#' @param id A vector of column names that uniquely identifies observations if the
#' data is in long format. Defaults to "image".
#' @param value A column name for the values to be drawn from in long data.
#' Defaults to "value".
#' @param raiseError Logical, should warnings/errors be raised for potentially large output?
#'  It is easy to ask for very many comparisons with this function so the goal of this argument
#'  is to catch a few of those and give estimates of how much time something may take.
#'   If the function is expected to take very long then a warning or an error is raised.
#'    If this is set to FALSE then no time estimates are made.
#' @param method Which method to use (one of "emd" or "euc"). Defaults to "emd".
#' @import ggplot2
#' @import parallel
#' @return A dataframe/matrix (if plot=FALSE) or a list with a dataframe/matrix and\
#' a ggplot (if plot=TRUE). The returned data contains pairwise EMD values.
#'
#' @keywords emd earth-mover's-distance multi-value histogram
#' @examples
#'
#'
#' set.seed(123)
#' test <- mvSim(
#'   dists = list(
#'     runif = list(min = 0, max = 100),
#'     rnorm = list(mean = 90, sd = 20)
#'   ),
#'   n_samples = 10
#' )
#' test$meta1 <- rep(LETTERS[1:3], length.out = nrow(test))
#' test$meta2 <- rep(LETTERS[4:5], length.out = nrow(test))
#'
#' x <- pcv.emd(
#'   df = test, cols = "sim", reorder = "group",
#'   include = c("meta1", "meta2"), mat = FALSE,
#'   plot = FALSE, parallel = 1
#' )
#' head(x)
#' x2 <- pcv.emd(
#'   df = test, cols = "sim", reorder = "group",
#'   include = c("meta1", "meta2"), mat = FALSE,
#'   plot = FALSE, parallel = 1, method = "euc"
#' )
#' head(x2)
#'
#' \donttest{
#' tryCatch(
#'   {
#'     library(data.table)
#'     file <- paste0(
#'       "https://media.githubusercontent.com/media/joshqsumner/",
#'       "pcvrTestData/main/pcv4-multi-value-traits.csv"
#'     )
#'     df1 <- read.pcv(file, "wide", reader = "fread")
#'
#'     df1$genotype <- substr(df1$barcode, 3, 5)
#'     df1$genotype <- ifelse(df1$genotype == "002", "B73",
#'       ifelse(df1$genotype == "003", "W605S",
#'         ifelse(df1$genotype == "004", "MM", "Mo17")
#'       )
#'     )
#'     df1$fertilizer <- substr(df1$barcode, 8, 8)
#'     df1$fertilizer <- ifelse(df1$fertilizer == "A", "100",
#'       ifelse(df1$fertilizer == "B", "50", "0")
#'     )
#'
#'     w <- pcv.emd(df1,
#'       cols = "hue_frequencies", reorder = c("fertilizer", "genotype"),
#'       mat = FALSE, plot = TRUE, parallel = 1
#'     )
#'   },
#'   error = function(err) {
#'     message(err)
#'   }
#' )
#'
#' # Note on computational complexity
#' # This scales as O^2, see the plot below for some idea
#' # of the time for different input data sizes.
#' emdTime <- function(x, n = 1) {
#'   x^2 / n * 0.0023
#' }
#' plot(
#'   x = c(18, 36, 54, 72, 108, 135), y = c(0.74, 2.89, 6.86, 10.99, 26.25, 42.44),
#'   xlab = "N Input Images", ylab = "time (seconds)"
#' ) # benchmarked test data
#' lines(x = 1:150, y = emdTime(1:150)) # exponential function
#'
#' plot(
#'   x = 1:1000, y = emdTime(1:1000), type = "l",
#'   xlab = "N Input Images", ylab = "time (seconds)"
#' )
#' }
#'
#' @export
#'
pcv.emd <- function(df, cols = NULL, reorder = NULL, include = reorder, mat = FALSE, plot = TRUE,
                    parallel = getOption("mc.cores", 1), trait = "trait", id = "image",
                    value = "value", raiseError = TRUE, method = "emd") {
  if (all(c(trait, value) %in% colnames(df))) {
    long <- TRUE
    traitCol <- trait
  } else {
    long <- FALSE
  }
  if (!is.null(reorder)) {
    df <- df[order(interaction(df[, reorder])), ]
  }
  if (long) {
    out_data <- .longEmdCalculation(
      df, cols, traitCol, raiseError,
      parallel, id, value, include, mat, method
    )
  } else { # wide input

    out_data <- .wideEmdCalculation(df, cols, raiseError, parallel, id, include, mat, method)
  }
  if (plot) {
    if (mat) {
      p <- stats::heatmap(out_data)
    } else {
      p <- ggplot2::ggplot(out_data, ggplot2::aes(x = .data$i, y = .data$j, fill = .data[[method]])) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::labs(fill = method) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.line.x.bottom = ggplot2::element_line(), axis.line.y.left = ggplot2::element_line(),
          legend.position = "bottom"
        )
    }
    outList <- list("data" = out_data, "plot" = p)
  } else {
    outList <- out_data
  }
  return(outList)
}


#' Error Raising function for very long EMD calculation times
#' @keywords internal
#' @noRd

.emdRaiseError <- function(raiseError, df, parallel, trait = NULL) {
  n <- ifelse(is.null(trait), nrow(df), nrow(df) / length(unique(df[[trait]])))
  if (raiseError) {
    et_sec <- 0.00125 * ((n / parallel)^2)
    et_min <- et_sec / 60
    et_hour <- et_min / 60
    if (et_sec <= 300) {
      message(paste0(
        "Estimated time of calculation is roughly ", round(et_sec, 1),
        " seconds using ", parallel, " cores in parallel."
      ))
    } else if (et_min < 60) {
      warning(paste0(
        "Estimated time of calculation is roughly ", round(et_min, 2),
        " minutes using ", parallel, " cores in parallel."
      ))
    } else if (et_min > 60) {
      stop(paste0(
        "Stopping, estimated time of calculation is roughly ", round(et_hour, 2),
        " hours using ", parallel, " cores in parallel.",
        "\nIf you wish to proceed then rerun this command with raiseError=FALSE"
      ))
    }
  }
  return(invisible(NULL))
}

#' Long EMD calculation
#' @keywords internal
#' @noRd

.longEmdCalculation <- function(df, cols, traitCol, raiseError, parallel, id,
                                value, include, mat, method) {
  dist_1d <- match.fun(paste0(method, "1d"))
  df <- df[grepl(cols, df[[traitCol]]), ]
  #* `raise error`
  .emdRaiseError(raiseError, df, parallel, trait = traitCol)
  #* `calculate emd`
  df$INNER_ID_EMD <- interaction(df[, id], drop = TRUE)
  if (mat) { # make dist matrix
    mat_obj <- matrix(0,
      nrow = length(unique(df$INNER_ID_EMD)),
      ncol = length(unique(df$INNER_ID_EMD))
    )
    values <- unlist(parallel::mclapply(seq_along(unique(df$INNER_ID_EMD)), function(i_n) {
      i_res <- lapply(seq_along(unique(df$INNER_ID_EMD)), function(j_n) {
        i <- unique(df$INNER_ID_EMD)[i_n]
        j <- unique(df$INNER_ID_EMD)[j_n]
        if (i_n < j_n) {
          d <- dist_1d(
            as.numeric(df[df$INNER_ID_EMD == as.character(i), value]),
            as.numeric(df[df$INNER_ID_EMD == as.character(j), value])
          )
          return(d)
        }
      })
      return(i_res)
    }, mc.cores = parallel))
    mat_obj[lower.tri(mat_obj)] <- values
    tmat_obj <- t(mat_obj)
    mat_obj[upper.tri(mat_obj)] <- tmat_obj[upper.tri(tmat_obj)]
    rownames(mat_obj) <- colnames(mat_obj) <- unique(df$INNER_ID_EMD)
    out_data <- mat_obj
  } else { # make long data
    out_data <- do.call(rbind, parallel::mclapply(seq_along(unique(df$INNER_ID_EMD)), function(i_n) {
      i_df <- do.call(rbind, lapply(seq_along(unique(df$INNER_ID_EMD)), function(j_n) {
        i <- unique(df$INNER_ID_EMD)[i_n]
        j <- unique(df$INNER_ID_EMD)[j_n]
        emdOut <- NULL
        if (i_n == j_n) {
          emdOut <- 0
        } else if (i_n < j_n) {
          emdOut <- dist_1d(
            as.numeric(df[df$INNER_ID_EMD == as.character(i), value]),
            as.numeric(df[df$INNER_ID_EMD == as.character(j), value])
          )
        }
        if (!is.null(emdOut)) {
          if (!is.null(include)) {
            x <- rbind(
              data.frame(
                i = i, j = j, emd = emdOut,
                df[df$INNER_ID_EMD == as.character(i), include][1, ],
                df[df$INNER_ID_EMD == as.character(j), include][1, ]
              ),
              data.frame(
                i = j, j = i, emd = emdOut,
                df[df$INNER_ID_EMD == as.character(j), include][1, ],
                df[df$INNER_ID_EMD == as.character(i), include][1, ]
              )
            )
            colnames(x) <- c("i", "j", method, paste0(include, "_i"), paste0(include, "_j"))
          } else {
            x <- data.frame(i = c(i, j), j = c(j, i), emd = emdOut)
          }
          return(x)
        }
      }))
      return(i_df)
    }, mc.cores = parallel))
  }
  return(out_data)
}

#* wide EMD calculation
#' @keywords internal
#' @noRd

.wideEmdCalculation <- function(df, cols, raiseError, parallel, id, include, mat, method) {
  dist_1d <- match.fun(paste0(method, "1d"))
  if (is.null(cols)) {
    cols <- colnames(df)
  } else if (is.character(cols) && length(cols) == 1) {
    cols <- grepl(cols, colnames(df))
  }
  #* `raise error`
  .emdRaiseError(raiseError, df, parallel, trait = NULL)
  #* `calculate emd`
  if (mat) { # make dist matrix
    mat_obj <- matrix(0, nrow = nrow(df), ncol = nrow(df))
    values <- unlist(parallel::mclapply(seq_len(nrow(df)), function(i) {
      i_res <- lapply(seq_len(nrow(df)), function(j) {
        if (i < j) {
          return(dist_1d(as.numeric(df[i, cols]), as.numeric(df[j, cols])))
        }
      })
      return(i_res)
    }, mc.cores = parallel))
    mat_obj[lower.tri(mat_obj)] <- values
    tmat_obj <- t(mat_obj)
    mat_obj[upper.tri(mat_obj)] <- tmat_obj[upper.tri(tmat_obj)]
    rownames(mat_obj) <- colnames(mat_obj) <- seq_len(nrow(df))
    out_data <- mat_obj
    if (!is.null(include)) {
      rownames(out_data) <- interaction(df[, include])
    }
  } else { # make long dataframe
    out_data <- do.call(rbind, parallel::mclapply(seq_len(nrow(df)), function(i) {
      i_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(j) {
        emdOut <- NULL
        if (i == j) {
          emdOut <- 0
        } else if (i < j) {
          emdOut <- dist_1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))
        }
        if (!is.null(emdOut)) {
          if (!is.null(include)) {
            x <- rbind(
              data.frame(i = i, j = j, emd = emdOut, df[i, include], df[j, include]),
              data.frame(i = j, j = i, emd = emdOut, df[j, include], df[i, include])
            )
            colnames(x) <- c("i", "j", method, paste0(include, "_i"), paste0(include, "_j"))
          } else {
            x <- data.frame(i = c(i, j), j = c(j, i), emd = emdOut)
          }
          return(x)
        }
      }))
      return(i_df)
    }, mc.cores = parallel))
  }
  return(out_data)
}

#' Earth Mover's Distance between spectral histograms
#'
#' @description emd1d computes 1 dimension EMD for two samples.
#'
#' @param s1 Histogram as a numeric vector of counts per position.
#' @param s2 Histogram as a numeric vector of counts per position. Must be the same length as s1.
#'
#' @keywords internal
#' @return Returns EMD between two samples as a numeric.
#' @examples
#'
#' set.seed(123)
#' s1 <- hist(rnorm(10000, 50, 10), breaks = seq(1, 100, 1))$counts
#' s2 <- hist(rlnorm(9000, log(30), 0.25), breaks = seq(1, 100, 1))$counts
#' plot(s2, type = "l"); lines(s1)
#' emd1d(s1, s2)
#'
#' @noRd
#'
emd1d <- function(s1, s2) {
  if (length(s1) != length(s2)) {
    stop("Samples must be from the same histogram and be of the same length")
  }
  s1 <- s1 / sum(s1)
  s2 <- s2 / sum(s2)
  emd_iter <- numeric(length(s1) + 1)
  for (i in seq_along(s1)) {
    emd_iter[i + 1] <- s1[i] - s2[i] + emd_iter[i]
  }
  return(sum(abs(emd_iter)))
}

#' Euclidean Distance between spectral histograms
#'
#' @description euc1d computes euclidean distance between two samples.
#'
#' @param s1 Histogram as a numeric vector of counts per position.
#' @param s2 Histogram as a numeric vector of counts per position. Must be the same length as s1.
#'
#' @importFrom stats dist
#' @keywords internal
#' @return Returns euclidean distance as a numeric
#'
#' @noRd

euc1d <- function(s1, s2) {
  if (length(s1) != length(s2)) {
    stop("Samples must be from the same histogram and be of the same length")
  }
  s1 <- s1 / sum(s1)
  s2 <- s2 / sum(s2)
  mat <- matrix(c(s1, s2), nrow = 2, byrow = TRUE)
  euc <- as.numeric(stats::dist(mat, method = "euclidean"))
  return(euc)
}

#' @rdname pcv.emd
#' @export
#'
pcv.euc <- function(df, cols = NULL, reorder = NULL, include = reorder, mat = FALSE, plot = TRUE,
                    parallel = getOption("mc.cores", 1), trait = "trait", id = "image",
                    value = "value", raiseError = TRUE, method = "euc") {
  out <- pcv.emd(df, cols, reorder, include, mat, plot, parallel, trait, id, value, raiseError, method)
  return(out)
}

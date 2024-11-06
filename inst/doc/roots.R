## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(pcvr)
library(ggplot2)
theme_set(pcv_theme())

## -----------------------------------------------------------------------------
rRhyzoDist <- function(n, theta = 0.3, u1_max = 20, u2_max = 5500, sd = 200, abs_max = 5500) {
  #* split n_pixels based on theta into background and gaussians
  n_unif_pixels <- ceiling(n * theta)
  n_gauss_pixels <- floor(n * (1 - theta))
  #* background is uniform
  background <- runif(n_unif_pixels, 0, u2_max)
  #* simulate a number of gaussians randomly between 1 and u1_max
  n_gaussians <- runif(1, 1, u1_max)
  #* each gaussian has a mean that is uniform between 1 and u2_max
  mu_is <- lapply(seq_len(n_gaussians), function(i) {
    runif(1, 1, u2_max)
  })
  #* each gaussian has a sigma that is half-normal based on sd
  sd_is <- lapply(seq_len(n_gaussians), function(i) {
    extraDistr::rhnorm(1, sd)
  })
  #* assign pixels randomly to gaussians
  index <- sample(seq_len(n_gaussians), size = n_gauss_pixels, replace = TRUE)
  px_is <- lapply(seq_len(n_gaussians), function(i) {
    sum(index == i)
  })
  #* draws n_pixels time from each gaussian
  d <- unlist(lapply(seq_len(n_gaussians), function(i) {
    rnorm(px_is[[i]], mu_is[[i]], sd_is[[i]])
  }))
  #* combine data
  d <- c(d, background)
  #* make sure no gaussians return data out of bounds
  d[d < 0] <- runif(sum(d < 0), 0, abs_max)
  d[d >= abs_max] <- runif(sum(d >= abs_max), 0, abs_max)
  return(d)
}
lastNonZeroBin <- function(d) {
  max(d[d$value > 0, "label"])
}
tubeAngleToDepth <- function(x, theta) {
  sin(theta) * x
}
mv_mean <- function(d) {
  weighted.mean(d$label, d$value)
}
mv_median <- function(d) {
  median(rep(d$label, d$value))
}
mv_std <- function(d) {
  sd(rep(d$label, d$value))
}
sv_from_mv <- function(df, theta) { # note this should also return mean/median/std
  metaCols <- colnames(df)[-which(grepl("value|label|trait", colnames(df)))]
  out <- aggregate(
    as.formula(paste0(
      "value ~ ",
      paste(metaCols, collapse = "+")
    )),
    data = df, sum, na.rm = TRUE
  )
  colnames(out)[ncol(out)] <- "area"
  out$max_pixel <- unlist(lapply(split(df, interaction(df[, metaCols])), lastNonZeroBin))
  out$height <- tubeAngleToDepth(out$max_pixel, 0.35)
  out$mean_x_frequencies <- unlist(lapply(split(df, interaction(df[, metaCols])), mv_mean))
  out$median_x_frequencies <- unlist(lapply(split(df, interaction(df[, metaCols])), mv_median))
  out$std_x_frequencies <- unlist(lapply(split(df, interaction(df[, metaCols])), mv_std))
  return(out)
}

## -----------------------------------------------------------------------------
set.seed(123)
ex <- do.call(rbind, lapply(1:20, function(rep) {
  n_total_pixels <- runif(1, 100, 3000)
  x <- rRhyzoDist(n = n_total_pixels)
  h <- hist(x, plot = FALSE, breaks = seq(0, 5500, 20))
  breaks <- h$breaks[-1]
  counts <- h$counts
  data.frame(
    rep = as.character(rep),
    value = counts, label = breaks,
    trait = "x_frequencies"
  )
}))
pcv.joyplot(ex, "x_frequencies", group = c("rep"))

## -----------------------------------------------------------------------------
n_times <- 5
parameters <- data.frame(
  time = c(1:n_times),
  n_min = rep(0, n_times),
  n_max = seq(2000, 5500, length.out = n_times),
  theta = rep(0.3, n_times),
  u1_max = seq(10, 20, length.out = n_times),
  u2_max = seq(2000, 5500, length.out = n_times),
  u2_max_noise = seq(400, 200, length.out = n_times),
  sd = seq(150, 250, length.out = n_times)
)
set.seed(123)
df <- do.call(rbind, lapply(seq_len(nrow(parameters)), function(time) {
  pars <- parameters[parameters$time == time, ]
  do.call(rbind, lapply(1:10, function(rep) {
    n_total_pixels <- runif(1, pars$n_min, pars$n_max)
    u2_max_iter <- ceiling(rnorm(1, pars$u2_max, pars$u2_max_noise))
    x <- rRhyzoDist(
      n = n_total_pixels, theta = pars$theta,
      u1_max = pars$u1_max,
      u2_max = u2_max_iter, sd = pars$sd
    )
    h <- hist(x, plot = FALSE, breaks = seq(0, 5500, 20))
    breaks <- h$breaks[-1]
    counts <- h$counts
    data.frame(
      rep = as.character(rep),
      time = as.character(time),
      value = counts, label = breaks,
      trait = "x_frequencies"
    )
  }))
}))
df$rep <- factor(df$rep, levels = seq_along(unique(df$rep)), ordered = TRUE)
sv <- sv_from_mv(df)

## -----------------------------------------------------------------------------
n_times <- 5
parameters <- data.frame(
  time = c(1:n_times),
  n_min_new_px = rep(50, n_times),
  n_max_new_px = seq(1000, 2500, length.out = n_times),
  theta = rep(0.3, n_times),
  u1_max = seq(10, 20, length.out = n_times),
  mean_added_depth = seq(log(200), log(1200), length.out = n_times),
  added_depth_noise = rep(0.1, n_times),
  sd = seq(150, 250, length.out = n_times)
)
set.seed(567)
df2 <- do.call(rbind, lapply(1:10, function(rep) {
  dList <- list()
  add_area <- 2000
  previous_max_depth <- 2000
  d <- numeric(0)
  for (time in seq_len(nrow(parameters))) {
    pars <- parameters[parameters$time == time, ]
    n_total_pixels <- add_area + runif(1, pars$n_min_new_px, pars$n_max_new_px)
    max_depth <- ceiling(previous_max_depth + rlnorm(1, pars$mean_added_depth, pars$added_depth_noise))
    x <- rRhyzoDist(
      n = n_total_pixels, theta = pars$theta,
      u1_max = pars$u1_max,
      u2_max = max_depth, sd = pars$sd
    )
    d <- append(d, x)
    h <- hist(d, plot = FALSE, breaks = seq(0, 5500, 20))
    breaks <- h$breaks[-1]
    counts <- h$counts
    dList[[time]] <- data.frame(
      rep = as.character(rep),
      time = as.character(time),
      value = counts, label = breaks,
      trait = "x_frequencies"
    )
    add_area <- 0
    previous_max_depth <- max_depth
  }
  do.call(rbind, dList)
}))
df2$rep <- factor(df2$rep, levels = seq_along(unique(df2$rep)), ordered = TRUE)
sv2 <- sv_from_mv(df2)

## -----------------------------------------------------------------------------
ggplot(sv, aes(x = time, y = area, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Area (px)", title = "Assuming roots can leave the image")

ggplot(sv2, aes(x = time, y = area, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Area (px)", title = "Assuming roots can only enter the image")

## -----------------------------------------------------------------------------
ggplot(sv, aes(x = time, y = mean_x_frequencies, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Mean Depth", title = "Assuming roots can leave the image")

ggplot(sv2, aes(x = time, y = mean_x_frequencies, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Mean Depth", title = "Assuming roots can only enter the image")

## -----------------------------------------------------------------------------
ggplot(sv, aes(x = time, y = median_x_frequencies, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Median Depth", title = "Assuming roots can leave the image")

ggplot(sv2, aes(x = time, y = median_x_frequencies, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Median Depth", title = "Assuming roots can only enter the image")

## -----------------------------------------------------------------------------
ggplot(sv, aes(x = time, y = height, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Max Depth", title = "Assuming roots can leave the image")

ggplot(sv2, aes(x = time, y = height, group = rep)) +
  geom_point() +
  geom_line() +
  labs(x = "Sampling Time", y = "Max Depth", title = "Assuming roots can only enter the image")

## -----------------------------------------------------------------------------
sv$geno <- "a"
sv2$geno <- "b"
ex <- rbind(sv, sv2)
ex$time <- as.numeric(ex$time)

## -----------------------------------------------------------------------------
ss <- growthSS("int_linear", area ~ time | rep / geno, df = ex, type = "nls")
m1 <- fitGrowth(ss)

## -----------------------------------------------------------------------------
growthPlot(m1, form = ss$pcvrForm, df = ss$df)

## -----------------------------------------------------------------------------
testGrowth(ss, m1, test = "I")$anova

## -----------------------------------------------------------------------------
testGrowth(ss, m1, test = "A")$anova

## -----------------------------------------------------------------------------
table(ss$df$geno, ss$df$geno_numericLabel)
testGrowth(ss, m1, test = "A1*1.1 - A2")

## -----------------------------------------------------------------------------
s1 <- ex[ex$geno == "a" & ex$time == max(ex$time), "area"]
s2 <- ex[ex$geno == "b" & ex$time == max(ex$time), "area"]

conj_ex <- conjugate(
  s1, s2, # specify data, here two samples
  method = "t", # use the "T" distribution
  priors = list(mu = 3000, n = 1, s2 = 1000), # prior distribution, here it is the same for both samples
  plot = TRUE, # return a plot
  rope_range = c(-500, 500), # differences of <500 pixels deemed not meaningful
  rope_ci = 0.89, cred.int.level = 0.89, # default credible interval lengths
  hypothesis = "equal", # hypothesis to test
  support = NULL # optionally you can provide support, this is generally unnecessary.
)

## -----------------------------------------------------------------------------
conj_ex$summary

## -----------------------------------------------------------------------------
do.call(rbind, conj_ex$posterior)

## -----------------------------------------------------------------------------
conj_ex$plot

## -----------------------------------------------------------------------------
pcv.joyplot(df, "x_frequencies", group = c("rep", "time"))

## -----------------------------------------------------------------------------
getPeaks <- function(d = NULL, intensity = 20, duration = 3) {
  binwidth <- as.numeric(unique(diff(d$label)))
  if (length(binwidth) > 1) {
    stop("label column should have constant bin size")
  }
  labels <- sort(d[d$value >= intensity, "label"])
  r <- rle(diff(labels))
  peaks <- sum(r$lengths[r$values == binwidth] >= duration)
  return(peaks)
}

## -----------------------------------------------------------------------------
d <- split(df, interaction(df[, c("rep", "time")]))
peak_df <- data.frame(peaks = unlist(lapply(d, getPeaks)))
rownames(peak_df) <- NULL
peak_df$rep <- unlist(lapply(names(d), function(nm) {
  strsplit(nm, "[.]")[[1]][[1]]
}))
peak_df$time <- unlist(lapply(names(d), function(nm) {
  strsplit(nm, "[.]")[[1]][[2]]
}))

## -----------------------------------------------------------------------------
s1 <- peak_df[peak_df$time == min(peak_df$time), "peaks"]
s2 <- peak_df[peak_df$time == max(peak_df$time), "peaks"]

conj_ex2 <- conjugate(
  s1, s2, # specify data, here two samples
  method = "poisson", # use the Poisson distribution
  priors = list(a = c(0.5, 0.5), b = c(0.5, 0.5)), # prior distributions for gamma on lambda
  plot = TRUE, # return a plot
  rope_range = c(-1, 1), # differences of <500 pixels deemed not meaningful
  rope_ci = 0.89, cred.int.level = 0.89, # default credible interval lengths
  hypothesis = "equal", # hypothesis to test
  support = NULL # optionally you can provide support, this is generally unnecessary.
)

## -----------------------------------------------------------------------------
conj_ex2$summary

## -----------------------------------------------------------------------------
do.call(rbind, conj_ex2$posterior)

## -----------------------------------------------------------------------------
conj_ex2$plot

## -----------------------------------------------------------------------------
set.seed(123)

simFreqs <- function(vec, group) {
  s1 <- hist(vec, breaks = seq(1, 181, 1), plot = FALSE)$counts
  s1d <- as.data.frame(cbind(data.frame(group), matrix(s1, nrow = 1)))
  colnames(s1d) <- c("group", paste0("sim_", 1:180))
  s1d
}

sim_df <- rbind(
  do.call(rbind, lapply(1:10, function(i) {
    simFreqs(rnorm(200, 50, 10), group = "normal")
  })),
  do.call(rbind, lapply(1:10, function(i) {
    simFreqs(rlnorm(200, log(30), 0.25), group = "lognormal")
  })),
  do.call(rbind, lapply(1:10, function(i) {
    simFreqs(c(rlnorm(125, log(15), 0.25), rnorm(75, 75, 5)), group = "bimodal")
  })),
  do.call(rbind, lapply(1:10, function(i) {
    simFreqs(c(rlnorm(100, log(15), 0.25), rnorm(50, 50, 5), rnorm(50, 90, 5)), group = "trimodal")
  })),
  do.call(rbind, lapply(1:10, function(i) {
    simFreqs(runif(200, 1, 180), group = "uniform")
  }))
)

sim_df_long <- as.data.frame(data.table::melt(data.table::as.data.table(sim_df), id.vars = "group"))
sim_df_long$bin <- as.numeric(sub("sim_", "", sim_df_long$variable))

ggplot(sim_df_long, aes(x = bin, y = value, fill = group), alpha = 0.25) +
  geom_col(position = "identity", show.legend = FALSE) +
  pcv_theme() +
  facet_wrap(~group)

## ----class.source="simulated", class.output="simulated"-----------------------
sim_emd <- pcv.emd(
  df = sim_df, cols = "sim_", reorder = c("group"),
  mat = FALSE, plot = TRUE, parallel = 1, raiseError = TRUE
)
sim_emd$plot

## ----class.source="simulated", class.output="simulated"-----------------------
n <- pcv.net(sim_emd$data, filter = "0.5")
net.plot(n, fill = "group")

## -----------------------------------------------------------------------------
pcv.joyplot(df, "x_frequencies", group = c("rep", "time"))

## -----------------------------------------------------------------------------
df1_emd <- pcv.emd(
  df = df, cols = "x_frequencies", reorder = c("rep", "time"),
  id = c("rep", "time"),
  mat = FALSE, plot = TRUE, parallel = 1, raiseError = FALSE
)

## ----class.source="simulated", class.output="simulated"-----------------------
n <- pcv.net(df1_emd$data, filter = "0.75")
net.plot(n, fill = "time")

## -----------------------------------------------------------------------------
pcv.joyplot(df2, "x_frequencies", group = c("rep", "time"))

df2_emd <- pcv.emd(
  df = df2, cols = "x_frequencies", reorder = c("rep", "time"),
  id = c("rep", "time"),
  mat = FALSE, plot = TRUE, parallel = 1, raiseError = FALSE
)
n2 <- pcv.net(df2_emd$data, filter = "0.75")
net.plot(n2, fill = "time")


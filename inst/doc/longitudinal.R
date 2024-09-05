## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message = FALSE----------------------------------------------------------
library(pcvr)
library(data.table) # for fread
library(ggplot2)
library(patchwork) # for easy ggplot manipulation/combination
library(brms)

## ----eval = FALSE-------------------------------------------------------------
#  if (!"cmdstanr" %in% installed.packages()) {
#    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#  }
#  if (!"brms" %in% installed.packages()) {
#    install.packages("brms")
#  }
#  library(brms)
#  library(cmdstanr)
#  cmdstanr::install_cmdstan()

## -----------------------------------------------------------------------------
simdf <- growthSim("logistic", n = 20, t = 25, params = list(
  "A" = c(200, 160),
  "B" = c(13, 11),
  "C" = c(3, 3.5)
))
l <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Logistic") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("gompertz", n = 20, t = 25, params = list(
  "A" = c(200, 160),
  "B" = c(13, 11),
  "C" = c(0.2, 0.25)
))
g <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Gompertz") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("monomolecular", n = 20, t = 25, params = list("A" = c(200, 160), "B" = c(0.08, 0.1)))
m <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Monomolecular") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("exponential", n = 20, t = 25, params = list("A" = c(15, 20), "B" = c(0.095, 0.095)))
e <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Exponential") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("linear", n = 20, t = 25, params = list("A" = c(1.1, 0.95)))
ln <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Linear") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("power law", n = 20, t = 25, params = list("A" = c(16, 11), "B" = c(0.75, 0.7)))
pl <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Power Law") +
  theme_minimal() +
  theme(legend.position = "none")

patch <- (l + g + m) / (e + ln + pl)
patch

## -----------------------------------------------------------------------------
set.seed(345)
gomp <- growthSim("gompertz", n = 20, t = 35, params = list(
  "A" = c(200, 180, 160),
  "B" = c(20, 22, 18),
  "C" = c(0.15, 0.2, 0.1)
))


sigma_df <- aggregate(y ~ group + time, data = gomp, FUN = sd)

ggplot(sigma_df, aes(x = time, y = y, group = group)) +
  geom_line(color = "gray60") +
  pcv_theme() +
  labs(y = "SD of y", title = "Gompertz Sigma")

## ----message = FALSE----------------------------------------------------------
draw_gomp_sigma <- function(x) {
  23 * exp(-21 * exp(-0.22 * x))
}
ggplot(sigma_df, aes(x = time, y = y)) +
  geom_line(aes(group = group), color = "gray60") +
  geom_hline(aes(yintercept = 12, color = "Homoskedastic"),
    linetype = 5,
    key_glyph = draw_key_path
  ) +
  geom_abline(aes(slope = 0.8, intercept = 0, color = "Linear"),
    linetype = 5,
    key_glyph = draw_key_path
  ) +
  geom_smooth(
    method = "gam", aes(color = "Spline"), linetype = 5, se = FALSE,
    key_glyph = draw_key_path
  ) +
  geom_function(fun = draw_gomp_sigma, aes(color = "Gompertz"), linetype = 5) +
  scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, linetype = 1))) +
  pcv_theme() +
  theme(legend.position = "bottom") +
  labs(y = "SD of y", title = "Gompertz Sigma", color = "")

## -----------------------------------------------------------------------------
ss <- growthSS(
  model = "gompertz", form = y ~ time | id / group, sigma = "int",
  df = gomp, start = list("A" = 130, "B" = 15, "C" = 0.25)
)

## ----eval=FALSE---------------------------------------------------------------
#  fit_h <- fitGrowth(ss, iter = 1000, cores = 4, chains = 4, silent = 0)
#  
#  brmPlot(fit_h, form = ss$pcvrForm, df = ss$df)

## -----------------------------------------------------------------------------
ss <- growthSS(
  model = "gompertz", form = y ~ time | id / group, sigma = "linear",
  df = gomp, start = list("A" = 130, "B" = 15, "C" = 0.25)
)

## ----eval=FALSE---------------------------------------------------------------
#  fit_l <- fitGrowth(ss,
#    iter = 1000, cores = 4, chains = 4, silent = 0,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )
#  
#  p1 <- brmPlot(fit_l, form = ss$pcvrForm, df = ss$df)
#  p2 <- p1 + coord_cartesian(ylim = c(0, 300))
#  p <- p1 / p2
#  p

## -----------------------------------------------------------------------------
ss <- growthSS(
  model = "gompertz", form = y ~ time | id / group, sigma = "spline",
  df = gomp, start = list("A" = 130, "B" = 15, "C" = 0.25)
)

## ----eval=FALSE---------------------------------------------------------------
#  fit_s <- fitGrowth(ss,
#    iter = 2000, cores = 4, chains = 4, silent = 0,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )
#  
#  brmPlot(fit_s, form = ss$pcvrForm, df = ss$df)

## -----------------------------------------------------------------------------
ss <- growthSS(
  model = "gompertz", form = y ~ time | id / group, sigma = "gompertz",
  df = gomp, start = list(
    "A" = 130, "B" = 15, "C" = 0.25,
    "sigmaA" = 15, "sigmaB" = 15, "sigmaC" = 0.25
  ),
  type = "brms"
)

## ----eval=FALSE---------------------------------------------------------------
#  fit_g <- fitGrowth(ss,
#    iter = 2000, cores = 4, chains = 4, silent = 0,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )
#  
#  brmPlot(fit_g, form = ss$pcvrForm, df = ss$df)

## ----message = FALSE----------------------------------------------------------
draw_gomp_sigma <- function(x) {
  23 * exp(-21 * exp(-0.22 * x))
}
draw_logistic_sigma <- function(x) {
  20 / (1 + exp((15 - x) / 2))
}
draw_logistic_exp <- function(x) {
  2.5 * exp(0.08 * x)
}
draw_logistic_quad <- function(x) {
  (0.3 * x) + (0.02 * x^2)
}

ggplot(sigma_df, aes(x = time, y = y)) +
  geom_line(aes(group = group), color = "gray60", linetype = 5) +
  geom_hline(aes(yintercept = 12, color = "Homoskedastic"), linetype = 1) +
  geom_abline(aes(slope = 0.8, intercept = 0, color = "Linear"),
    linetype = 1,
    key_glyph = draw_key_path
  ) +
  geom_smooth(
    method = "gam", aes(color = "Spline"), linetype = 1, se = FALSE,
    key_glyph = draw_key_path
  ) +
  geom_function(fun = draw_gomp_sigma, aes(color = "Gompertz"), linetype = 1) +
  geom_function(fun = draw_logistic_sigma, aes(color = "Logistic"), linetype = 1) +
  geom_function(fun = draw_logistic_exp, aes(color = "Exponential"), linetype = 1) +
  geom_function(fun = draw_logistic_quad, aes(color = "Quadratic"), linetype = 1) +
  scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, linetype = 1))) +
  pcv_theme() +
  theme(legend.position = "bottom") +
  labs(y = "SD of y", title = "Gompertz Sigma", color = "")

## ----eval = FALSE-------------------------------------------------------------
#  loo_spline <- add_criterion(fit_s, "loo")
#  loo_homo <- add_criterion(fit_h, "loo")
#  loo_linear <- add_criterion(fit_l, "loo")
#  loo_gomp <- add_criterion(fit_g, "loo")
#  
#  h <- loo_homo$criteria$loo$estimates[3, 1]
#  s <- loo_spline$criteria$loo$estimates[3, 1]
#  l <- loo_linear$criteria$loo$estimates[3, 1]
#  g <- loo_gomp$criteria$loo$estimates[3, 1]
#  
#  loodf <- data.frame(loo = c(h, s, l, g), model = c("Homosked", "Spline", "Linear", "Gompertz"))
#  loodf$model <- factor(loodf$model, levels = unique(loodf$model[order(loodf$loo)]), ordered = TRUE)
#  
#  ggplot(
#    loodf,
#    aes(x = model, y = loo, fill = model)
#  ) +
#    geom_col() +
#    scale_fill_viridis_d() +
#    labs(y = "LOO Information Criteria", x = "Sub Model of Sigma") +
#    theme_minimal() +
#    theme(legend.position = "none")

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(345)
#  ln <- growthSim("linear", n = 5, t = 10, params = list("A" = c(2, 3, 10)))
#  
#  strongPrior <- prior(student_t(3, 0, 5), dpar = "sigma", class = "b") +
#    prior(gamma(2, 0.1), class = "nu", lb = 0.001) +
#    prior(normal(10, .05), nlpar = "A", lb = 0)
#  
#  ss <- growthSS(
#    model = "linear", form = y ~ time | id / group, sigma = "homo",
#    df = ln, priors = strongPrior
#  )
#  
#  fit <- fitGrowth(ss, iter = 1000, cores = 2, chains = 2, silent = 0)
#  
#  brmPlot(fit, form = ss$pcvrForm, df = ss$df) +
#    coord_cartesian(ylim = c(0, 100))

## ----eval = FALSE-------------------------------------------------------------
#  weakPrior <- prior(student_t(3, 0, 5), dpar = "sigma", class = "b") +
#    prior(gamma(2, 0.1), class = "nu", lb = 0.001) +
#    prior(lognormal(log(10), 0.25), nlpar = "A", lb = 0)
#  
#  ss <- growthSS(
#    model = "linear", form = y ~ time | id / group, sigma = "homo",
#    df = ln, priors = weakPrior
#  )
#  
#  fit <- fitGrowth(ss, iter = 1000, cores = 2, chains = 2, silent = 0)
#  
#  brmPlot(fit, form = ss$pcvrForm, df = ss$df) +
#    coord_cartesian(ylim = c(0, 100))

## -----------------------------------------------------------------------------
priors <- list("A" = 130, "B" = 10, "C" = 0.2)
priorPlots <- plotPrior(priors)
priorPlots[[1]] / priorPlots[[2]] / priorPlots[[3]]

## -----------------------------------------------------------------------------
twoPriors <- list("A" = c(100, 130), "B" = c(6, 12), "C" = c(0.5, 0.25))
plotPrior(twoPriors, "gompertz", n = 100)[[1]]

## -----------------------------------------------------------------------------
set.seed(123)
simdf <- growthSim("logistic", n = 20, t = 25, params = list(
  "A" = c(200, 160),
  "B" = c(13, 11),
  "C" = c(3, 3.5)
))

## -----------------------------------------------------------------------------
nls_ss <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  df = simdf, type = "nls"
)

## -----------------------------------------------------------------------------
nlrq_ss <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  df = simdf, type = "nlrq",
  tau = seq(0.01, 0.99, 0.04)
)

## -----------------------------------------------------------------------------
nlme_ss <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  df = simdf, sigma = "power", type = "nlme"
)

## -----------------------------------------------------------------------------
mgcv_ss <- growthSS(
  model = "gam", form = y ~ time | id / group,
  df = simdf, type = "mgcv"
)

## ----eval=FALSE---------------------------------------------------------------
#  brms_ss <- growthSS(
#    model = "logistic", form = y ~ time | id / group,
#    sigma = "spline", df = simdf,
#    start = list("A" = 130, "B" = 10, "C" = 1)
#  )

## -----------------------------------------------------------------------------
ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group))

## ----warning=FALSE, message=FALSE---------------------------------------------
nls_fit <- fitGrowth(nls_ss)
nlrq_fit <- fitGrowth(nlrq_ss)
nlme_fit <- fitGrowth(nlme_ss)
mgcv_fit <- fitGrowth(mgcv_ss)

## ----eval=FALSE---------------------------------------------------------------
#  brms_fit <- fitGrowth(brms_ss,
#    iter = 500, cores = 1, chains = 1,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )

## ----warning=FALSE, message=FALSE---------------------------------------------
growthPlot(nls_fit, form = nls_ss$pcvrForm, df = nls_ss$df)
growthPlot(nlrq_fit, form = nlrq_ss$pcvrForm, df = nlrq_ss$df)
growthPlot(nlme_fit, form = nlme_ss$pcvrForm, df = nlme_ss$df)
growthPlot(mgcv_fit, form = mgcv_ss$pcvrForm, df = mgcv_ss$df)

## ----eval=FALSE---------------------------------------------------------------
#  growthPlot(brms_fit, form = brms_ss$pcvrForm, df = brms_ss$df)

## -----------------------------------------------------------------------------
testGrowth(nls_ss, nls_fit, test = "A")$anova

## -----------------------------------------------------------------------------
testGrowth(nlrq_ss, nlrq_fit, test = "A")[["0.49"]]

## ----warning=FALSE, message=FALSE---------------------------------------------
testGrowth(nlme_ss, nlme_fit, test = "A")$anova

## -----------------------------------------------------------------------------
testGrowth(mgcv_ss, mgcv_fit)$anova

## -----------------------------------------------------------------------------
testGrowth(fit = nls_fit, test = list(
  "A1 - A2 *1.1",
  "(B1+1) - B2",
  "C1 - (C2-0.5)",
  "A1/B1 - (1.1 * A2/B2)"
))

## -----------------------------------------------------------------------------
testGrowth(fit = nlme_fit, test = list(
  "(A.groupa / A.groupb) - 0.9",
  "1 + (B.groupa - B.groupb)",
  "C.groupa/C.groupb - 1"
))

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  (hyp <- brms::hypothesis(brms_fit, "(A_groupa) > 1.1 * (A_groupb)"))

## ----echo=FALSE, eval=TRUE----------------------------------------------------
structure(list(
  Hypothesis = "((A_groupa))-(1.1*(A_groupb)) > 0",
  Estimate = 16.7880224, Est.Error = 3.13518501598807, CI.Lower = 12.012705,
  CI.Upper = 21.952505, Evid.Ratio = Inf, Post.Prob = 1, Star = "*"
), row.names = c(NA, -1L), class = "data.frame")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  simdf <- growthSim(
#    model = "linear + linear",
#    n = 20, t = 25,
#    params = list("linear1A" = c(15, 12), "changePoint1" = c(8, 6), "linear2A" = c(3, 5))
#  )
#  
#  ss <- growthSS(
#    model = "linear + linear", form = y ~ time | id / group, sigma = "spline",
#    start = list("linear1A" = 10, "changePoint1" = 5, "linear2A" = 2),
#    df = simdf, type = "brms"
#  )
#  
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  simdf <- growthSim("linear + logistic",
#    n = 20, t = 25,
#    params = list(
#      "linear1A" = c(15, 12), "changePoint1" = c(8, 6),
#      "logistic2A" = c(100, 150), "logistic2B" = c(10, 8),
#      "logistic2C" = c(3, 2.5)
#    )
#  )
#  
#  ss <- growthSS(
#    model = "linear + logistic", form = y ~ time | id / group, sigma = "spline",
#    list(
#      "linear1A" = 10, "changePoint1" = 5,
#      "logistic2A" = 100, "logistic2B" = 10, "logistic2C" = 3
#    ),
#    df = simdf, type = "brms"
#  )
#  
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  ss <- growthSS(
#    model = "linear + gam", form = y ~ time | id / group, sigma = "int",
#    list("linear1A" = 10, "changePoint1" = 5),
#    df = simdf, type = "brms"
#  )
#  
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  simdf <- growthSim("linear + linear + linear",
#    n = 25, t = 50,
#    params = list(
#      "linear1A" = c(10, 12), "changePoint1" = c(8, 6),
#      "linear2A" = c(1, 2), "changePoint2" = c(25, 30), "linear3A" = c(20, 24)
#    )
#  )
#  
#  ss <- growthSS(
#    model = "linear + linear + linear", form = y ~ time | id / group, sigma = "spline",
#    list(
#      "linear1A" = 10, "changePoint1" = 5,
#      "linear2A" = 2, "changePoint2" = 15,
#      "linear3A" = 5
#    ), df = simdf, type = "brms"
#  )
#  
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  set.seed(123)
#  noise <- do.call(rbind, lapply(1:30, function(i) {
#    chngpt <- rnorm(2, 18, 2)
#    rbind(
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[1], group = "a",
#        y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
#      ),
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[2], group = "b",
#        y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
#      )
#    )
#  }))
#  noise2 <- do.call(rbind, lapply(1:30, function(i) {
#    start1 <- max(noise[noise$id == paste0("id_", i) & noise$group == "a", "time"])
#    start2 <- max(noise[noise$id == paste0("id_", i) & noise$group == "b", "time"])
#    rbind(
#      data.frame(
#        id = paste0("id_", i), time = start1:40, group = "a",
#        y = c(runif(length(start1:40), 15, 50))
#      ),
#      data.frame(
#        id = paste0("id_", i), time = start2:40, group = "b",
#        y = c(runif(length(start2:40), 15, 50))
#      )
#    )
#  }))
#  simdf <- rbind(noise, noise2)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  ss <- growthSS(
#    model = "int + int", form = y ~ time | id / group, sigma = "int + int",
#    list(
#      "int1" = 10, "changePoint1" = 10, "int2" = 20, # main model
#      "sigmaint1" = 10, "sigmachangePoint1" = 10, "sigmaint2" = 10
#    ), # sub model
#    df = simdf, type = "brms"
#  )
#  
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  set.seed(123)
#  noise <- do.call(rbind, lapply(1:30, function(i) {
#    chngpt <- rnorm(2, 18, 2)
#    rbind(
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[1], group = "a",
#        y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
#      ),
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[2], group = "b",
#        y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
#      )
#    )
#  }))
#  signal <- growthSim("linear",
#    n = 30, t = 20,
#    params = list("A" = c(3, 5))
#  )
#  signal <- do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int) {
#    noisesub <- noise[paste0(noise$id, noise$group) == int, ]
#    signalSub <- signal[paste0(signal$id, signal$group) == int, ]
#    y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
#    signalSub$time <- signalSub$time + max(noisesub$time)
#    signalSub$y <- y_end + signalSub$y
#    signalSub
#  }))
#  simdf <- rbind(noise, signal)

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  ss <- growthSS(
#    model = "int + linear", form = y ~ time | id / group, sigma = "int + linear",
#    list(
#      "int1" = 10, "changePoint1" = 10, "linear2A" = 20,
#      "sigmaint1" = 10, "sigmachangePoint1" = 10, "sigmalinear2A" = 10
#    ),
#    df = simdf, type = "brms"
#  )
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  set.seed(123)
#  noise <- do.call(rbind, lapply(1:30, function(i) {
#    chngpt <- rnorm(2, 18, 2)
#    rbind(
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[1], group = "a",
#        y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
#      ),
#      data.frame(
#        id = paste0("id_", i), time = 1:chngpt[2], group = "b",
#        y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
#      )
#    )
#  }))
#  signal <- growthSim("logistic",
#    n = 20, t = 30,
#    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#  )
#  signal <- do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int) {
#    noisesub <- noise[paste0(noise$id, noise$group) == int, ]
#    signalSub <- signal[paste0(signal$id, signal$group) == int, ]
#    y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
#    signalSub$time <- signalSub$time + max(noisesub$time)
#    signalSub$y <- y_end + signalSub$y
#    signalSub
#  }))
#  simdf <- rbind(noise, signal)
#  simdf <- simdf[simdf$time < 45, ]

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  ss <- growthSS(
#    model = "int+logistic", form = y ~ time | id / group, sigma = "int + spline",
#    list(
#      "int1" = 5, "changePoint1" = 10,
#      "logistic2A" = 130, "logistic2B" = 10, "logistic2C" = 3,
#      "sigmaint1" = 5, "sigmachangePoint1" = 15
#    ),
#    df = simdf, type = "brms"
#  )
#  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)

## ----eval=FALSE---------------------------------------------------------------
#  print(load(url("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata")))
#  from3to25 <- list(
#    fit_3, fit_5, fit_7, fit_9, fit_11, fit_13,
#    fit_15, fit_17, fit_19, fit_21, fit_23, fit_25
#  )
#  
#  distributionPlot(fits = from3to25, form = y ~ time | id / group, params = c("A", "B", "C"), d = simdf)


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libraries, message=FALSE-------------------------------------------------
library(pcvr)
library(brms) # for rvon_mises
library(ggplot2)
library(patchwork) # for easy ggplot manipulation/combination

## ----mv gaussian sim----------------------------------------------------------
mv_gauss <- mvSim(
  dists = list(
    rnorm = list(mean = 50, sd = 10),
    rnorm = list(mean = 60, sd = 12)
  ),
  n_samples = c(30, 40)
)

## ----run conjugate------------------------------------------------------------
vm_ex1 <- conjugate(
  s1 = mv_gauss[1:30, -1],
  s2 = mv_gauss[31:70, -1],
  method = "vonmises",
  priors = list(mu = 45, kappa = 1, boundary = c(0, 180), known_kappa = 1, n = 1),
  plot = TRUE, rope_range = c(-5, 5), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

## ----show summary-------------------------------------------------------------
vm_ex1$summary

## ----heavy conj plot, eval = FALSE--------------------------------------------
#  vm_ex1$plot # not printed due to being a very dense ggplot

## ----run conjugate 2----------------------------------------------------------
vm_ex1_1 <- conjugate(
  s1 = rnorm(30, 50, 10),
  s2 = rnorm(40, 60, 12),
  method = "vonmises",
  priors = list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1),
  plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)
do.call(rbind, vm_ex1_1$posterior)

## ----run conjugate 3----------------------------------------------------------
set.seed(42)
vm_ex2 <- conjugate(
  s1 = brms::rvon_mises(100, -3.1, 2),
  s2 = brms::rvon_mises(100, 3.1, 2),
  method = "vonmises",
  priors = list(mu = 0, kappa = 1, known_kappa = 2),
  plot = TRUE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

## ----show conjugate 3 out data------------------------------------------------
vm_ex2$summary
do.call(rbind, vm_ex2$posterior)

## ----conjugate 3 plot, eval = FALSE-------------------------------------------
#  vm_ex2$plot # not printed due to being a very dense ggplot

## ----polar coord conj1, eval = FALSE------------------------------------------
#  p <- vm_ex2$plot
#  p[[1]] <- p[[1]] +
#    ggplot2::coord_polar() +
#    ggplot2::scale_y_continuous(limits = c(-pi, pi))

## ----polar coord conj2, eval=FALSE, echo=FALSE--------------------------------
#  coord_polar_linear <- function(theta = "x", start = 0, direction = 1, clip = "on") {
#    theta <- match.arg(theta, c("x", "y"))
#    r <- if (theta == "x") {
#      "y"
#    } else {
#      "x"
#    }
#    ggproto(NULL, CoordPolar, theta = theta, r = r, start = start,
#            direction = sign(direction), clip = clip,
#            is_linear = function() {TRUE})
#  }
#  
#  p <- vm_ex2$plot
#  p[[1]] <- p[[1]] +
#    coord_polar_linear() +
#    ggplot2::scale_y_continuous(limits = c(-pi, pi))
#  p

## ----run conjugate 4----------------------------------------------------------
vm2_ex1 <- conjugate(
  s1 = mv_gauss[1:30, -1],
  s2 = mv_gauss[31:70, -1],
  method = "vonmises2",
  priors = list(mu = 45, kappa = 1, boundary = c(0, 180), n = 1),
  plot = TRUE, rope_range = c(-5, 5), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

## ----show conj 4 posterior----------------------------------------------------
do.call(rbind, vm2_ex1$posterior)

## ----run conjugate 5----------------------------------------------------------
set.seed(42)
vm2_ex2 <- conjugate(
  s1 = brms::rvon_mises(100, -3.1, 2),
  s2 = brms::rvon_mises(100, 3.1, 2),
  method = "vonmises2",
  priors = list(mu = 0, kappa = 1),
  plot = TRUE, rope_range = c(-0.75, 0.75), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

## ----plot from conjugate 5----------------------------------------------------
vm2_ex2$plot # much lighter to print this since it is in radians

## ----posterior from conjugate 5-----------------------------------------------
do.call(rbind, vm2_ex2$posterior)

## ----pcvr model setup---------------------------------------------------------
nReps <- 25
time <- 1:20
muTrend1 <- -2 + (0.25 * time)
muTrend2 <- -1 + (0.2 * time)
kappaTrend1 <- (0.5 * time)
kappaTrend2 <- (0.3 * time)
set.seed(123)
vm2 <- do.call(rbind, lapply(1:nReps, function(rep) {
  do.call(rbind, lapply(time, function(ti) {
    v1 <- brms::rvon_mises(1, muTrend1[ti], kappaTrend1[ti])
    v2 <- brms::rvon_mises(1, muTrend2[ti], kappaTrend2[ti])
    data.frame(y = c(v1, v2), x = ti, group = c("a", "b"), rep = rep)
  }))
}))

ss <- growthSS(
  model = "von_mises: int_linear", form = y ~ x | rep / group, sigma = "int", df = vm2,
  start = NULL, type = "brms"
)
ss$prior # default priors
ss$formula # formula specifies kappa based on sigma argument

## ----brms example 1, eval = FALSE---------------------------------------------
#  set.seed(123)
#  n <- 1000
#  vm1 <- data.frame(
#    x = c(brms::rvon_mises(n, 1.5, 3), brms::rvon_mises(n, 3, 2)),
#    y = rep(c("a", "b"), each = n)
#  )
#  
#  basePlot <- ggplot(vm1, aes(x = x, fill = y)) +
#    geom_histogram(binwidth = 0.1, alpha = 0.75, position = "identity") +
#    labs(fill = "Group") +
#    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
#    scale_fill_viridis_d() +
#    theme_minimal() +
#    theme(legend.position = "bottom")
#  
#  basePlot +
#    coord_polar() +
#    scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3.1415), labels = c(-2, -1, 0, 1, 2, "Pi"))
#  basePlot + scale_x_continuous(breaks = c(-round(pi, 2), -1.5, 0, 1.5, round(pi, 2)))
#  
#  prior1 <- set_prior("student_t(3,0,2.5)", coef = "ya") +
#    set_prior("student_t(3,0,2.5)", coef = "yb") +
#    set_prior("normal(5.0, 0.8)", coef = "ya", dpar = "kappa") +
#    set_prior("normal(5.0, 0.8)", coef = "yb", dpar = "kappa")
#  
#  fit1 <- brm(bf(x ~ 0 + y, kappa ~ 0 + y),
#    family = von_mises,
#    prior = prior1,
#    data = vm1,
#    iter = 1000, cores = 2, chains = 2, backend = "cmdstanr", silent = 0, init = 0,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )
#  fit1
#  
#  x <- brmsfamily("von_mises")
#  pars <- colMeans(as.data.frame(fit))
#  mus <- pars[grepl("b_y", names(pars))]
#  x$linkinv(mus) # inverse half tangent function
#  # should be around 1.5, 3
#  kappas <- pars[grepl("kappa", names(pars))]
#  exp(kappas) # kappa is log linked
#  # should be around 3, 2
#  pred_draws <- as.data.frame(predict(fit1, newdata = data.frame(y = c("a", "b")), summary = FALSE))
#  preds <- data.frame(
#    draw = c(pred_draws[, 1], pred_draws[, 2]),
#    y = rep(c("a", "b"), each = nrow(pred_draws))
#  )
#  predPlot <- ggplot(preds, aes(x = draw, fill = y)) +
#    geom_histogram(binwidth = 0.1, alpha = 0.75, position = "identity") +
#    labs(fill = "Group", y = "Predicted Draws") +
#    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
#    scale_fill_viridis_d() +
#    theme_minimal() +
#    theme(legend.position = "bottom")
#  predPlot + scale_x_continuous(breaks = c(-round(pi, 2), -1.5, 0, 1.5, round(pi, 2)))
#  predPlot +
#    coord_polar() +
#    scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3.1415), labels = c(-2, -1, 0, 1, 2, "Pi"))

## ----brms example 2, eval = FALSE---------------------------------------------
#  nReps <- 25
#  time <- 1:20
#  muTrend1 <- -2 + (0.25 * time)
#  muTrend2 <- -1 + (0.2 * time)
#  kappaTrend1 <- (0.5 * time)
#  kappaTrend2 <- (0.3 * time)
#  set.seed(123)
#  vm2 <- do.call(rbind, lapply(1:nReps, function(rep) {
#    do.call(rbind, lapply(time, function(ti) {
#      v1 <- rvon_mises(1, muTrend1[ti], kappaTrend1[ti])
#      v2 <- rvon_mises(1, muTrend2[ti], kappaTrend2[ti])
#      data.frame(y = c(v1, v2), x = ti, group = c("a", "b"), rep = rep)
#    }))
#  }))
#  
#  ggplot(vm2, aes(x = x, y = y, color = group, group = interaction(group, rep))) +
#    geom_line() +
#    labs(y = "Y (Von Mises)") +
#    theme_minimal()
#  
#  ggplot(vm2, aes(y = x, x = y, color = group, group = interaction(group, rep), alpha = x)) +
#    geom_line() +
#    labs(y = "Time", x = "Von Mises") +
#    theme_minimal() +
#    guides(alpha = "none") +
#    coord_polar() +
#    scale_x_continuous(
#      breaks = c(-2, -1, 0, 1, 2, 3.1415),
#      limits = c(-pi, pi),
#      labels = c(-2, -1, 0, 1, 2, "Pi")
#    )
#  
#  prior2 <- set_prior("normal(5,0.8)", nlpar = "K") +
#    set_prior("student_t(3, 0, 2.5)", nlpar = "I") +
#    set_prior("student_t(3, 0, 2.5)", nlpar = "M")
#  
#  fit2 <- brm(
#    bf(y ~ I + M * x,
#      nlf(kappa ~ K * x),
#      I + M ~ 0 + group,
#      K ~ 0 + group,
#      autocor = ~ arma(x | rep:group, 1, 1),
#      nl = TRUE
#    ),
#    family = von_mises,
#    prior = prior2,
#    data = vm2,
#    iter = 2000, cores = 4, chains = 4, backend = "cmdstanr", silent = 0, init = 0,
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  )
#  fit2
#  pars <- colMeans(as.data.frame(fit2))
#  pars[grepl("^b_", names(pars))]
#  
#  outline <- data.frame(
#    group = rep(c("a", "b"), each = 20),
#    x = rep(1:20, 2)
#  )
#  probs <- seq(0.01, 0.99, 0.02)
#  preds <- cbind(outline, predict(fit2, newdata = outline, probs = probs))
#  
#  pal <- viridis::plasma(n = length(probs))
#  p2 <- ggplot(preds, aes(y = x)) +
#    facet_wrap(~group) +
#    lapply(seq(1, 49, 2), function(lower) {
#      geom_ribbon(aes(xmin = .data[[paste0("Q", lower)]], xmax = .data[[paste0("Q", 100 - lower)]]),
#        fill = pal[lower]
#      )
#    }) +
#    theme_minimal() +
#    coord_polar() +
#    scale_x_continuous(
#      breaks = c(-2, -1, 0, 1, 2, 3.1415),
#      limits = c(-pi, pi),
#      labels = c(-2, -1, 0, 1, 2, "Pi")
#    )


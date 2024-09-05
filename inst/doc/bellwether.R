## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(pcvr)
library(data.table) # for fread
library(ggplot2)
library(patchwork) # for easy ggplot manipulation/combination

## ----class.source="static-code", eval=F---------------------------------------
#  devtools::install_github("joshqsumner/pcvr", build_vignettes = TRUE)
#  library(pcvr)

## ----class.source="static-code", eval=F---------------------------------------
#  complicatedFunction("syntax") # do not run this style

## -----------------------------------------------------------------------------
1 + 1 # run this style

## ----class.source="simulated", eval=T-----------------------------------------
support <- seq(0, 1, 0.0001) # this style is simulated data
plot(support, dbeta(support, 5, 5), type = "l", main = "simulated example")

## -----------------------------------------------------------------------------
base_url <- "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/"
base_url2 <- "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/"
sv <- read.pcv(
  filepath = paste0(base_url, "pcv4-single-value-traits.csv"),
  reader = "fread"
)

## -----------------------------------------------------------------------------
key <- read.csv(paste0(base_url, "smallPhenotyperRun_key.csv"))
head(key)
sv <- merge(sv, key, by = "barcode")
table(sv$genotype, sv$fertilizer)

## -----------------------------------------------------------------------------
genotype <- substr(sv$barcode, 3, 5)
genotype <- ifelse(genotype == "002", "B73",
  ifelse(genotype == "003", "W605S",
    ifelse(genotype == "004", "MM", "Mo17")
  )
)
fertilizer <- substr(sv$barcode, 8, 8)
fertilizer <- ifelse(fertilizer == "A", 100,
  ifelse(fertilizer == "B", 50, 0)
)
table(genotype, fertilizer)

## -----------------------------------------------------------------------------
chip_size_px <- mean(c(sv$median_color_chip_height_median, sv$median_color_chip_width_median)) # ~52

px_per_cm <- chip_size_px / 1.2 # ~ 43.5
pixels_per_cmsq <- px_per_cm^2 # ~ 1890

sv$area_cm2 <- sv$area_pixels / pixels_per_cmsq
sv$height_cm <- sv$height_pixels / px_per_cm

## ----class.source="static-code", eval=F---------------------------------------
#  example <- read.pcv("prohibitivelyLargeFile.csv",
#    filters = list(
#      "trait in area_pixels, area_above_reference_pixels, area_below_reference_pixels",
#      "sample in default"
#    )
#  )

## ----class.source="static-code", eval=F---------------------------------------
#  sv <- read.pcv(paste0(base_url2, "smallPhenotyperRun.csv"),
#    mode = "wide",
#    reader = "fread"
#  )
#  
#  if (TRUE) { # we can parse barcodes for the metadata that we need
#    sv$genotype <- substr(sv$barcode, 3, 5)
#    sv$genotype <- ifelse(sv$genotype == "002", "B73",
#      ifelse(sv$genotype == "003", "W605S",
#        ifelse(sv$genotype == "004", "MM", "Mo17")
#      )
#    )
#    sv$fertilizer <- substr(sv$barcode, 8, 8)
#    sv$fertilizer <- ifelse(sv$fertilizer == "A", "100",
#      ifelse(sv$fertilizer == "B", "50", "0")
#    )
#  } else { # or we might use a key file and join it to our data
#    key <- read.csv(paste0(base_url, "smallPhenotyperRun_key.csv"))
#    sv <- merge(sv, key, by = "barcode")
#  }

## ----eval=F-------------------------------------------------------------------
#  onlyPhenos <- read.pcv.3(file = paste0(base_url, "pcv3Phenos.csv"), metaCol = NULL)
#  colnames(onlyPhenos)
#  
#  phenosAndMeta <- read.pcv.3(
#    file = paste0(base_url, "pcv3Phenos.csv"), metaCol = "meta",
#    metaForm = "vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",
#    joinSnapshot = "id"
#  )
#  colnames(phenosAndMeta)
#  
#  all <- read.pcv.3(
#    file = paste0(base_url, "pcv3Phenos.csv"),
#    snapshotFile = paste0(base_url, "pcv3Snapshot.csv"),
#    designFile = paste0(base_url, "pcv3Design.csv"),
#    metaCol = "meta", metaForm = "vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",
#    joinSnapshot = "id", conversions = list(area = 13.2 * 3.7 / 46856)
#  )
#  colnames(all)

## -----------------------------------------------------------------------------
out <- bw.time(sv,
  plantingDelay = 0, phenotype = "area_pixels", cutoff = 10,
  timeCol = "timestamp", group = c("barcode", "rotation"), plot = TRUE
)
out$plot
sv <- out$data
dim(sv)

## -----------------------------------------------------------------------------
checkGroups(sv, c("DAS", "barcode", "rotation", "genotype", "fertilizer"))

## -----------------------------------------------------------------------------
phenotypes <- colnames(sv)[c(19:35, 43:45, 48:49)]
phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
groupForm <- "DAS+timestamp+barcode+genotype+fertilizer"
form <- as.formula(paste0(phenoForm, "~", groupForm))
sv_ag_with_outliers <- aggregate(form, data = sv, mean, na.rm = TRUE)
dim(sv_ag_with_outliers)

## -----------------------------------------------------------------------------
out <- bw.outliers(
  df = sv_ag_with_outliers, phenotype = "area_pixels",
  group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode")
)
sv_ag <- out$data
out$plot
dim(sv_ag)

## -----------------------------------------------------------------------------
checkGroups(sv_ag, c("DAS", "barcode", "genotype", "fertilizer"))

## -----------------------------------------------------------------------------
water <- bw.water(paste0(base_url, "metadata.json"))

water$genotype <- substr(water$barcode, 3, 5)
water$genotype <- ifelse(water$genotype == "002", "B73",
  ifelse(water$genotype == "003", "W605S",
    ifelse(water$genotype == "004", "MM", "Mo17")
  )
)
water$fertilizer <- substr(water$barcode, 8, 8)
water$fertilizer <- ifelse(water$fertilizer == "A", "100",
  ifelse(water$fertilizer == "B", "50", "0")
)

ggplot(water[water$weight_after != -1, ], aes(
  x = DAS,
  y = water_amount, group = barcode, color = genotype
)) +
  facet_wrap(~ factor(fertilizer, levels = c("0", "50", "100"))) +
  geom_line() +
  pcv_theme() +
  guides(color = guide_legend(title = "Condition", override.aes = list(linewidth = 5))) +
  labs(y = "Watering Amount (g)") +
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
test <- pwue(df = sv_ag, w = water, pheno = "area_pixels", time = "timestamp", id = "barcode")

ggplot(test, aes(x = DAS, y = pWUE, color = genotype, group = barcode)) +
  geom_line() +
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  labs(y = expression("Pseudo WUE (" ~ frac(
    Delta ~ textstyle("Area")[" pixels"],
    Delta ~ textstyle("Weight")[" g"]
  ) ~ ")")) +
  pcv_theme()

## ----warning = FALSE, message = FALSE-----------------------------------------
frem(sv_ag,
  des = c("genotype", "fertilizer"),
  phenotypes = c(
    "area_pixels", "area_above_reference_pixels", "area_below_reference_pixels",
    "convex_hull_area_pixels", "convex_hull_vertices_none", "ellipse_angle_degrees",
    "ellipse_eccentricity_none", "ellipse_major_axis_pixels", "ellipse_minor_axis_pixels",
    "height_pixels", "height_above_reference_pixels", "height_below_reference_pixels",
    "horizontal_reference_position_none", "hue_circular_mean_degrees", "hue_circular_std_degrees",
    "hue_median_degrees", "perimeter_pixels", "solidity_none", "width_pixels"
  ),
  timeCol = "DAS", cor = TRUE, returnData = FALSE, combine = FALSE, markSingular = TRUE, time = NULL
)

## ----warning = FALSE, message = FALSE-----------------------------------------
frem(sv_ag,
  des = c("genotype", "fertilizer"),
  phenotypes = c(
    "area_pixels", "area_above_reference_pixels", "area_below_reference_pixels",
    "convex_hull_area_pixels", "convex_hull_vertices_none", "ellipse_angle_degrees",
    "ellipse_eccentricity_none", "ellipse_major_axis_pixels",
    "ellipse_minor_axis_pixels",
    "height_pixels", "height_above_reference_pixels", "height_below_reference_pixels",
    "horizontal_reference_position_none", "hue_circular_mean_degrees", "hue_circular_std_degrees",
    "hue_median_degrees", "perimeter_pixels", "solidity_none", "width_pixels"
  ),
  timeCol = "DAS", cor = FALSE, returnData = FALSE, combine = FALSE, markSingular = FALSE, time = "all"
)

## -----------------------------------------------------------------------------
ggplot(sv_ag, aes(
  x = DAS, y = area_cm2, group = interaction(genotype, fertilizer, lex.order = TRUE),
  color = genotype
)) +
  facet_wrap(~ factor(fertilizer, levels = c("0", "50", "100"))) +
  geom_smooth(method = "loess", se = TRUE, fill = "gray90") +
  geom_line(aes(group = barcode), linewidth = 0.15) +
  labs(
    y = expression("Area" ~ "(cm"^2 ~ ")"),
    color = "Genotype"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  pcv_theme() +
  theme(axis.text.x.bottom = element_text(angle = 0))

## -----------------------------------------------------------------------------
mo17_area <- sv_ag[sv_ag$genotype == "Mo17" & sv_ag$DAS > 18 & sv_ag$fertilizer == 100, "area_cm2"]
b73_area <- sv_ag[sv_ag$genotype == "B73" & sv_ag$DAS > 18 & sv_ag$fertilizer == 100, "area_cm2"]

area_res_t <- conjugate(s1 = mo17_area, s2 = b73_area, method = "t", plot = TRUE, rope_range = c(-5, 5))

## -----------------------------------------------------------------------------
rt <- relativeTolerance(sv_ag,
  phenotypes = c("area_cm2", "height_cm"),
  grouping = c("fertilizer", "genotype", "DAS"), control = "fertilizer", controlGroup = "100"
)


ggplot(
  rt[rt$phenotype == "area_cm2" & rt$DAS %in% c(10:12), ],
  aes(x = DAS, y = mu_rel, fill = interaction(fertilizer, genotype))
) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mu_rel - 1.96 * se_rel, ymax = mu_rel + 1.96 * se_rel),
    position = position_dodge(width = 0.9), width = 0.3
  ) +
  pcv_theme() +
  labs(y = "Relative Tolerance", fill = "Fertilizer\nand Genotype")

## -----------------------------------------------------------------------------
pd <- rt[rt$phenotype == "area_cm2" & rt$DAS %in% c(5:19) & rt$fertilizer == "0", ]
pd$upper_2se <- pd$mu_rel + 2 * pd$se_rel
pd$lower_2se <- pd$mu_rel - 2 * pd$se_rel

ggplot(pd, aes(x = DAS, y = mu_rel, fill = genotype)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(
    ymin = lower_2se, ymax = upper_2se,
    group = genotype
  ), position = "dodge") +
  pcv_theme() +
  labs(y = "Relative Tolerance")

## -----------------------------------------------------------------------------
cp <- cumulativePheno(sv_ag, phenotypes = c("area_cm2", "height_cm"), group = "barcode", timeCol = "DAS")

## -----------------------------------------------------------------------------
ggplot(cp, aes(x = DAS, y = area_cm2_csum, color = genotype, group = barcode)) +
  facet_wrap(~ factor(fertilizer, levels = c("0", "50", "100"))) +
  geom_line() +
  pcv_theme() +
  labs(
    y = expression("Cumulative Sum of Area" ~ "(cm"^2 ~ ")"),
    color = "Genotype"
  )

## ----class.source="simulated", class.output="simulated"-----------------------
simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
l <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Logistic") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("gompertz",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.2, 0.25))
)
g <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Gompertz") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("monomolecular",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(0.08, 0.1))
)
m <- ggplot(simdf, aes(time, y, group = interaction(group, id))) +
  geom_line(aes(color = group)) +
  labs(title = "Monomolecular") +
  theme_minimal() +
  theme(legend.position = "none")

simdf <- growthSim("exponential",
  n = 20, t = 25,
  params = list("A" = c(15, 20), "B" = c(0.095, 0.095))
)
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

(l + g + m) / (e + ln + pl)

## -----------------------------------------------------------------------------
priors <- list("A" = 130, "B" = 10, "C" = 0.2)
priorPlots <- plotPrior(priors)
priorPlots[[1]] / priorPlots[[2]] / priorPlots[[3]]

## -----------------------------------------------------------------------------
twoPriors <- list("A" = c(100, 130), "B" = c(6, 12), "C" = c(0.5, 0.25))
plotPrior(twoPriors, "gompertz", n = 100)[[1]]

## -----------------------------------------------------------------------------
sv_ag$group <- interaction(sv_ag$fertilizer, sv_ag$genotype)

## ----eval=FALSE---------------------------------------------------------------
#  library(brms)
#  library(cmdstanr)
#  cmdstanr::install_cmdstan()

## ----eval=FALSE---------------------------------------------------------------
#  ss <- growthSS(
#    model = "gompertz", form = area_cm2 ~ DAS | barcode / group, sigma = "spline", df = sv_ag,
#    start = list("A" = 130, "B" = 10, "C" = 0.5), type = "brms"
#  )

## -----------------------------------------------------------------------------
ggplot(sv_ag, aes(x = DAS, y = area_cm2, group = barcode, color = group)) +
  geom_line() +
  theme_minimal() +
  labs(
    y = expression("Area" ~ "(cm"^2 ~ ")"),
    color = "Genotype\nand Soil"
  )

## ----class.source="static-code", eval=FALSE-----------------------------------
#  fit <- fitGrowth(ss,
#    iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",
#    control = list(adapt_delta = 0.999, max_treedepth = 20)
#  ) # options to increase performance

## ----eval=FALSE---------------------------------------------------------------
#  brmPlot(fit, form = area_cm2 ~ DAS | barcode / group, df = ss$df) +
#    labs(y = expression("Area" ~ "(cm"^2 ~ ")"))

## ----eval=FALSE---------------------------------------------------------------
#  brmViolin(
#    model = fit, params = NULL,
#    hyp = "num/denom>1.05", compareX = c("0.B73", "50.B73", "100.B73"), againstY = "0.B73",
#    group_sep = "[.]", groups_into = c("soil", "genotype"), x = "soil", facet = "genotype",
#    returnData = FALSE
#  )

## -----------------------------------------------------------------------------
ggplot(sv_ag[sv_ag$DAS == 18, ], aes(
  x = fertilizer, y = hue_circular_mean_degrees,
  fill = as.character(fertilizer)
)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 0.5) +
  scale_fill_manual(values = c(viridis::viridis(3, 1, 0.1)), breaks = c("0", "50", "100")) +
  pcv_theme() +
  theme(legend.position = "none") +
  facet_wrap(~genotype, scales = "free_x") +
  scale_x_discrete(limits = c("0", "50", "100")) +
  labs(y = "Hue Circular Mean (degrees)", x = "Soil and Genotype")

## -----------------------------------------------------------------------------
hue_wide <- read.pcv(paste0(base_url2, "pcv4-multi-value-traits.csv"),
  mode = "wide", reader = "fread"
)
hue_wide$genotype <- substr(hue_wide$barcode, 3, 5)
hue_wide$genotype <- ifelse(hue_wide$genotype == "002", "B73",
  ifelse(hue_wide$genotype == "003", "W605S",
    ifelse(hue_wide$genotype == "004", "MM", "Mo17")
  )
)
hue_wide$fertilizer <- substr(hue_wide$barcode, 8, 8)
hue_wide$fertilizer <- ifelse(hue_wide$fertilizer == "A", "100",
  ifelse(hue_wide$fertilizer == "B", "50", "0")
)
hue_wide <- bw.time(hue_wide, timeCol = "timestamp", group = "barcode", plot = FALSE)

phenotypes <- colnames(hue_wide)[grepl("hue_frequencies", colnames(hue_wide))]
phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
groupForm <- "DAS+barcode+genotype+fertilizer"
form <- as.formula(paste0(phenoForm, "~", groupForm))
hue_wide <- aggregate(form, data = hue_wide, mean, na.rm = TRUE)

## ----message=FALSE------------------------------------------------------------
p <- pcv.joyplot(hue_wide[hue_wide$DAS == 18, ],
  index = "hue_frequencies",
  group = c("fertilizer", "genotype")
)
p + scale_fill_gradientn(colors = scales::hue_pal(l = 65)(360)) +
  scale_y_discrete(limits = c("0", "50", "100"))

## ----message=FALSE------------------------------------------------------------
p <- pcv.joyplot(hue_wide[hue_wide$DAS %in% c(5, 10, 15), ],
  index = "hue_frequencies", group = c("fertilizer", "genotype"),
  y = "DAS", id = NULL
)
p + scale_fill_gradientn(colors = scales::hue_pal(l = 65)(360)) +
  scale_y_discrete(limits = c("5", "10", "15"))

## -----------------------------------------------------------------------------
mo17_sample <- hue_wide[
  hue_wide$genotype == "Mo17" & hue_wide$DAS > 18 & hue_wide$fertilizer == 100,
  grepl("hue_freq", colnames(hue_wide))
]
b73_sample <- hue_wide[
  hue_wide$genotype == "B73" & hue_wide$DAS > 18 & hue_wide$fertilizer == 100,
  grepl("hue_freq", colnames(hue_wide))
]

hue_res_ln <- conjugate(
  s1 = mo17_sample, s2 = b73_sample, method = "lognormal",
  plot = TRUE, rope_range = c(-10, 10), hypothesis = "equal"
)

## -----------------------------------------------------------------------------
pcadf(hue_wide, cols = "hue_frequencies", color = "genotype", returnData = FALSE) +
  facet_wrap(~ factor(fertilizer, levels = c("0", "50", "100")))

## ----class.source="simulated", class.output="simulated"-----------------------
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
n <- pcv.net(sim_emd$data, filter = 0.5)
net.plot(n, fill = "group")

## ----class.source="simulated", class.output="simulated"-----------------------
n <- pcv.net(sim_emd$data, filter = "0.5")
net.plot(n, fill = "group")

## ----eval=FALSE---------------------------------------------------------------
#  EMD <- pcv.emd(
#    df = hue_wide[hue_wide$DAS %in% c(5, 12, 19), ], cols = "hue_frequencies",
#    reorder = c("fertilizer", "genotype", "DAS"),
#    mat = FALSE, plot = TRUE, parallel = 12, raiseError = TRUE
#  )

## ----eval=FALSE, include=FALSE------------------------------------------------
#  head(EMD$data)
#  EMD$plot

## -----------------------------------------------------------------------------
hue_ag1 <- mv_ag(df = hue_wide, group = c("DAS", "genotype", "fertilizer"), n_per_group = 2)
dim(hue_ag1)

hue_ag2 <- mv_ag(hue_wide, group = c("DAS", "genotype", "fertilizer"), n_per_group = 1)
dim(hue_ag2)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(456)
#  net <- pcv.net(EMD$data, meta = c("fertilizer", "genotype", "DAS"), filter = 0.5)

## ----eval=FALSE---------------------------------------------------------------
#  net.plot(net, fill = "DAS", shape = "fertilizer", size = 2)


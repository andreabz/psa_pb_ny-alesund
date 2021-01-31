#### Potential source areas for atmospheric lead reaching Ny-Ålesund from 
#### 2010 to 2018
#### Andrea Bazzano, Stefano Bertinetti, Francisco Ardini, David Cappelletti and
#### Marco Grotti
#### corresponding author: andrea.bazzano@edu.unige.it
#### Supplementary material: script to reproduce the data analysis described in
#### the manuscript for chemical and isotopic data.
############################### Preamble ######################################

#### Loading of specialised libraries ----
# data manipulation
library(data.table)
library(dplyr)
library(lubridate)
# statistics
library(summarytools)
library(fitdistrplus)
library(dunn.test)
library(mclust)
library(QuantPsyc)
library(energy)
library(MASS)
# graphics
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggforce)
library(patchwork)
library(scales)

# color-blind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#CC79A7", "#D55E00")

#### Additional functions ----
# opposite of %in%
`%!in%` = Negate(`%in%`)


# Deming regression: noise_ratio is the ratio of errors on y and x
# https://en.wikipedia.org/wiki/Deming_regression
# https://stackoverflow.com/questions/26995923/
#   ggplot2-how-to-plot-an-orthogonal-regression-line
deming.fit <- function(x, y, noise_ratio = sd(y)/sd(x)) {
  if(missing(noise_ratio) || 
     is.null(noise_ratio)) noise_ratio <- 
      eval(formals(sys.function(0))$noise_ratio) 
  # this is just a complicated way to write `sd(y)/sd(x)`
  delta <-  noise_ratio^2
  x_name <- deparse(substitute(x))
  
  s_yy <- var(y)
  s_xx <- var(x)
  s_xy <- cov(x, y)
  beta1 <- (s_yy - delta*s_xx + 
              sqrt((s_yy - delta*s_xx)^2 + 4*delta*s_xy^2)) / (2*s_xy)
  beta0 <- mean(y) - beta1 * mean(x) 
  
  res <- c(beta0 = beta0, beta1 = beta1)
  names(res) <- c("(Intercept)", x_name)
  class(res) <- "Deming"
  res
}

deming <- function(formula, data, R = 100, noise_ratio = NULL, ...){
  ret <- boot::boot(
    data = model.frame(formula, data), 
    statistic = function(data, ind) {
      data <- data[ind, ]
      args <- rlang::parse_exprs(colnames(data))
      names(args) <- c("y", "x")
      rlang::eval_tidy(rlang::expr(deming.fit(!!!args, 
                                              noise_ratio = noise_ratio)), 
                       data, env = rlang::current_env())
    },
    R=R
  )
  class(ret) <- c("Deming", class(ret))
  ret  
}

predictdf.Deming <- function(model, xseq, se, level) {
  pred <- as.vector(tcrossprod(model$t0, cbind(1, xseq)))
  if(se) {
    preds <- tcrossprod(model$t, cbind(1, xseq))
    data.frame(
      x = xseq,
      y = pred,
      ymin = apply(preds, 2, function(x) quantile(x, probs = (1-level)/2)),
      ymax = apply(preds, 2, function(x) quantile(x, probs = 1-((1-level)/2)))
    )
  } else {
    return(data.frame(x = xseq, y = pred))
  }
}


# percentage of end-member a and b in Pb isotope ratio values in x
# x, a and b must have two columns named pb208206 and pb207206.
# b has higher isotope ratio values than a.
# x has isotope ratio values between a and b.
bperc <- function(x, a, b){
  ax1 <- x$pb208206 - a$pb208206
  ax2 <- x$pb207206 - a$pb207206
  ab1 <- b$pb208206 - a$pb208206
  ab2 <- b$pb207206 - a$pb207206
  
  100 * (1 - sqrt((ax1^2 + ax2^2) / (ab1^2 + ab2^2)))
}

# statistics for the median
median_stats <- function(x) {
                    median_x = median(x, na.rm = TRUE)
                    iqr_x = IQR(x, na.rm = TRUE)
                    mad_x = mad(x, na.rm = TRUE)
                    q1 = quantile(x, probs = 0.25, na.rm = TRUE)
                    q2 = quantile(x, probs = 0.75, na.rm = TRUE)
              
                    list(
                      median_x = median_x,
                      iqr_x = iqr_x,
                      mad_x = mad_x,
                      q1 = q1,
                      q2 = q2)}

############################### Datasets ######################################

#### PM10 samples collected at Ny-Ålesund from 2010 to 2018 ----
# date reports the sampling data in YYYY-MM-DD format.
# volume is the sampling volume in m3.
# pb_sign is = for Pb concentrarion data above limit of quantification (LoQ)
# and < for data below LoQ.
# pb_val is numeric and it is the measured Pb concentration or LoQ in pg/m3.
# pb is text and it is the measured Pb concentration or <LoQ in pg/m3.
# al_ef is the enrichment factor (EF) EF(Pb/Al)c in comparison to the upper
# continental crust (UCC, Wedepohl 1995).
# pb20x20y is the value measured for 20xPb / 20yPb isotope ratio.
# u20x20y is the 95-confidence level uncertainty for the measured 20xPb / 20yPb
# isotope ratio value.
# Missing values are reported as NA
nya_pm <- fread(input = "dataset/nyalesund_2010-2018.csv", skip = 15)[, 
            date := as.POSIXct(date)][,                     # convert to date
                `:=` (year_col = year(date),                # extract the year
                      month_col = month(date,               # extract the month
                                        label = TRUE,
                                        locale = "en_GB.UTF8"))][     
                      !is.na(pb208206) & !is.na(pb207206),  # number of sample
                      `:=` (n208206 = .N,                       
                            n207206 = .N),
                      by = .(year_col, month_col)][, 
            season := dplyr::case_when(                     # group by season
              month_col %in% c("Feb", "Mar", "Apr", "May") ~ "spring",
              month_col %in% c("Jun", "Jul", "Aug", "Sep", "Oct") ~ "summer")]

#### Pb isotope ratios for soils and sediments near Ny-Ålesund ----
# pb20x20y columns should be read as 20xPb / 20yPb.
# see reference column for author, year citations.
# Missing values are reported as NA.
mineral <- fread("dataset/soils_sediments.csv", skip = 5)

#### Pb isotope ratios for Atmospheric particulate in some relevant areas ----
# pb20x20y columns should be read as 20xPb / 20yPb.
# see reference column for author, year citations.
# Missing values are reported as NA.
atm_nh <- fread("dataset/atm_nh.csv", skip = 8)
atm_nh[, zone := ifelse(area == "russia" | area == "europe",
                        "eurasia",
                        area)]

#### Pb isotope ratios for Chinese coals ----
# pb20x20y columns should be read as 20xPb / 20yPb.
# data from Bi et al, (2017).
# Missing values are reported as NA.
chinese_coals <- fread("dataset/chinese_coals.csv", skip = 3)

#### Pb isotope ratios for possible end-members ----
# Some values resulted from the following elaborations.
# pb20x20y columns should be read as 20xPb / 20yPb.
# see reference column for author, year citations.
# nudge_x and nudge_y will be used in plots as graphical parameters.
# Missing values are reported as NA.
endmembers <- fread("dataset/endmembers.csv", skip = 8)

############################# Section 3.1 #####################################
# Only the data mentioned in this section of the manuscript are reported.
# Data for tables and figures are reported in a separate section at the end
# of this file.

#### Pb concentration (pg/m3) ----
# mean and median
nya_pm[, .(
  mean(pb_val, na.rm = TRUE),
  median(pb_val, na.rm = TRUE))]

# median of values above the third-quartile
nya_pm[pb_val > quantile(pb_val, 0.75, na.rm = TRUE), .(
  median(pb_val, na.rm = TRUE))]

# estimate density function to describe the data distribution
qqnorm(nya_pm$pb_val)       # not normal
qqnorm(log(nya_pm$pb_val))  # lognormal
pb_ft <- fitdist(nya_pm$pb_val, distr = "lnorm") # parameters in log scale
# goodness of fit test with Kolmogorov-Smirnov test
ks.test(nya_pm$pb_val, "plnorm",
        meanlog = pb_ft$estimate[1],
        sdlog = pb_ft$estimate[2]) # p-value = 0.51.

#### Enrichment factors (EF) EF(Pb/Al)c ----
# mean
nya_pm[, .(
  mean(al_ef, na.rm = TRUE))]

# estimate density function to describe the data distribution
qqnorm(nya_pm$al_ef)       # not normal
qqnorm(log(nya_pm$al_ef))  # lognormal
ef_ft <- fitdist(nya_pm[!is.na(al_ef), al_ef], distr = "lnorm") # in log scale
# goodness of fit test with Kolmogorov-Smirnov test
ks.test(nya_pm$al_ef, "plnorm",
        meanlog = ef_ft$estimate[1],
        sdlog = ef_ft$estimate[2]) # p-value = 0.95.

#### 208Pb/206Pb and 207Pb/206Pb isotope ratio values ----
# correlation
cor.test(nya_pm$pb208206, nya_pm$pb207206, use = "complete.obs")
# Mardia's test for multivariate normality
mult.norm(nya_pm[!is.na(pb208206), .(pb208206, pb207206)]) # p-value < 0.001
# E-test for multivariate normality
mvnorm.etest(nya_pm[!is.na(pb208206),
                            .(pb208206, pb207206)],
                     R = 1000)  # p-value < 0.0001

#### Seasonality ----
# Mann-Whitney test on medians
nya_pm[, lapply(.SD, function(x) wilcox.test(x ~ season)),
       .SDcols = c("pb_val", "al_ef", "pb208206", "pb207206")]
# p-value < 0.001 for all investigated variables

############################# Section 3.2 #####################################
# Only the data mentioned in this section of the manuscript are reported.
# Data for tables and figures are reported in a separate section at the end
# of this file.

#### 208Pb/206Pb and 207Pb/206Pb vs 1/EF correlation ----
# correlation
cor.test(nya_pm$pb208206, 1/nya_pm$al_ef) # p-value < 0.001
cor.test(nya_pm$pb207206, 1/nya_pm$al_ef) # p-value < 0.001

# divide data into EF classes
nya_pm[, al_ef_bin := dplyr::case_when(
  al_ef < 10 ~ "<10",
  al_ef >= 10 & al_ef < 25 ~ "[10, 25)",
  al_ef >= 25 & al_ef < 50 ~ "[25, 50)",
  al_ef >= 50 & al_ef < 100 ~ "[50, 100)",
  al_ef >= 100 & al_ef < 1000 ~ "≥100")][,
    al_ef_bin := factor(al_ef_bin,
      levels = c(
                 "<10",
                 "[10, 25)",
                 "[25, 50)",
                 "[50, 100)",
                 "≥100"),
      labels = c(
        expression(EF(Pb/Al)[c]*" < 10"),
        expression("10 ≤ "*EF(Pb/Al)[c]*" < 25"),
        expression("25 ≤ "*EF(Pb/Al)[c]*" < 50"),
        expression("50 ≤ "*EF(Pb/Al)[c]*" < 100"),
        expression(EF(Pb/Al)[c]*" ≥ 100")),
        ordered = TRUE)]

# linear models
nya_pm[, lapply(.SD, median, na.rm = TRUE),             # median by EF classes
       .SDcols = c("pb208206", "pb207206", "al_ef"), 
       by = al_ef_bin] %>% 
  lm(data =., cbind(pb208206, pb207206) ~ I(1/al_ef)) %>% 
  summary

nya_pm[, lapply(.SD, median, na.rm = TRUE),             # median by EF classes
       .SDcols = c("pb208206", "pb207206", "al_ef"), 
       by = al_ef_bin] %>% 
  rlm(data =., pb208206 ~ I(1/al_ef)) %>% 
  summary

nya_pm[, lapply(.SD, median, na.rm = TRUE),             # median by EF classes
       .SDcols = c("pb208206", "pb207206", "al_ef"), 
       by = al_ef_bin] %>% 
  rlm(data =., pb207206 ~ I(1/al_ef)) %>% 
  summary 
# same result by ordinary least-squares and iterated re-weighted least squares

#### Predict possible end-members for EF = 1 and EF ~ Inf ----
nya_pm[, lapply(.SD, median, na.rm = TRUE),             # median by EF classes
       .SDcols = c("pb208206", "pb207206", "al_ef"), 
       by = al_ef_bin] %>% 
  lm(data =., cbind(pb208206) ~ I(1/al_ef)) %>%
  predict(newdata = data.frame(al_ef = c(1, 10^6)),
          interval = "confidence")
# EF = 1: 208Pb/206Pb = 2.041 +- 0.011
# EF ~ Inf: 208Pb/206Pb = 2.100 +- 0.003

nya_pm[, lapply(.SD, median, na.rm = TRUE),             # median by EF classes
       .SDcols = c("pb208206", "pb207206", "al_ef"), 
       by = al_ef_bin] %>% 
  rlm(data =., cbind(pb207206) ~ I(1/al_ef)) %>%
  predict(newdata = data.frame(al_ef = c(1, 10^6)),
          interval = "confidence")
# EF = 1: 207Pb/206Pb = 0.8098 +- 0.012
# EF ~ Inf: 207Pb/206Pb = 0.863 +- 0.001
# values for calculated anthropogenic end-members (EF ~ Inf) and natural
# end-member (EF = 1) in "endmembers" data.frame.

#### Estimate mean crustal contribution ----
nya_pm[, lapply(.SD, mean, na.rm = TRUE),
       .SDcols = c("pb208206", "pb207206")] %>%
  bperc(b = endmembers[3, ],
        a = endmembers[4, ])

max_nat <- nya_pm[, lapply(.SD, mean, na.rm = TRUE),
            .SDcols = c("pb208206", "pb207206")] %>%
              bperc(b = endmembers[3, c(2, 3)]+ c(0.003, 0.001),
                    a = endmembers[4, c(2, 3)] + c(0.011 , 0.012))

min_nat <- nya_pm[, lapply(.SD, mean, na.rm = TRUE),
            .SDcols = c("pb208206", "pb207206")] %>%
              bperc(b = endmembers[3, c(2, 3)] - c(0.003, 0.001),
                    a = endmembers[4, c(2, 3)] - c(0.011 , 0.012))
# 5-16% with a mean contribution of ~10%

#### Estimate crustal contribution to Pb concentration (pg/m3) ----
nya_pm[, median(pb_val, na.rm = TRUE)] * c(max_nat, min_nat) / 100
# 1-4 pg/m3

############################# Section 3.3 #####################################
# Only the data mentioned in this section of the manuscript are reported.
# Data for tables and figures are reported in a separate section at the end
# of this file.

#### Percentage of significantly enriched samples ----
nya_pm[al_ef > 10, .N]/nya_pm[!is.na(al_ef), .N]*100
# 83% of data with reported EF(Pb/Al)c has EF > 10

#### Median +- IQR Pb isotope ratio values for the two seasons ----
nya_pm[, .(median = lapply(.SD, median_stats, na.rm = TRUE),
           iqr = lapply(.SD, IQR, na.rm = TRUE)),
       .SDcols = c("pb208206", "pb207206")]
# 208Pb/206Pb = 2.098+-0.016 and 207Pb/206Pb = 0.861+-0.009 (median +- IQR)

nya_pm[, .(median = lapply(.SD, median, na.rm = TRUE),
          iqr = lapply(.SD, IQR, na.rm = TRUE)),
.SDcols = c("pb208206", "pb207206"),
by = "season"]
# ratio       season median   iqr
# 208Pb/206Pb spring  2.101  0.01
# 207Pb/206Pb spring  0.862 0.004
# 208Pb/206Pb summer  2.091 0.021
# 207Pb/206Pb summer  0.857 0.011

#### Gaussian Mixture Modeling (GMM) ----
# modeling
nya_gmm <- Mclust(
            as.matrix(
              nya_pm[!is.na(pb208206) & 
                     !is.na(al_ef), 
                        .(pb207206,
                          pb208206,
                          log_ef = log10(al_ef),     # log EFs for normality
                          log_pb = log10(pb_val))]   # log Pb  for normality
              )
            )

# summary of the model
nya_gmm_summary <-summary(nya_gmm, parameters = TRUE)
# cluster 1 is cluster A: Central Asia
# cluster 2 is cluster B: North America

# bootstrap for confidence intervals of the parameters (slow)
nya_gmm_boot <- MclustBootstrap(nya_gmm, nboot = 10^3)
nya_gmm_boot_summary <- summary(nya_gmm_boot, "ci")
# cluster 1 is cluster A: Central Asia
# cluster 2 is cluster B: North America

# saving parameters in a data.frame
nya_gmm_parameters <- data.table(
                        parameters = c(rep("mean", 2),             # mean
                                       rep("sd", 2),               # sd
                                       rep("ci95_sup", 2), # 95conf. lvl. upper 
                                       rep("ci95_inf", 2), # 95conf. lvl. lower
                                       rep("df_mean", 2)), # density f. at mean
                        cluster = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2),
                        rbind(
                          t(nya_gmm_summary$mean),
                          sqrt(diag(nya_gmm_summary$variance[,,1])),
                          sqrt(diag(nya_gmm_summary$variance[,,2])),
                          summary(nya_gmm_boot, "ci")$mean[,,1][2,],
                          summary(nya_gmm_boot, "ci")$mean[,,2][2,],
                          summary(nya_gmm_boot, "ci")$mean[,,1][1,],
                          summary(nya_gmm_boot, "ci")$mean[,,2][1,],
                          dnorm(nya_gmm_summary$mean[1, ],
                                mean = nya_gmm_summary$mean[1, ],
                                sd = sqrt(diag(nya_gmm_summary$variance[,,1]))),
                          dnorm(nya_gmm_summary$mean[2, ],
                                mean = nya_gmm_summary$mean[2, ],
                                sd = sqrt(diag(nya_gmm_summary$variance[,,2])))
                          )
                      )

# add clusters to the PM10 samples dataset
nya_pm[!is.na(pb208206) & 
       !is.na(al_ef), 
       cluster := nya_gmm_summary$classification]

#### Characterizing the clusters -----

# proportion of samples from spring and summer in the two clusters
nya_pm[!is.na(cluster), .(season, cluster)] %>% 
  table %>% 
  prop.table(margin = 2)

# test for significance of the different number of samples
nya_pm[!is.na(cluster), .(season, cluster)] %>% 
  table %>% 
  chisq.test
# p-value < 0.001

# Pb concentrations and EFs for the two clusters
nya_pm[!is.na(cluster),
       .(median = lapply(.SD, median, na.rm = TRUE),
            iqr = lapply(.SD, mad, na.rm = TRUE)),
        .SDcols = c("pb_val", "al_ef"),
       by = "cluster"]

# testing for differences
wilcox.test(data = nya_pm, pb_val ~ cluster, na.action = na.exclude)
wilcox.test(data = nya_pm, al_ef ~ cluster, na.action = na.exclude)
# p-value < 0.001

# proportion of anthropogenic Pb associated to each cluster
nya_pm[!is.na(cluster) & al_ef > 10,
       .(pb = sum(pb_val, na.rm = TRUE)), 
       by = cluster][,
                     .(pb_prop = pb/sum(pb),
                       cluster)]

############################# Section 3.4 #####################################
# Only the data mentioned in this section of the manuscript are reported.
# Data for tables and figures are reported in a separate section at the end
# of this file.

#### Temporal variations in Pb concentration, EF and isotope ratios ----

#### Seasonality ----

wilcox.test(data = nya_pm, pb_val ~ season, na.action = na.exclude)
wilcox.test(data = nya_pm, al_ef ~ season, na.action = na.exclude)
wilcox.test(data = nya_pm, pb208206 ~ season, na.action = na.exclude)
wilcox.test(data = nya_pm, pb207206 ~ season, na.action = na.exclude)
# p-value < 0.001

# variations by month
kruskal.test(data = nya_pm, pb_val ~ month_col, na.action = na.exclude)
kruskal.test(data = nya_pm, al_ef ~ month_col, na.action = na.exclude)
kruskal.test(data = nya_pm, pb208206 ~ month_col, na.action = na.exclude)
kruskal.test(data = nya_pm, pb207206 ~ month_col, na.action = na.exclude)
# p-value < 0.001

# Post-Hoc Dunn's test for pair comparisons
with(nya_pm, dunn.test(pb_val, month_col, method = "bh"))
with(nya_pm, dunn.test(al_ef, month_col, method = "bh"))
with(nya_pm, dunn.test(pb208206, month_col, method = "bh"))
with(nya_pm, dunn.test(pb207206, month_col, method = "bh"))
# significant differences only between spring and summer months

#### Inter-annual variations ----

kruskal.test(data = nya_pm, pb_val ~ year_col, na.action = na.exclude)
p-value = 0.2
kruskal.test(data = nya_pm, al_ef ~ year_col, na.action = na.exclude)
kruskal.test(data = nya_pm, pb208206 ~ year_col, na.action = na.exclude)
kruskal.test(data = nya_pm, pb207206 ~ year_col, na.action = na.exclude)
# p-value < 0.001

# Post-Hoc Dunn's test for pair comparisons
with(nya_pm, dunn.test(al_ef, year_col, method = "bh"))

nya_pm[year_col == "2018",
       .(median = median(al_ef, na.rm = TRUE),
         iqr = mad(al_ef, na.rm = TRUE))]
# EF = 21 +- 18
nya_pm[year_col != "2018",
       .(median = median(al_ef, na.rm = TRUE),
         iqr = mad(al_ef, na.rm = TRUE))]
# EF = 38 +- 37
# Pb sampled in 2018 is less anthropogenically enriched compared to other years

with(nya_pm, dunn.test(pb208206, year_col, method = "bh"))
with(nya_pm, dunn.test(pb207206, year_col, method = "bh"))

nya_pm[, .(median = lapply(.SD, median, na.rm = TRUE),
           iqr = lapply(.SD, IQR, na.rm = TRUE)),
       .SDcols = c("pb208206", "pb207206"),
       by = "year_col"]
# 2015 and 2018 have lower isotope ratio values compared to other years

# taking into account only summer samples
with(nya_pm[season == "summer"], 
     dunn.test(pb208206, year_col, method = "bh"))
with(nya_pm[season == "summer"],
     dunn.test(pb207206, year_col, method = "bh"))

nya_pm[season == "summer",
       .(median = lapply(.SD, median, na.rm = TRUE),
           iqr = lapply(.SD, IQR, na.rm = TRUE),
         n = sum(!is.na(pb208206))),
       .SDcols = c("pb208206", "pb207206"),
       by = "year_col"]
# differences confirmed

#### Contributions of the two clusters ----

# version 1
# nya_pm[!is.na(cluster), 
#        .(pb = sum(pb_val, na.rm = TRUE), 
#                           pbm = mean(pb_val, na.rm = TRUE)),
#        by = .(cluster, year_col, season)][,
#        .(pb = pb / sum(pb) * pbm,
#          pb_prop = pb/sum(pb),
#           cluster),
#        by = .(year_col, season)][
#         order(year_col, season, cluster)][, 
#        .(median = median(pb_prop),
#          q1 = quantile(pb_prop, 0.25),
#          q2 = quantile(pb_prop, 0.75)),
#        by = c("cluster", "season")]

# Relative contribution: version 2
nya_pm[season == "spring",.(year_col, cluster)] %>% 
  table %>% 
  prop.table(margin = 1) %>%
  apply(2, median_q1q3)

nya_pm[season == "spring",. (year_col, cluster)] %>% 
  table %>% 
  pairwise.prop.test(p.adjust.method = "BH")

nya_pm[season == "summer",. (year_col, cluster)] %>% 
  table %>% 
  prop.table(margin = 1) %>%
  apply(2, median_q1q3)

nya_pm[season == "summer",. (year_col, cluster)] %>% 
  table %>% 
  pairwise.prop.test(p.adjust.method = "BH")
# 2015 and 2018 have higher proportion of cluster 2 compared with 2011, 2012

# Contribution to Pb concentration (pb/m3)
nya_pm[!is.na(cluster),
       lapply(.SD, median_stats), 
       .SDcols = c("pb_val"),
       by = c("season", "cluster")]

wilcox.test(data = nya_pm[cluster == 1], pb_val ~ season)
wilcox.test(data = nya_pm[cluster == 2], pb_val ~ season)

################################ Tables #######################################

#### Table 1 ----
descr(nya_pm[, .(pb_val, al_ef, pb208206, pb207206)], round.digits = 3)
# n > LoQ for Pb
nya_pm[pb_sign == "=", .N]

#### Table 2 ----
nya_gmm_parameters[1:4, ]
nya_gmm_summary$variance

#### Table S1 ----
nya_pm[, .(id,
           date,
           pb_conc = ifelse(pb_sign == "<", paste0(pb_sign, pb_val), pb_val),
           efs = al_ef,
           pb208206,
           u208206,
           pb207206,
           u207206)]

################################ Figures #######################################

### Setting the frame for a three isotope plot with data for atmoshperic Pb
theme.size = 20
geom.text.size = (theme.size - 4) * 0.352777778

# labels
isotext <- data.table(label = c("China", 
                                "USA", 
                                "Canada", 
                                "Europe and\n Russia"),
                      f207206 = c(0.862, 0.825, 0.86, 0.89),
                      f208206 = c(2.1216, 2.025, 2.075, 2.125),
                      f208204 = c(38.5, 39, 38.5, 37),
                      f207204 = c(15.65, 15.8, 15.75, 15.5))

# common base plot
iso208207206_atm <- ggplot() +
  stat_ellipse(data = atm_nh,
               aes(x = pb207206,
                   y = pb208206,
                   group = zone),
               geom = "polygon",
               type = "norm",
               fill = "white",
               col = "black") +
  geom_text_repel(data = isotext,
                  aes(x = f207206,
                      y = f208206,
                      label = label),
                  nudge_x = c(-0.02, +0.015, 0.03, 0.03),
                  nudge_y = c(0.025, -0.015, -0.025, -0.05),
                  size = geom.text.size) +
  scale_y_continuous(breaks = seq(from = 1.7,
                                  to = 2.3,
                                  by = 0.02),
                     limits = c(1.7, 2.3)) +
  scale_x_continuous(breaks = seq(from = 0.7,
                                  to = 1.0,
                                  by = 0.01),
                     limits = c(0.7, 1.0)) +
  labs(x = expression({}^207*"Pb/"*{}^206*"Pb"),
       y = expression({}^208*"Pb/"*{}^206*"Pb"),
       col = element_blank()) +
  theme_bw(base_size = theme.size)

#### Figure 1 ----
iso208207206_atm +
  # 95% confidence interval ellipses for PM10 data sorted by EF classes
  stat_ellipse(data = nya_pm[!is.na(al_ef_bin)],
               aes(
                 y = pb208206,
                 x = pb207206,
                 fill = al_ef_bin,
               ),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.7) +
  # 95% confidence interval ellipses for sediment Pb isotope ratio values
  stat_ellipse(data = mineral[site %like% "sediments"],
               aes(
                 x = pb207206,
                 y = pb208206),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.5) +
  # confidence ellipses around mean isotope ratio for EF classes
  stat_conf_ellipse(data = nya_pm[!is.na(al_ef_bin)],
                    aes(y = pb208206,
                        x = pb207206,
                        fill = al_ef_bin),
                    size = 0.5,
                    col = "black",
                    linetype = "solid",
                    geom = "polygon") +
  # Pb isotope ratio values for PM10 a Ny-Ålesund
  geom_point(data = nya_pm,
             aes(
               x = pb207206,
               y = pb208206,
               shape = "PM10"
             ),
             size = 2,
             col = "gray20") +
  # Pb isotope ratio values for possibile mineral contributions
  geom_point(data = mineral,
             aes(
               x = pb207206,
               y = pb208206,
               shape = site
             ),
             size = 2,
             col = "black") +
  # regression line through EF medians
  geom_smooth(data = nya_pm[, 
                      lapply(.SD, median, na.rm = TRUE),
                      .SDcols = c("pb_val", "pb208206",
                                  "pb207206", "al_ef"),
                      by = .(al_ef_bin)][!is.na(al_ef_bin)],
              aes(
                x = pb207206,
                y = pb208206
              ),
              col = "black",
              linetype = "longdash",
              method = "lm",
              fullrange = TRUE) +
  # medians for EF classes
  geom_point(data = nya_pm[,
                      lapply(.SD, median, na.rm = TRUE),
                      .SDcols = c("pb_val", "pb208206",
                                  "pb207206", "al_ef"),
                      by = .(al_ef_bin)][!is.na(al_ef_bin)],
             aes(
               x = pb207206,
               y = pb208206,
               fill = al_ef_bin
             ),
             shape = 23,
             size = 4) +
  # possible end-members
  geom_point(data = endmembers[1:4, ],
             aes(
               x = pb207206,
               y = pb208206
             ),
             size = 4,
             col = "red") +
  # labels for the end-members
  geom_text_repel(data = endmembers[1:4, ],
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           nudge_x = c(-0.005, -0.002, -0.025, -0.005),
                           nudge_y = c(0.03, 0.03, 0.022, -0.025),
                           min.segment.length = unit(0, 'lines'),
                  size = geom.text.size) +
  # adding label on top of regression line
  geom_text(data = data.table(
    pb208206 = c(2.142),
    pb207206 = c(0.8963),
    slope = c(1.09653),
    r2 = c(0.999),
    label = c("through~EF~medians")
  ),
  aes(
    x = pb207206,
    y = pb208206,
    angle = 180 * atan(slope / 2.80) / pi,
    label = paste(label, "~(R^2 ==", r2, ")")
  ),
  size = geom.text.size-0.6444,
  parse = TRUE) +
  # adding explanatory ellipses and points
  geom_ellipse(aes(
    x0 = 0.879,
    y0 = 1.99,
    a = 0.01,
    b = 0.01,
    angle = 0),
    col = NA,
    fill = "gray80",
    alpha = 0.5
  ) +
  geom_ellipse(aes(
    x0 = 0.879,
    y0 = 1.99,
    a = 0.003,
    b = 0.003,
    angle = 0),
    col = "black",
    fill = "gray80") +
  geom_point(aes(
    x = 0.880,
    y = 1.9915
  ),
  col = "black",
  fill = "gray80",
  shape = 23,
  size = 3) +
  # explanatory labels
  geom_text_repel(data = 
      data.frame(pb208206 = c(1.99, 1.99, 1.9915),
                 pb207206 = c(0.871, 0.879, 0.880),
                 label = c("95%-confindence interval\n of the data",
                           "95%-confidence interval\n of the mean",
                           "median")),
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           nudge_x = c(-0.025, 0, +0.01),
                           nudge_y = c(0, 0.025, 0.01),
                           min.segment.length = unit(0, 'lines'),
      size = geom.text.size) +
  # uncertainty on y axis
  geom_errorbar(data = nya_pm[,
                              .(u208206_med = median(u208206, na.rm = TRUE),
                                u207206_med = median(u207206, na.rm = TRUE))],
                aes(x = 0.905,
                    ymin = 1.99 - u208206_med,
                    ymax = 1.99 + u208206_med),
                width = 0.002) +
  # uncertainty on x axis
  geom_errorbarh(data = nya_pm[,
                               .(u208206_med = median(u208206, na.rm = TRUE),
                                 u207206_med = median(u207206, na.rm = TRUE))],
                 aes(y = 1.99,
                     xmin = 0.905 - u207206_med,
                     xmax = 0.905 + u207206_med),
                 height = 0.005) +
  # specifing median uncertainty
  geom_text(
    aes(x = 0.905,
        y = 2.008),
    label = "median\n uncertainty",
    size = geom.text.size) +
  # selecting shapes and labels for data points
  scale_shape_manual(breaks = c("PM10",
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "road"),
                     labels = c(expression("PM"[10]~ "at Ny-Ålesund"),
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "Ny-Ålesund road dust"),
                     values = c(3, 15, 20, 17, 4)) +
  # selecting colours and labels for EF classses
  scale_color_manual(name = expression("PM"[10]~"Enrichment factors"),
                     labels = scales::parse_format(),
                     values = cbPalette[-c(1,2, 6)]) +
  scale_fill_manual(name = expression("PM"[10]~"Enrichment factors"),
                    labels = scales::parse_format(),
                    values = cbPalette[-c(1,2, 6)]) +
  # adjusting the layout
  theme(legend.text.align = 0,
        legend.justification = c(0, 1),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal",
        text = element_text(size = theme.size)) +
  # ordering the legend guides
  guides(shape = guide_legend(order = 1,
                              title = element_blank()), 
         colour = guide_legend(order = 2)) +
  # dixed aspect ratio
  coord_fixed(ratio = 1/2.80,
              xlim = c(0.80, 0.91),
              ylim = c(1.982, 2.185)) +
  # saving the plot
  ggsave(file = "output/figure1.pdf", device = cairo_pdf,
         width = 297, height = 297, unit = "mm",
         scale = 0.9) +
  ggsave(file = "output/figure1.png", dpi = 300,
         width = 297, height = 297, unit = "mm",
         scale = 0.9)

#### Figure 2 ----

# Upper main panel: 207Pb/206Pb vs 208Pb/206Pb sorted by cluster
nya_cluster.plot <- 
  iso208207206_atm +
  # weigheted average for chinese coals
  stat_conf_ellipse(data = chinese_coals,
                    aes(
                      x = pb207206+0.0054, # adjustment for matching Bi 2017
                      y = pb208206+0.0133
                    ),
                    geom = "polygon",
                    fill = "white",
                    col = "black") +
  # PM10 data divided by cluster
  geom_point(data = nya_pm,
             aes(
               x = pb207206,
               y = pb208206,
               col = factor(cluster),
               group = factor(cluster)
             )) +
  # confidence ellipses for clusters
  stat_ellipse(
    data = nya_pm,
    aes(
      x = pb207206,
      y = pb208206,
      col = factor(cluster),
      group = factor(cluster)
    ),
    type = "t",
    geom = "path",
    size = 1) +
  # adjusting color for clusters
  scale_color_manual(name = element_blank(),
                     labels = c("cluster A: Central Asia",
                                "cluster B: North-America"),
                     na.translate = FALSE,
                     drop = TRUE,
                     values = cbPalette[c(8, 4)]) +
  # line by deming regression and u(208Pb/206Pb) = 2 u(207Pb/206)
  geom_smooth(data = nya_pm,
              aes(x = pb207206,
                  y = pb208206),
              method = deming,
              method.args = list(noise_ratio = 2),
              col = "black",
              linetype = "longdash",
              fullrange = TRUE) +
  # Pb isotope ratio values for possible end-members
  geom_point(data = endmembers[c(5,6),],
            aes(
                x = pb207206,
                y = pb208206
                ),
            size = 4,
            col = "red") +
  # labels for the possible end-members
  geom_text_repel(data = endmembers[c(5,6),],
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           nudge_x = c(0.0, -0.02),
                           nudge_y = c(0.07, 0.01),
                           min.segment.length = unit(0, 'lines'),
                  size = geom.text.size) +
  # label for regression line
  geom_text(data = 
              data.table(
                pb208206 = 2.11,
                pb207206 = 0.875,
                slope = 1.9,
                r2 = 0.683,
                label = "through~PM[10]~data"),
  aes(
    x = pb207206,
    y = pb208206,
    angle = 180 * atan(slope / 2.80) / pi,
    label = paste(label, "~(R^2 ==", r2, ")")
  ),
  size = geom.text.size - 0.64444,
  hjust = 0,
  vjust = 0,
  parse = TRUE) +
  # theme layout
  coord_fixed(ratio = 1/2.80,
              xlim = c(0.80, 0.91),
              ylim = c(1.982, 2.185)) +
  theme_bw() +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0),
        text = element_text(size = theme.size))

# upper marginal histogram and density plot
nya207206_cluster <-
  ggplot() +
  # histogram for 207Pb/206Pb values sorted by cluster
  geom_histogram(data = nya_pm[!is.na(cluster)],
                 aes(
                   x = pb207206,
                   fill = factor(cluster),
                   group = factor(cluster),
                   y = ..density..),
                 position = "identity",
                 binwidth = 0.003,
                 alpha = 0.6,
                 col = "white") +
  stat_function(fun = function(x) dnorm(x, 
                                        mean = nya_gmm_parameters$pb207206[1],
                                        sd = nya_gmm_parameters$pb207206[3]),
                n = 1000,
                col = cbPalette[c(8)],
                size = 1) +
  stat_function(fun = function(x) dnorm(x,
                                        mean = nya_gmm_parameters$pb207206[2],
                                        sd = nya_gmm_parameters$pb207206[4]),
                n = 1000,
                col = cbPalette[c(4)],
                size = 1) +
  # segment for the mean value
  geom_segment(
      aes(
    x = nya_gmm_parameters$pb207206[1:2],
    xend = nya_gmm_parameters$pb207206[1:2],
    y = 0,
    yend = nya_gmm_parameters$pb207206[9:10]),
    col = cbPalette[c(8, 4)],
    linetype = "dashed",
    size = 1) +
  # add standard error of the mean
  geom_pointrange(aes(
    xmin = nya_gmm_parameters$pb207206[7:8],
    xmax = nya_gmm_parameters$pb207206[5:6],
    x = nya_gmm_parameters$pb207206[1:2],
    y = 0
  ),
  col = cbPalette[c(8, 4)],
  fatten = 2,
  size = 1) +
  scale_fill_manual(name = element_blank(),
                    guide = NULL,
                    values = cbPalette[c(8, 4)]) +
  scale_x_continuous(
    name = bquote(phantom(.) ^ 207 * Pb ~ "/" * phantom(.) ^ 206 * Pb),
    breaks = seq(0.80, 0.91, by = 0.01),
    limits = c(0.80, 0.91)) +
  theme_void()

# upper right marginal plot 
nya208206_cluster <-
  ggplot() +
  # histogram for 208Pb/206Pb values sorted by cluster
  geom_histogram(data = nya_pm,
                 aes(
                   x = pb208206,
                   fill = factor(cluster),
                   group = factor(cluster),
                   y = ..density..),
                 position = "identity",
                 binwidth = 0.006,
                 alpha = 0.6,
                 col = "white") +
  stat_function(fun = function(x) dnorm(x, 
                                        mean = nya_gmm_parameters$pb208206[1],
                                        sd = nya_gmm_parameters$pb208206[3]),
                n = 1000,
                col = cbPalette[c(8)],
                size = 1) +
  stat_function(fun = function(x) dnorm(x,
                                        mean = nya_gmm_parameters$pb208206[2],
                                        sd = nya_gmm_parameters$pb208206[4]),
                n = 1000,
                col = cbPalette[c(4)],
                size = 1) +
  # segment for the mean value
  geom_segment(
    aes(
      x = nya_gmm_parameters$pb208206[1:2],
      xend = nya_gmm_parameters$pb208206[1:2],
      y = 0,
      yend = nya_gmm_parameters$pb208206[9:10]),
    col = cbPalette[c(8, 4)],
    linetype = "dashed",
    size = 1) +
  # add standard error of the mean
  geom_pointrange(aes(
    xmin = nya_gmm_parameters$pb208206[7:8],
    xmax = nya_gmm_parameters$pb208206[5:6],
    x = nya_gmm_parameters$pb208206[1:2],
    y = 0
  ),
  col = cbPalette[c(8, 4)],
  fatten = 2,
  size = 1) +
  scale_fill_manual(name = element_blank(),
                    guide = NULL,
                    values = cbPalette[c(8, 4)]) +
  scale_x_continuous(
    name = bquote(phantom(.) ^ 208 * Pb ~ "/" * phantom(.) ^ 206 * Pb),
    breaks = seq(1.98, 2.19, by = 0.01),
    limits = c(1.98, 2.19)) +
  theme_void() +
  coord_flip()

# central main panel: 207Pb/206Pb vs log(EF)
nya_cluster.efplot <-
  ggplot() +
  # scatterplot 207Pb/206Pb vs log(EF) for PM10 data sorted by cluster
  geom_point(data = nya_pm,
             aes(
               x = pb207206,
               y = al_ef,
               col = factor(cluster),
               group = factor(cluster)
             )) +
  # ellipses for clusters
  stat_ellipse(
    data = nya_pm,
    aes(
      x = pb207206,
      y = al_ef,
      col = factor(cluster),
      group = factor(cluster)
    ),
    type = "t",
    geom = "path",
    size = 1) +
  # adjust the colors
  scale_color_manual(name = element_blank(),
                     guide = NULL,
                     values = cbPalette[c(8, 4)]) +
  # add log10 y scale
  scale_y_log10(breaks = c(1, 10, 25, 50, 100, 1000),
                limits = c(1, 1000)) +
  annotation_logticks(base = 10, sides = "l") +
  # theme layout 
  scale_x_continuous(breaks = seq(from = 0.7,
                                  to = 1.0,
                                  by = 0.01),
                     limits = c(0.80, 0.91)) +
  labs(x = expression({}^207*"Pb/"*{}^206*"Pb"),
       y = expression(EF(Pb/Al)[c])) +
  coord_fixed(ratio = 0.066666/2.80) +
  theme_bw() +
  theme(text = element_text(size = theme.size))

# bottom main panel: 207Pb/206Pb vs log(Pb concentration)
nya_cluster.pbplot <-
  ggplot() +
  # scatterplot 207Pb/206Pb vs log(Pb) for PM10 data sorted by cluster
  geom_point(data = nya_pm,
             aes(
               x = pb207206,
               y = pb_val,
               col = factor(cluster),
               group = factor(cluster)
             )) +
  # ellipses for the two clusters
  stat_ellipse(
    data = nya_pm,
    aes(
      x = pb207206,
      y = pb_val,
      col = factor(cluster),
      group = factor(cluster)
    ),
    type = "t",
    geom = "path",
    size = 1) +
  # adjust colours
  scale_color_manual(name = element_blank(),
                     guide = NULL,
                     values = cbPalette[c(8, 4)]) +
  # add log10 scale on y axis
  scale_y_log10(
    limits = c(1, 1500)) +
  annotation_logticks(base = 10, sides = "l") +
  # theme layout
  scale_x_continuous(breaks = seq(from = 0.7,
                                  to = 1.0,
                                  by = 0.01),
                     limits = c(0.80, 0.91)) +
  labs(x = expression({}^207*"Pb/"*{}^206*"Pb"),
       y = expression("Pb " * (pg/m^3))) +
  coord_fixed(ratio = 0.06297/2.80) +
  theme_bw() +
  theme(text = element_text(size = theme.size))

# marginal plot for log(EFs)
nya_cluster.efy <-
  ggplot() +
  geom_histogram(data = nya_pm[!is.na(cluster)],
                 aes(
                   x = log10(al_ef),
                   fill = factor(cluster),
                   group = factor(cluster),
                   y = ..density..),
                 position = "identity",
                 binwidth = 0.1,
                 alpha = 0.6,
                 col = "white",
                 na.rm = TRUE) +
  stat_function(fun = function(x) dnorm(x, 
                                        mean = nya_gmm_parameters$log_ef[1],
                                        sd = nya_gmm_parameters$log_ef[3]),
                n = 1000,
                col = cbPalette[c(8)],
                size = 1) +
  stat_function(fun = function(x) dnorm(x,
                                        mean = nya_gmm_parameters$log_ef[2],
                                        sd = nya_gmm_parameters$log_ef[4]),
                n = 1000,
                col = cbPalette[c(4)],
                size = 1) +
  # segment for the mean value
  geom_segment(
    aes(
      x = nya_gmm_parameters$log_ef[1:2],
      xend = nya_gmm_parameters$log_ef[1:2],
      y = 0,
      yend = nya_gmm_parameters$log_ef[9:10]),
    col = cbPalette[c(8, 4)],
    linetype = "dashed",
    size = 1) +
  # add standard error of the mean
  geom_pointrange(aes(
    xmin = nya_gmm_parameters$log_ef[7:8],
    xmax = nya_gmm_parameters$log_ef[5:6],
    x = nya_gmm_parameters$log_ef[1:2],
    y = 0
  ),
  col = cbPalette[c(8, 4)],
  fatten = 2,
  size = 1) +
  scale_fill_manual(name = element_blank(),
                    guide = NULL,
                    values = cbPalette[c(8, 4)]) +
  xlim(c(0, 3)) +
  theme_void() +
  coord_flip()

nya_cluster.pby <-
  ggplot() +
  geom_histogram(data = nya_pm[!is.na(cluster)],
                 aes(
                   x = log10(pb_val),
                   fill = factor(cluster),
                   group = factor(cluster),
                   y = ..density..),
                 position = "identity",
                 binwidth = 0.1,
                 alpha = 0.6,
                 col = "white",
                 na.rm = TRUE) +
  stat_function(fun = function(x) dnorm(x, 
                                        mean = nya_gmm_parameters$log_pb[1],
                                        sd = nya_gmm_parameters$log_pb[3]),
                n = 1000,
                col = cbPalette[c(8)],
                size = 1) +
  stat_function(fun = function(x) dnorm(x,
                                        mean = nya_gmm_parameters$log_pb[2],
                                        sd = nya_gmm_parameters$log_pb[4]),
                n = 1000,
                col = cbPalette[c(4)],
                size = 1) +
  # segment for the mean value
  geom_segment(
    aes(
      x = nya_gmm_parameters$log_pb[1:2],
      xend = nya_gmm_parameters$log_pb[1:2],
      y = 0,
      yend = nya_gmm_parameters$log_pb[9:10]),
    col = cbPalette[c(8, 4)],
    linetype = "dashed",
    size = 1) +
  # add standard error of the mean
  geom_pointrange(aes(
    xmin = nya_gmm_parameters$log_pb[7:8],
    xmax = nya_gmm_parameters$log_pb[5:6],
    x = nya_gmm_parameters$log_pb[1:2],
    y = 0
  ),
  col = cbPalette[c(8, 4)],
  fatten = 2,
  size = 1) +
  scale_fill_manual(name = element_blank(),
                    guide = NULL,
                    values = cbPalette[c(8, 4)]) +
  xlim(c(0, 3.176091)) +
  theme_void() +
  coord_flip()

# arranging the plot in a single page
nya207206_cluster + plot_spacer() + nya_cluster.plot + nya208206_cluster + 
  nya_cluster.efplot + nya_cluster.efy +
  nya_cluster.pbplot + nya_cluster.pby +
  plot_layout(ncol = 2, nrow = 4, widths = c(6, 1.8), heights = c(1, 3, 3, 3))

ggsave(file = "output/figure2.pdf",
       width = 200, 
       height = 297,
       scale = 1.6,
       unit = "mm") +
  ggsave(file = "output/figure2.png",
         width = 200, 
         height = 297,
         dpi = 300,
         scale = 1.6,
         unit = "mm")

#### Figure 3 ----

### time-series for Pb concentrations, EFs and Pb isotope ratio values
# setting the language to English
i1 <- Sys.getlocale()
i2 <- Sys.getenv()

Sys.setlocale("LC_TIME", "en_GB.UTF-8")

# new column of fake sampling dates, all in the same year (2010)
nya_pm[, date_ddmm := as.POSIXct(paste0("2010-", 
                                        format(date, "%m-%d", 
                                               locale = "en_GB.UTF8")))]
# divide the years in group of three
nya_pm[, year_classes := dplyr::case_when(
  year_col %in% c(2010, 2011, 2012) ~ "2010, 2011, 2012",
  year_col %in% c(2013, 2014, 2015) ~ "2013, 2014, 2015",
  year_col %in% c(2016, 2017, 2018) ~ "2016, 2017, 2018")][,
    yaer_classes := factor(year_classes,
                          levels = c("2010, 2011, 2012",
                                     "2013, 2014, 2015",
                                     "2016, 2017, 2018"),
                          labels = c("2010, 2011, 2012",
                                     "2013, 2014, 2015",
                                     "2016, 2017, 2018"),
                          ordered = TRUE)]

### Time-series for Pb concentration for 2010-2012
# the code is commented only for this plot, the following plots are similar
pb_ts_2010_2012 <- ggplot() +
  # adding Pb concentration values sorted by year
  geom_point(data = nya_pm[year_classes == "2010, 2011, 2012"],
             aes(x = date_ddmm,
                 y = pb_val,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  # adding lines by locally estimated scatterplot smoothing (loess)
  geom_smooth(data = nya_pm[year_classes == "2010, 2011, 2012"],
              aes(x = date_ddmm,
                  y = pb_val,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  # adding a loess lines averaged for years from 2010 to 2018
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb_val,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  # annotate the loess line
  annotate(x = as.POSIXct("2010-06-03"),
           xend = as.POSIXct("2010-05-25"),
           y = 20,
           yend = 3,
           geom = "segment") +
  annotate(x = as.POSIXct("2010-05-25"),
           y = 2,
           geom = "text",
           size = geom.text.size,
           label = "2010 – 2018") +
  # log scale on y axis
  scale_y_log10(limits = c(1, 1500),
                breaks = c(1, 10, 100, 1000)) +
  # formatting the month on x axis
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  annotation_logticks(side = "l") +
  # adjusting the colours
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  # divide years by panel
  facet_wrap(vars(year_classes)) +
  # theme layout
  theme_bw(base_size = theme.size) +
  labs(y = expression("[Pb]"~(pg/m^3)),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        legend.key.size = unit(0.8, "line"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for Pb concentration for 2013-2015
pb_ts_2013_2015 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2013, 2014, 2015"],
             aes(x = date_ddmm,
                 y = pb_val,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2013, 2014, 2015"],
              aes(x = date_ddmm,
                  y = pb_val,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb_val,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_log10(limits = c(1, 1500),
                breaks = c(1, 10, 100, 1000)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  annotation_logticks(side = "l") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for Pb concentration for 2016-2018
pb_ts_2016_2018 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2016, 2017, 2018"],
             aes(x = date_ddmm,
                 y = pb_val,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2016, 2017, 2018"],
              aes(x = date_ddmm,
                  y = pb_val,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb_val,
                  group = 1),
              col = "black",
              fill = "black",
              span = 0.8,
              show.legend = FALSE,
              n = 100) +
  scale_y_log10(limits = c(1, 1500),
                breaks = c(1, 10, 100, 1000)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  annotation_logticks(side = "l") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

# aggregate time-series for Pb concentrations
pb_ts <- ggarrange(plotlist = list(pb_ts_2010_2012,
                                   pb_ts_2013_2015,
                                   pb_ts_2016_2018),
                   ncol = 3,
                   widths = c(1, 0.82, 0.82))

### Time-series for EFs for 2010-2012
ef_ts_2010_2012 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2010, 2011, 2012"],
             aes(x = date_ddmm,
                 y = al_ef,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2010, 2011, 2012"],
              aes(x = date_ddmm,
                  y = al_ef,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              span = 0.8,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = al_ef,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              span = 0.8,
              n = 100) +
  geom_hline(yintercept = c(1, 10),
             linetype = "dashed",
             col = "black") +
  scale_y_log10(limits = c(-100, 800),
                breaks = c(1, 10, 25, 50, 100, 200, 400)) +
  annotation_logticks(side = "l") +
  coord_cartesian(ylim = c(1, 500)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = expression(EF(Pb/Al)[c]),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for EFs for 2013-2015
ef_ts_2013_2015 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2013, 2014, 2015"],
             aes(x = date_ddmm,
                 y = al_ef,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = 
                nya_pm[year_classes == "2013, 2014, 2015" &
                       ifelse(month_col %in% c("Jun", "Jul", "Aug"),
                              al_ef < 300, al_ef) ],
              aes(x = date_ddmm,
                  y = al_ef,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              span = 0.8,
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = al_ef,
                  group = 1),
              col = "black",
              fill = "black",
              span = 0.8,
              show.legend = FALSE,
              n = 100) +
  geom_hline(yintercept = c(1, 10),
             linetype = "dashed",
             col = "black") +
  scale_y_log10(limits = c(-100, 800),
                breaks = c(1, 10, 25, 50, 100, 200, 400)) +
  annotation_logticks(side = "l") +
  coord_cartesian(ylim = c(1, 500)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for EFs for 2016-2018
ef_ts_2016_2018 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2016, 2017, 2018"],
             aes(x = date_ddmm,
                 y = al_ef,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2016, 2017, 2018"],
              aes(x = date_ddmm,
                  y = al_ef,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = al_ef,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  geom_hline(yintercept = c(1, 10),
             linetype = "dashed",
             col = "black") +
  scale_y_log10(limits = c(-100, 800),
                breaks = c(1, 10, 25, 50, 100, 200, 400)) +
  annotation_logticks(side = "l") +
  coord_cartesian(ylim = c(1, 500)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Aggregate plot for EFs
ef_ts <- ggarrange(plotlist = list(ef_ts_2010_2012,
                                   ef_ts_2013_2015,
                                   ef_ts_2016_2018),
                   ncol = 3,
                   widths = c(1, 0.82, 0.82))

### Time-series for 208Pb/206Pb for 2010-2012
pb208206_ts_2010_2012 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2010, 2011, 2012"],
             aes(x = date_ddmm,
                 y = pb208206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2010, 2011, 2012"],
              aes(x = date_ddmm,
                  y = pb208206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb208206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_continuous(limits = c(2.045, 2.125),
                     breaks = seq(2.045, 2.125, by = 0.02)) +
  coord_cartesian(ylim = c(2.045, 2.125)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = bquote(phantom(.) ^ 208 * Pb ~ "/" * phantom(.) ^ 206 * Pb),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for 208Pb/206Pb for 2013-2015
pb208206_ts_2013_2015 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2013, 2014, 2015"],
             aes(x = date_ddmm,
                 y = pb208206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2013, 2014, 2015"],
              aes(x = date_ddmm,
                  y = pb208206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb208206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_continuous(limits = c(2.045, 2.125),
                     breaks = seq(2.045, 2.125, by = 0.02)) +
  coord_cartesian(ylim = c(2.045, 2.125)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for 208Pb/206Pb for 2016-2018
pb208206_ts_2016_2018 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2016, 2017, 2018"],
             aes(x = date_ddmm,
                 y = pb208206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2016, 2017, 2018"],
              aes(x = date_ddmm,
                  y = pb208206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              span = 0.8,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb208206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              span = 0.8,
              n = 100) +
  scale_y_continuous(limits = c(2.045, 2.125),
                     breaks = seq(2.045, 2.125, by = 0.02)) +
  coord_cartesian(ylim = c(2.045, 2.125)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

# aggregate plots for 208Pb/206Pb
pb208206_ts <- ggarrange(plotlist = list(pb208206_ts_2010_2012,
                                         pb208206_ts_2013_2015,
                                         pb208206_ts_2016_2018),
                         ncol = 3,
                         widths = c(1, 0.82, 0.82))

### Time-series for 207Pb/206Pb for 2010-2012
pb207206_ts_2010_2012 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2010, 2011, 2012"],
             aes(x = date_ddmm,
                 y = pb207206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2010, 2011, 2012"],
              aes(x = date_ddmm,
                  y = pb207206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb207206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_continuous(limits = c(0.835, 0.880),
                     breaks = seq(0.835, 0.880, by = 0.01)) +
  coord_cartesian(ylim = c(0.835, 0.880)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = bquote(phantom(.) ^ 207 * Pb ~ "/" * phantom(.) ^ 206 * Pb),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for 207Pb/206Pb for 2013-2015
pb207206_ts_2013_2015 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2013, 2014, 2015"],
             aes(x = date_ddmm,
                 y = pb207206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2013, 2014, 2015"],
              aes(x = date_ddmm,
                  y = pb207206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb207206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_continuous(limits = c(0.835, 0.880),
                     breaks = seq(0.835, 0.880, by = 0.01)) +
  coord_cartesian(ylim = c(0.835, 0.880)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Time-series for 207Pb/206Pb for 2016-2018
pb207206_ts_2016_2018 <- ggplot() +
  geom_point(data = nya_pm[year_classes == "2016, 2017, 2018"],
             aes(x = date_ddmm,
                 y = pb207206,
                 col = factor(year_col),
                 fill = factor(year_col),
                 group = factor(year_col))) +
  geom_smooth(data = nya_pm[year_classes == "2016, 2017, 2018"],
              aes(x = date_ddmm,
                  y = pb207206,
                  col = factor(year_col),
                  fill = factor(year_col),
                  group = factor(year_col)),
              level = 0.95,
              n = 100) +
  geom_smooth(data = within(nya_pm, year_classes <- NULL),
              aes(x = date_ddmm,
                  y = pb207206,
                  group = 1),
              col = "black",
              fill = "black",
              show.legend = FALSE,
              n = 100) +
  scale_y_continuous(limits = c(0.835, 0.880),
                     breaks = seq(0.835, 0.880, by = 0.01)) +
  coord_cartesian(ylim = c(0.835, 0.880)) +
  scale_x_datetime(date_labels = "%b",
                   date_breaks = "1 month",
                   date_minor_breaks = "1 week") +
  scale_color_manual(values = c(cbPalette[c(5, 6, 7)])) +
  scale_fill_manual(values = c(cbPalette[c(5, 6, 7)])) +
  facet_wrap(vars(year_classes)) +
  theme_bw(base_size = theme.size) +
  labs(y = element_blank(),
       x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        #legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(ncol = 1))

### Aggregate plots for 207Pb/206Pb
pb207206_ts <- ggarrange(plotlist = list(pb207206_ts_2010_2012,
                                         pb207206_ts_2013_2015,
                                         pb207206_ts_2016_2018),
                         ncol = 3,
                         widths = c(1, 0.82, 0.82))

### Aggregate plots for Pb concentration, EF and isotope ratio values
ggarrange(plotlist = list(pb_ts,
                          ef_ts,
                          pb208206_ts,
                          pb207206_ts),
          ncol = 1) %>% 
### saving the plot
  ggsave(file = "output/figure3.pdf", device = cairo_pdf,
         scale = 1.6,
         height = 297, width = 210, unit = "mm") +
  ggsave(file = "output/figure3.png",
         scale = 1.6, dpi = 300,
         height = 297, width = 210, unit = "mm")

#### Figure 4 ----
## inter-annual variations by season

## Spring
# only this plot is commented, the following plot is very similar
iso_spring <-
  iso208207206_atm +
  # confidence ellipses for PM10 data sorted by year
  stat_ellipse(data = nya_pm[season == "spring"],
               aes(
                 y = pb208206,
                 x = pb207206,
                 fill = as.factor(year_col),
               ),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.5) +
  # confidence ellipse for Kongsjorden sediments
  stat_ellipse(data = mineral[site %like% "sediments"],
               aes(
                 x = pb207206,
                 y = pb208206),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.5) +
  # confidence ellipses of the mean for each year
  stat_conf_ellipse(data = nya_pm[season == "spring"],
                    aes(y = pb208206,
                        x = pb207206,
                        fill = as.factor(year_col)),
                    size = 0.5,
                    col = "black",
                    linetype = "solid",
                    geom = "polygon") +
  # PM10 data points
  geom_point(data = nya_pm[season == "spring"],
             aes(
               x = pb207206,
               y = pb208206,
               shape = "PM10",
             ),
             size = 2,
             col = "gray20") +
  # data points for Kongsfjorden sediments
  geom_point(data = mineral,
             aes(
               x = pb207206,
               y = pb208206,
               shape = site
             ),
             size = 2) +
  # deming regression line for PM10 data
  geom_smooth(data = nya_pm,
              aes(x = pb207206,
                  y = pb208206),
              method = deming,
              method.args = list(noise_ratio = 2),
              col = "black",
              linetype = "longdash",
              fullrange = TRUE) +
  # regression line through median of EF classes
  geom_smooth(data = nya_pm[, lapply(.SD, median, na.rm = TRUE),
                            .SDcols = c("pb_val", "pb208206",
                                        "pb207206", "al_ef"),
                            by = .(al_ef_bin)][
                            !is.na(al_ef_bin)],
              aes(
                x = pb207206,
                y = pb208206
              ),
              col = "black",
              linetype = "longdash",
              method = "lm",
              fullrange = TRUE) +
  # median values for each year
  geom_point(data = nya_pm[season == "spring",
                         lapply(.SD, median, na.rm = TRUE),
                         .SDcols = c("pb_val", "pb208206",
                                     "pb207206", "al_ef"),
                         by = .(year_col)],
           aes(
             x = pb207206,
             y = pb208206,
             fill = as.factor(year_col)
           ),
           shape = 23,
           size = 4) +
  # adding possible end-members
  geom_point(data = endmembers[1:5,],
             aes(
               x = pb207206,
               y = pb208206
             ),
             size = 4,
             col = "red") +
  # labels for the end-members
  geom_text_repel(data = endmembers[1:5,],
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           size = geom.text.size,
                           nudge_x = c(-0.005, -0.002, -0.025, -0.005, 0.0),
                           nudge_y = c(0.03, 0.03, 0.012, -0.025, 0.07),
                           min.segment.length = unit(0, 'lines')) +
  # labels for regression lines
  geom_text(data = data.table(
                    pb208206 = c(2.142, 2.154),
                    pb207206 = c(0.8963, 0.886),
                    slope = c(1.09653, 2.2),
                    r2 = c(0.999, 0.683),
                    label = c("through~EF~medians",
                              "through~PM[10]~data")),
            aes(
              x = pb207206,
              y = pb208206,
              angle = 180 * atan(slope / 2.80) / pi,
              label = paste(label, "~(R^2 ==", r2, ")")
                ),
             parse = TRUE,
            size = geom.text.size - 0.64444) +
  # adjusting symbols and colours
  scale_shape_manual(breaks = c("PM10",
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "road"),
                     labels = c(expression("PM"[10]~ "at Ny-Ålesund"),
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "Ny-Ålesund road dust"),
                     values = c(3, 15, 20, 17, 4)) +
  scale_color_manual(name = element_blank(),
              values = RColorBrewer::brewer.pal(name = "Paired", n = 9)) +
  scale_fill_manual(name = element_blank(),
              values = RColorBrewer::brewer.pal(name = "Paired", n = 9)) +
  # annotate the plot with the season
  labs(title = "Spring") +
  # theme layout
  theme(legend.text.align = 0,
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99)) +
  # sorting the legend guides
  guides(shape = guide_legend(order = 1,
                              title = element_blank()),
         color = FALSE,
         fill = FALSE) +
  # fixed aspect ratio for the plot
  coord_fixed(ratio = 1/2.80,
              xlim = c(0.80, 0.91),
              ylim = c(1.982, 2.185))

# Summer
iso_summer <-
  iso208207206_atm +
  stat_ellipse(data = nya_pm[season == "summer"],
               aes(
                 y = pb208206,
                 x = pb207206,
                 fill = as.factor(year_col),
               ),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.5) +
  stat_ellipse(data = mineral[site %like% "sediments"],
               aes(
                 x = pb207206,
                 y = pb208206),
               geom = "polygon",
               linetype = "solid",
               alpha = 0.5) +
  stat_conf_ellipse(data = nya_pm[season == "summer"],
                    aes(y = pb208206,
                        x = pb207206,
                        fill = as.factor(year_col)),
                    size = 0.5,
                    col = "black",
                    linetype = "solid",
                    geom = "polygon") +
  geom_point(data = nya_pm[season == "summer"],
             aes(
               x = pb207206,
               y = pb208206,
               shape = "PM10"
             ),
             size = 2,
             col = "gray20") +
  geom_point(data = mineral,
             aes(
               x = pb207206,
               y = pb208206,
               shape = site
             ), size = 2) +
  geom_smooth(data = nya_pm,
              aes(x = pb207206,
                  y = pb208206),
              method = deming,
              method.args = list(noise_ratio = 2),
              col = "black",
              linetype = "longdash",
              fullrange = TRUE) +
  geom_smooth(data = nya_pm[, lapply(.SD, median, na.rm = TRUE),
                            .SDcols = c("pb_val", "pb208206",
                                        "pb207206", "al_ef"),
                            by = .(al_ef_bin)][
                              !is.na(al_ef_bin)],
              aes(
                x = pb207206,
                y = pb208206
              ),
              col = "black",
              linetype = "longdash",
              method = "lm",
              fullrange = TRUE) +
geom_point(data = nya_pm[season == "summer", 
                         lapply(.SD, median, na.rm = TRUE),
                         .SDcols = c("pb_val", "pb208206",
                                     "pb207206", "al_ef"),
                         by = .(year_col)],
           aes(
             x = pb207206,
             y = pb208206,
             fill = as.factor(year_col)
           ),
           shape = 23,
           size = 4) +
  geom_point(data = endmembers[1:5, ],
             aes(
               x = pb207206,
               y = pb208206
             ),
             size = 4,
             col = "red") +
  geom_text_repel(data = endmembers[1:5,],
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           size = geom.text.size,
                           nudge_x = c(-0.005, -0.002, -0.025, -0.005, 0.0),
                           nudge_y = c(0.03, 0.03, 0.012, -0.025, 0.07),
                           min.segment.length = unit(0, 'lines')) +
  geom_text(data = data.table(
                    pb208206 = c(2.142, 2.154),
                    pb207206 = c(0.8963, 0.886),
                    slope = c(1.09653, 2.2),
                    r2 = c(0.999, 0.683),
                    label = c("through~EF~medians",
                              "through~PM[10]~data")),
            aes(
              x = pb207206,
              y = pb208206,
              angle = 180 * atan(slope / 2.80) / pi,
              label = paste(label, "~(R^2 ==", r2, ")")
            ),
            parse = TRUE,
            size = geom.text.size - 0.6444) +
  # adding explanatory ellipses and labelled points
  geom_ellipse(aes(
    x0 = 0.879,
    y0 = 1.99,
    a = 0.01,
    b = 0.01,
    angle = 0),
    col = NA,
    fill = "gray80",
    alpha = 0.5
  ) +
  geom_ellipse(aes(
    x0 = 0.879,
    y0 = 1.99,
    a = 0.003,
    b = 0.003,
    angle = 0),
    col = "black",
    fill = "gray80") +
  geom_point(aes(
    x = 0.880,
    y = 1.9915
  ),
  col = "black",
  fill = "gray80",
  shape = 23,
  size = 3) +
  geom_text_repel(data = data.frame(
    pb208206 = c(1.99, 1.99, 1.9915),
    pb207206 = c(0.871, 0.879, 0.880),
    label = c("95%-confindence interval\n of the data",
              "95%-confidence interval\n of the mean",
               "median")),
                           aes(
                             x = pb207206,
                             y = pb208206,
                             label = label),
                           size = geom.text.size,
                           nudge_x = c(-0.025, 0, +0.01),
                           nudge_y = c(0, 0.025, 0.01),
                           min.segment.length = unit(0, 'lines')) +
  geom_errorbar(data = nya_pm[,
                        .(u208206_med = median(u208206, na.rm = TRUE),
                          u207206_med = median(u207206, na.rm = TRUE))],
                aes(x = 0.905,
                    ymin = 1.99 - u208206_med,
                    ymax = 1.99 + u208206_med),
                width = 0.002) +
  geom_errorbarh(data = nya_pm[,
                          .(u208206_med = median(u208206, na.rm = TRUE),
                            u207206_med = median(u207206, na.rm = TRUE))],
                 aes(y = 1.99,
                     xmin = 0.905 - u207206_med,
                     xmax = 0.905 + u207206_med),
                 height = 0.005) +
  geom_text(
    aes(x = 0.905,
        y = 2.008),
    size = geom.text.size,
    label = "median\n uncertainty") +
  scale_shape_manual(breaks = c("PM10",
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "road"),
                     labels = c(expression("PM"[10]~ "at Ny-Ålesund"),
                                "Kongsfjorden sediments",
                                "Gruvebadet Laboratory",
                                "Zeppelin Observatory",
                                "Ny-Ålesund road dust"),
                     values = c(3, 15, 20, 17, 4)) +
  scale_color_manual(name = element_blank(),
              values = RColorBrewer::brewer.pal(name = "Paired", n = 9)) +
  scale_fill_manual(name = element_blank(),
              values = RColorBrewer::brewer.pal(name = "Paired", n = 9)[-7]) +
  labs(title = "Summer") +
  theme(legend.text.align = 0,
        #legend.title = element_text(size = 11), 
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99)) +
  guides(colour = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 2),
         shape = FALSE) +
  coord_fixed(ratio = 1/2.80,
              xlim = c(0.80, 0.91),
              ylim = c(1.982, 2.185))

# saving the plot
ggsave(file = "output/figure4.pdf",
       gridExtra::arrangeGrob(iso_spring, iso_summer),
       device = cairo_pdf,
       width = 210, 
       height = 297,
       scale = 1.3,
       unit = "mm") +
  ggsave(file = "output/figure4.png",
         gridExtra::arrangeGrob(iso_spring, iso_summer),
         width = 210, 
         height = 297,
         dpi = 300,
         scale = 1.3,
         unit = "mm")

#### Figure 5 ----
# Temporal trend for relative importance of the two cluster and their effect
# Pb concentration during spring and summer from 2010 to 2018

# Contribution to Pb concentration from each cluster
nya_sum <- nya_pm[al_ef > 10 & !is.na(cluster),
                  .(pb = sum(pb_val, na.rm = TRUE)),
                  by = .(season, year_col, cluster)]

nya_n <- nya_pm[al_ef > 10 & !is.na(cluster),
                .(n = .N),
                by = .(season, year_col)]

nya_pb_cluster <- merge(nya_sum, nya_n)[, pb_m := pb/n]

# Relative contribution to Pb concentration
nya_pb_cluster[, pb_prop := pb_m / sum(pb_m), by = .(season, year_col)]

# Relative contribution from each cluster considering sample numbers
nya_cluster_prop  <- nya_pm[al_ef > 10 & !is.na(cluster),
                                .(year_col, cluster, season)] %>% 
                                table %>% 
                                prop.table(margin = c(1, 3)) %>% 
                                data.table %>%
                                .[!is.nan(N) & N != 0, 
                                  .(year_col = as.numeric(year_col), 
                                    cluster, 
                                    season, 
                                    cluster_prop = N)]

nya_pb_cluster <- merge(nya_pb_cluster,
         nya_cluster_prop, 
         by = c("season", "year_col", "cluster"))

## Cluster relative contribution
nya_prop_plot <-

  # proportion data by year and season
  ggplot(data = nya_pb_cluster,
                  aes(x = interaction(year_col, season, 
                                      lex.order = TRUE, drop = TRUE),
                      y = cluster_prop,
                      fill = factor(cluster))) +
  # adjust colors
  geom_col(position = "fill",
           width = 1,
           col = "black") +
  # divide by year
  facet_grid(.~ year_col, space = 'free_x', scales = 'free_x', switch = 'x') +
  # extract season from interaction of year and season
  scale_x_discrete(labels = function(x) substr(x, 6, 11),
                   expand = c(0, 0)) +
  # adjust colors
  scale_fill_manual(labels = c("Central Asia", "North America"),
                    values = cbPalette[c(8, 4)],
                    name = NULL) +
  # theme layout
  scale_y_continuous(name = "Cluster relative contribution (%)",
                     breaks = seq(0, 1, by = 0.1)) +
  labs(x = "") +
  theme_bw(base_size = theme.size) +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  # remove facet spacing on x-direction
  theme(panel.spacing.x = unit(0,"line")) +
  # switch the facet strip label to outside 
  # remove background color
  theme(strip.placement = 'outside',
        strip.background.x = element_blank())

nya_conc_plot <-

  ggplot(data = nya_pb_cluster,
    aes(x = interaction(year_col, season,
                             lex.order = TRUE, drop = TRUE),
             y = pb_m,
             fill = factor(cluster))) +
  geom_col(position = "stack",
           width = 1,
           col = "black") +
  facet_grid(.~ year_col, space = 'free_x', scales = 'free_x', switch = 'x') +
  scale_x_discrete(labels = function(x) substr(x, 6, 11),
                   expand = c(0, 0)) +
  scale_fill_manual(labels = c("Central Asia", "North America"),
                    values = cbPalette[c(8, 4)],
                    name = NULL) +
  scale_y_continuous(name = "Pb contribution"~(pg/m^3),
                     breaks = seq(0, 300, 20)) +
  labs(x = "") +
  theme_bw(base_size = theme.size) +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  # remove facet spacing on x-direction
  theme(panel.spacing.x = unit(0,"line")) +
  # switch the facet strip label to outside 
  # remove background color
  theme(strip.placement = 'outside',
        strip.background.x = element_blank())

ggarrange(nya_prop_plot, nya_conc_plot,
          nrow = 2,
          labels = c("a", "b"),
          label.y = 0.15,
          legend = "top",
          common.legend = TRUE) +
  ggsave(file = "output/figure5.pdf",
         device = cairo_pdf,
         width = 200, 
         height = 297,
         scale = 1.2,
         unit = "mm") +
  ggsave(file = "output/figure5.png",
         width = 200, 
         height = 297,
         scale = 1.2,
         dpi = 300,
         unit = "mm")

#### Figure S1 ----
# Comparison of Pb isotope ratio values for PM10 data with Pb signatures for
# Chinese coals

ss1_plot <- 
  ggplot() +
  # Ny-Ålesund PM10 Pb isotope ratio values
  geom_point(data = nya_pm,
             aes(
               x = pb207206,
               y = pb208206,
               col = factor(cluster),
               group = factor(cluster)
             )) +
  # confidence ellipses for the mean of Pb isotope ratio values of Chinese
  # coals from different regions (Bi et al., 2017)
  stat_conf_ellipse(data = chinese_coals[, n := .N, by = region][n >= 8],
                    aes(x = pb207206,
                        y = pb208206,
                        group = region),
                    geom = "polygon",
                    fill = "white",
                    col = "black") +
  # confidence ellipse for Pb isotope ratio values of PM10 data from
  # Central Asian and North American sources
  stat_ellipse(data = nya_pm,
               aes(
                x = pb207206,
                y = pb208206,
                col = factor(cluster),
                group = factor(cluster)
                ),
               type = "t",
               geom = "path",
               size = 1) +
  # ajusting colors and names
  scale_color_manual(name = element_blank(),
                     labels = c("cluster A: Central Asia",
                                "cluster B: North-America"),
                     na.translate = FALSE,
                     drop = TRUE,
                     values = cbPalette[c(8, 4)]) +
  # deming regression line for PM10 Pb isotope ratio values
  geom_smooth(data = nya_pm,
              aes(x = pb207206,
                  y = pb208206),
              method = deming,
              method.args = list(noise_ratio = 2),
              col = "black",
              linetype = "longdash",
              fullrange = TRUE) +
  # labels for Chinese coals
  geom_label_repel(data = 
                        chinese_coals[n >= 8,
                                     lapply(.SD, median, na.rm = TRUE),
                                     .SDcols = c("pb208206", "pb207206"),
                                     by = "region"],
                            aes(
                              x = pb207206,
                              y = pb208206,
                              label = region
                            ),
                            size = geom.text.size,
                            box.padding = 1,
                            max.overlaps = Inf) +
  # label for the regression line
  geom_text(data = 
              data.table(
                pb208206 = 2.11,
                pb207206 = 0.875,
                slope = 1.9,
                r2 = 0.683,
                label = "through~PM[10]~data"),
            aes(
              x = pb207206,
              y = pb208206,
              angle = 180 * atan(slope / 2.80) / pi,
              label = paste(label, "~(R^2 ==", r2, ")")
            ),
            size = geom.text.size - 0.64444,
            hjust = 0,
            vjust = 0,
            parse = TRUE) +
  # theme layout
  scale_y_continuous(breaks = seq(from = 1.7,
                                  to = 2.3,
                                  by = 0.02),
                     limits = c(1.7, 2.3)) +
  scale_x_continuous(breaks = seq(from = 0.7,
                                  to = 1.0,
                                  by = 0.01),
                     limits = c(0.7, 1.0)) +
  labs(x = expression({}^207*"Pb/"*{}^206*"Pb"),
       y = expression({}^208*"Pb/"*{}^206*"Pb"),
       col = element_blank()) +
  coord_fixed(ratio = 1/2.80,
              xlim = c(0.80, 0.91),
              ylim = c(1.982, 2.185)) +
  theme_bw(base_size = theme.size) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1)) +
  # save the plot
  ggsave(file = "output/figure_s1.pdf",
         device = cairo_pdf,
         width = 297, 
         height = 200,
         scale = 1,
         unit = "mm") +
  ggsave(file = "output/figure_s1.png",
         width = 297, 
         height = 200,
         scale = 1,
         dpi = 300,
         unit = "mm")
############################### The end #######################################
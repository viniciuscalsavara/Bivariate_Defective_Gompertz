## Full application analysis workflow used in the manuscript
## This script reproduces: descriptive summaries, sex-only model,
## sex+age model, independence test, and copula model comparison.

## Packages ####################################
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)
library(mvtnorm)
library(patchwork)
library(boot)




## functions ##########
source("functions.R")

## Data preparation used across all application analyses
dados <- read.table("dados_tmo.txt", header = TRUE)
dados <- clean_names(dados)
dados$sexo1 <- ifelse(dados$sexo == 1, 0, 1)
dados$sex <- factor(dados$sexo1)
dados$age <- dados$idade
dados$tempdeag <- dados$tempdeag / 30.4
dados$tempdecr <- dados$tempdecr / 30.4

## Descriptive summaries used in the manuscript application section
n_total <- nrow(dados)
n_both <- sum(dados$deag == 1 & dados$decr == 1)
n_only_t1 <- sum(dados$deag == 1 & dados$decr == 0)
n_only_t2 <- sum(dados$deag == 0 & dados$decr == 1)
n_none <- sum(dados$deag == 0 & dados$decr == 0)

complete_idx <- which(dados$deag == 1 & dados$decr == 1)
spearman_complete <- suppressWarnings(cor(
  dados$tempdeag[complete_idx],
  dados$tempdecr[complete_idx],
  method = "spearman"
))

desc_tbl <- data.frame(
  metric = c("n_total", "n_both", "n_only_t1", "n_only_t2", "n_none", "spearman_complete_cases"),
  value = c(n_total, n_both, n_only_t1, n_only_t2, n_none, spearman_complete)
)
write.csv(desc_tbl,  "descriptive_summary.csv", row.names = FALSE)

# Baseline covariates summary (sex and age)
n_total <- nrow(dados)
n_male <- sum(dados$sex == 1, na.rm = TRUE)
pct_male <- 100 * n_male / n_total

age_mean <- mean(dados$age, na.rm = TRUE)
age_sd <- sd(dados$age, na.rm = TRUE)
age_min <- min(dados$age, na.rm = TRUE)
age_max <- max(dados$age, na.rm = TRUE)

cat("Total patients:", n_total, "\n")
cat("Male:", n_male, sprintf("(%.1f%%)", pct_male), "\n")
cat("Age mean (SD):", sprintf("%.1f (%.1f)", age_mean, age_sd), "\n")
cat("Age range:", sprintf("%d-%d", age_min, age_max), "\n\n")


# Median follow-up with 95% CI for each outcome
# (KM of follow-up time using reverse censoring)
# T1: event indicator = 1 - deag
fit_fu_t1 <- survfit(Surv(tempdeag, 1 - deag) ~ 1, data = dados)
tab_t1 <- summary(fit_fu_t1)$table
med_fu_t1 <- unname(tab_t1["median"])
lcl_fu_t1 <- unname(tab_t1["0.95LCL"])
ucl_fu_t1 <- unname(tab_t1["0.95UCL"])

# T2: event indicator = 1 - decr
fit_fu_t2 <- survfit(Surv(tempdecr, 1 - decr) ~ 1, data = dados)
tab_t2 <- summary(fit_fu_t2)$table
med_fu_t2 <- unname(tab_t2["median"])
lcl_fu_t2 <- unname(tab_t2["0.95LCL"])
ucl_fu_t2 <- unname(tab_t2["0.95UCL"])

cat("Median follow-up T1 (months):",
    sprintf("%.2f (95%%CI: %.2f to %.2f)", med_fu_t1, lcl_fu_t1, ucl_fu_t1), "\n")
cat("Median follow-up T2 (months):",
    sprintf("%.2f (95%%CI: %.2f to %.2f)", med_fu_t2, lcl_fu_t2, ucl_fu_t2), "\n\n")


# Figure 1: Kaplan-Meier curves by sex for acute and chronic GVHD

dados <- dados |>
  mutate( 
    sexo1 = factor(sex, levels = c(0, 1), labels = c("Female", "Male"))
  )

fit_t1 <- survfit(Surv(tempdeag, deag) ~ sexo1, data = dados)
fit_t2 <- survfit(Surv(tempdecr, decr) ~ sexo1, data = dados)

p1 <- ggsurvplot(
  fit_t1,
  data = dados,
  conf.int = FALSE,
  pval = TRUE,
  pval.coord = c(0, 5),
  pval.method = FALSE,
  censor = FALSE,
  risk.table = FALSE,
  legend.title = "Sex",
  legend.labs = c("Female", "Male"),
  xlab = "Time",
  ylab = "Survival probability (%)",
  fun = "pct",
  break.time.by = 25,
  palette = c("#1F77B4", "#D62728"),
  ggtheme = theme_minimal()
)$plot + ggtitle("Acute GVHD")

p2 <- ggsurvplot(
  fit_t2,
  data = dados,
  conf.int = FALSE,
  pval = TRUE,
  pval.coord = c(0, 5),
  pval.method = FALSE,
  censor = FALSE,
  risk.table = FALSE,
  legend.title = "Sex",
  legend.labs = c("Female", "Male"),
  xlab = "Time",
  ylab = "Survival probability (%)",
  fun = "pct",
  break.time.by = 25,
  palette = c("#1F77B4", "#D62728"),
  ggtheme = theme_minimal()
)$plot + ggtitle("Chronic GVHD")

fig1 <- p1 + p2 + plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

print(fig1)

ggsave(
  filename = "figure1_km_by_sex.pdf",
  plot = fig1,
  width = 11,
  height = 4.8
)

#################### Models fit ################################################
## Sex-only model (FGM defective Gompertz) for Tables/Figures in Section 7.1
form_a1_sex <- tempdeag ~ sex
form_b1_sex <- tempdeag ~ sex
form_a2_sex <- tempdecr ~ sex
form_b2_sex <- tempdecr ~ sex

fit_fgm_sex <- ajuste(
  linear_pred_a1 = form_a1_sex,
  linear_pred_b1 = form_b1_sex,
  linear_pred_a2 = form_a2_sex,
  linear_pred_b2 = form_b2_sex,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados
)

write.csv(fit_fgm_sex$estimates, "fgm_sex_estimates.csv", row.names = FALSE)
write.csv(fit_fgm_sex$fracoes_cura, "fgm_sex_cure_fractions.csv", row.names = FALSE)

ggsave("figure5_survival_t1_sex.pdf", fit_fgm_sex$plot_t1, width = 7, height = 5)
ggsave("figure6_survival_t2_sex.pdf", fit_fgm_sex$plot_t2, width = 7, height = 5)

## Sex+age model (FGM defective Gompertz) for Section 7.2
form_a1_age <- tempdeag ~ sex + age
form_b1_age <- tempdeag ~ sex + age
form_a2_age <- tempdecr ~ sex + age
form_b2_age <- tempdecr ~ sex + age

fit_fgm_age <- ajuste_covar_cont(
  linear_pred_a1 = form_a1_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b2 = form_b2_age,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados
)

frac_age <- fracoes_covar_cont(
  linear_pred_a1 = form_a1_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b2 = form_b2_age,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados,
  var_cont = "age",
  var_cat = "sex",
  fit_par = fit_fgm_age$estimates$emv[-c(length(fit_fgm_age$estimates$emv) - 1, length(fit_fgm_age$estimates$emv))],
  var_aux_covar = fit_fgm_age$var_covar,
)

write.csv(fit_fgm_age$estimates, "fgm_sex_age_estimates.csv", row.names = FALSE)
write.csv(frac_age$fracoes_cura, "fgm_sex_age_cure_fractions.csv", row.names = FALSE)

## Likelihood ratio test for independence (phi = 0)
fit_fgm_age0 <- ajuste_covar_cont_phi0(
  linear_pred_a1 = form_a1_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b2 = form_b2_age,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados
)

par_fgm <- fit_fgm_age$estimates$emv[1:13]
llik_fgm <- -lik_dg_fgm_mv(
  par = par_fgm,
  data = dados,
  cens1 = "deag",
  cens2 = "decr",
  linear_pred_a1 = form_a1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_b2 = form_b2_age,
  eps = 1e-12
)

llik_phi0 <- -lik_dg_phi0_mv(
  par = fit_fgm_age0$estimates$emv,
  data = dados,
  cens1 = "deag",
  cens2 = "decr",
  linear_pred_a1 = form_a1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_b2 = form_b2_age,
  eps = 1e-12
)

lrt_phi <- 2 * (llik_fgm - llik_phi0)
pval_phi <- pchisq(lrt_phi, df = 1, lower.tail = FALSE)

lrt_tbl <- data.frame(
  test = "H0: phi = 0",
  logLik_full = llik_fgm,
  logLik_null = llik_phi0,
  LR = lrt_phi,
  df = 1,
  p_value = pval_phi
)
write.csv(lrt_tbl, "lrt_independence_phi.csv", row.names = FALSE)


## Frank model (sex+age) used in model comparison
fit_frank_age <- ajuste_frank_covar_cont(
  linear_pred_a1 = form_a1_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b2 = form_b2_age,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados
)

par_frank <- head(fit_frank_age$estimates$emv, -1)
llik_frank <- - lik_dg_frank_mv(
  par = par_frank,
  data = dados,
  cens1 = "deag",
  cens2 = "decr",
  linear_pred_a1 = form_a1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_b2 = form_b2_age,
  eps = 1e-12
)

## Gaussian model (sex+age) used in model comparison
fit_gaussian_age <- ajuste_gaussian_covar_cont(
  linear_pred_a1 = form_a1_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b2 = form_b2_age,
  cens1 = "deag",
  cens2 = "decr",
  dados = dados
)

npar_gauss <- ncol(model.matrix(form_a1_age, dados)) +
  ncol(model.matrix(form_b1_age, dados)) +
  ncol(model.matrix(form_a2_age, dados)) +
  ncol(model.matrix(form_b2_age, dados)) + 1

par_gauss <- fit_gaussian_age$estimates$emv[1:npar_gauss]
llik_gauss <- -lik_dg_gaussian_mv(
  par = par_gauss,
  data = dados,
  cens1 = "deag",
  cens2 = "decr",
  linear_pred_a1 = form_a1_age,
  linear_pred_a2 = form_a2_age,
  linear_pred_b1 = form_b1_age,
  linear_pred_b2 = form_b2_age,
  eps = 1e-12
)

## AIC and AICc table for FGM, Frank, and Gaussian copulas
k_fgm <- length(par_fgm)
k_frank <- length(par_frank)
k_gauss <- length(par_gauss)
n_obs <- nrow(dados)

cmp_tbl <- data.frame(
  model = c("Defective Gompertz + FGM copula", "Defective Gompertz + Frank copula", "Defective Gompertz + Gaussian copula"),
  AIC = c(AIC_calc(llik_fgm, k_fgm), AIC_calc(llik_frank, k_frank), AIC_calc(llik_gauss, k_gauss)),
  AICc = c(AICc_calc(llik_fgm, k_fgm, n_obs), AICc_calc(llik_frank, k_frank, n_obs), AICc_calc(llik_gauss, k_gauss, n_obs))
)

write.csv(cmp_tbl, "model_comparison_aic_aicc.csv", row.names = FALSE)




#########

dados$tempplaq_meses <- dados$tempplaq / 30.4
dados$tempdeag_meses <- dados$tempdeag / 30.4
dados$tempdecr_meses <- dados$tempdecr / 30.4


dados <- dados %>%
  mutate(
    grupo = case_when(
      deag == 1 & decr == 0 ~ "Only acute GVHD",
      deag == 0 & decr == 1 ~ "Only chronic GVHD",
      deag == 1 & decr == 1 ~ "Both",
      deag == 0 & decr == 0 ~ "None"
    )
  )


p<-ggplot(dados, aes(x = tempdeag_meses, y = tempdecr_meses)) +
  geom_point(aes(shape = grupo, color = grupo), size = 3) +
  labs(
    x = "Time to acute GVHD (months)",
    y = "Time to chronic GVHD (months)",
    shape = "Event",
    color = "Event"
  ) +
  scale_x_log10(breaks = c(0.5,2,5,20,40,100)) +
  scale_y_log10(breaks = c(0.5,2,5,20,40,100))+
  #geom_abline(linetype = 2)+
  theme_minimal(base_size = 20)+
  theme(
    legend.position = "bottom"
  )


p


ggsave(
  filename = "scatter_bivariate_log10.pdf",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)



##Correlation: Spearman

cor.test(
  dados$tempdeag_meses,
  dados$tempdecr_meses,
  method = "spearman"
)


#Only events

dados_eventos <- subset(dados, deag == 1 & decr == 1)

cor.test(
  dados_eventos$tempdeag_meses,
  dados_eventos$tempdecr_meses,
  method = "spearman"
)






# bootstrap
spearman_boot <- function(data, indices){
  d <- data[indices, ]
  return(cor(d$tempdeag_meses, d$tempdecr_meses, method = "spearman"))
}


set.seed(123)

boot_res <- boot(
  data = dados_eventos,
  statistic = spearman_boot,
  R = 2000
)


boot.ci(boot_res, type = "perc")





### Hazard function and Hazard ratio - Figure 7


#===========================================================
# Time grid
#===========================================================
t_grid <- seq(0, 130, 0.01)

#===========================================================
# Application results
# 0 = Female
# 1 = Male
#===========================================================
scenario <- list(
  
  list(
    outcome = "Acute GVHD",
    beta0   = -1.2404,
    beta1   =  0.4399,
    alpha0  = -0.7077,
    alpha1  = -0.1669,
    se_beta0  = 0.3139,
    se_beta1  = 0.4621,
    se_alpha0 = 0.1840,
    se_alpha1 = 0.3120
  ),
  
  list(
    outcome = "Chronic GVHD",
    beta0   = -2.3630,
    beta1   = -0.0504,
    alpha0  = -0.0560,
    alpha1  = -0.0025,
    se_beta0  = 0.2226,
    se_beta1  = 0.3568,
    se_alpha0 = 0.0195,
    se_alpha1 = 0.0293
  )
)

#===========================================================
# Functions
#===========================================================

# log-hazard
loghaz_gomp <- function(t, x, beta0, beta1, alpha0, alpha1) {
  beta0 + beta1 * x + (alpha0 + alpha1 * x) * t
}

# hazard
haz_gomp <- function(t, x, beta0, beta1, alpha0, alpha1) {
  exp(loghaz_gomp(t, x, beta0, beta1, alpha0, alpha1))
}

# pointwise SE for log-hazard
# ignoring covariances
se_loghaz_gomp <- function(t, x, se_beta0, se_beta1, se_alpha0, se_alpha1) {
  sqrt(
    se_beta0^2 +
      (x^2) * se_beta1^2 +
      (t^2) * se_alpha0^2 +
      (x^2) * (t^2) * se_alpha1^2
  )
}

# HR: male vs female
loghr_gomp <- function(t, beta1, alpha1) {
  beta1 + alpha1 * t
}

hr_gomp <- function(t, beta1, alpha1) {
  exp(loghr_gomp(t, beta1, alpha1))
}

# pointwise SE for log-HR
# ignoring covariance between beta1 and alpha1
se_loghr_gomp <- function(t, se_beta1, se_alpha1) {
  sqrt(se_beta1^2 + (t^2) * se_alpha1^2)
}

#===========================================================
# Build hazard data with 95% CI
#===========================================================

haz_dat <- do.call(
  rbind,
  lapply(scenario, function(s) {
    
    do.call(
      rbind,
      lapply(c(0, 1), function(x) {
        
        group_label <- ifelse(x == 0, "Female", "Male")
        
        logh <- loghaz_gomp(
          t = t_grid, x = x,
          beta0 = s$beta0, beta1 = s$beta1,
          alpha0 = s$alpha0, alpha1 = s$alpha1
        )
        
        se_logh <- se_loghaz_gomp(
          t = t_grid, x = x,
          se_beta0 = s$se_beta0, se_beta1 = s$se_beta1,
          se_alpha0 = s$se_alpha0, se_alpha1 = s$se_alpha1
        )
        
        data.frame(
          t = t_grid,
          y = exp(logh),
          lower = exp(logh - 1.96 * se_logh),
          upper = exp(logh + 1.96 * se_logh),
          Group = group_label,
          Outcome = s$outcome
        )
      })
    )
  })
)

#===========================================================
# Build HR data with 95% CI
#===========================================================

hr_dat <- do.call(
  rbind,
  lapply(scenario, function(s) {
    
    loghr <- loghr_gomp(t_grid, s$beta1, s$alpha1)
    se_loghr <- se_loghr_gomp(t_grid, s$se_beta1, s$se_alpha1)
    
    data.frame(
      t = t_grid,
      y = exp(loghr),
      lower = exp(loghr - 1.96 * se_loghr),
      upper = exp(loghr + 1.96 * se_loghr),
      Outcome = s$outcome
    )
  })
)



#===========================================================
# Plot hazard with pointwise 95% CI
#===========================================================

p_haz <- ggplot(haz_dat, aes(x = t, y = y, linetype = Group, fill = Group)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, color = NA) +
  scale_color_manual(
    values = c("#00468B", "#ED0000")
  )+
  geom_line(linewidth = 1) +
  facet_wrap(~Outcome, scales = "free_y") +
  #scale_y_continuous(
  #  limits = c(0, 1.6),
  #  breaks = seq(0, 1.6, by = 0.25)
  #)+
  scale_x_sqrt(
    breaks = c(0,5, 25, 50, 75, 100, 125)
  ) +
  labs(
    x = "Time (months)",
    y = "Hazard function",
    linetype = "Sex",
    fill = "Sex"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )



p_haz <- ggplot(
  haz_dat,
  aes(x = t, y = y, color = Group, linetype = Group)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Outcome, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Female" = "#00468B",
      "Male"   = "#ED0000"
    )
  ) +
  scale_x_sqrt(
    breaks = c(0, 5, 25, 50, 75, 100, 125)
  ) +
  labs(
    x = "Time (months)",
    y = "Hazard function",
    color = "Sex",
    linetype = "Sex"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )

p_haz


#===========================================================
# Plot HR with pointwise 95% CI
#===========================================================
# 
# p_hr <- ggplot(hr_dat, aes(x = t, y = y)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
#   geom_line(linewidth = 1) +
#   geom_hline(yintercept = 1, linetype = 3) +
#   facet_wrap(~Outcome, scales = "free_y") +
#   theme_bw(base_size = 18) +
#   theme(
#     strip.background = element_blank(),
#     legend.position = "none"
#   ) +
#   labs(
#     x = "Time (months)",
#     y = "Hazard ratio (Male vs Female)"
#   ) +
#   theme_bw(base_size = 18) +
#   theme(
#     strip.background = element_blank(),
#     legend.position = "none"
#   )
# 
# p_hr

p_hr <- ggplot(hr_dat, aes(x = t, y = y, color = Outcome)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = 3, color = "black") +
  facet_wrap(~Outcome) +
  scale_y_continuous(
    limits = c(0, 1.6),
    breaks = seq(0, 1.6, by = 0.25)
  )+
  scale_x_sqrt(
    breaks = c(0,5, 25, 50, 75,100, 125)
  ) +
  scale_color_manual(values = c(
    "Acute GVHD" = "firebrick",
    "Chronic GVHD" = "steelblue"
  )) +
  labs(
    x = "Time (months)",
    y = "Hazard ratio (Male vs Female)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

p_hr
#===========================================================
# Combined plot
#===========================================================

combined_plot <- p_haz / p_hr

print(combined_plot)

ggsave(
  filename = "Hazard_and_HR_application_sex.pdf",
  plot = combined_plot,
  width = 12,
  height = 10
)
############################################################
## Bivariate Gompertz-FGM model with LR test for H0: phi=0 ##
## Parallel version with CSV output per scenario and final ##
## Updated to vary both n and phi                            ##
############################################################

rm(list = ls(all = TRUE))
start.time <- Sys.time()

safe_log <- function(z, eps = 1e-12) log(pmax(z, eps))
clip01   <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

############################
## 1) Gompertz functions  ##
############################
Sg <- function(a, scale, t, eps = 1e-10) {
  ifelse(abs(a) < eps,
         exp(-scale * t),
         exp(-(scale / a) * (exp(a * t) - 1)))
}

fg <- function(a, scale, t, eps = 1e-10) {
  ifelse(abs(a) < eps,
         scale * exp(-scale * t),
         scale * exp(a * t) * exp(-(scale / a) * (exp(a * t) - 1)))
}

Fg <- function(a, scale, t, eps = 1e-10) 1 - Sg(a, scale, t, eps)

pcure <- function(a, scale, eps = 1e-10) {
  ifelse(a < -eps, exp(scale / a), 0)
}

############################
## 2) Joint survival      ##
############################
S12 <- function(alpha0, alpha1, beta0, beta1,
                gamma0, gamma1, lambda0, lambda1, phi,
                t1, t2, x, eps = 1e-10) {
  a1 <- alpha0 + alpha1 * x
  s1 <- exp(beta0 + beta1 * x)
  a2 <- gamma0 + gamma1 * x
  s2 <- exp(lambda0 + lambda1 * x)
  
  S1 <- Sg(a1, s1, t1, eps)
  S2 <- Sg(a2, s2, t2, eps)
  
  S1 * S2 * (1 + phi * (1 - S1) * (1 - S2))
}

dS1 <- function(alpha0, alpha1, beta0, beta1,
                gamma0, gamma1, lambda0, lambda1, phi,
                t1, t2, x, eps = 1e-10) {
  a1 <- alpha0 + alpha1 * x
  s1 <- exp(beta0 + beta1 * x)
  a2 <- gamma0 + gamma1 * x
  s2 <- exp(lambda0 + lambda1 * x)
  
  -fg(a1, s1, t1, eps) * Sg(a2, s2, t2, eps) +
    phi * Sg(a2, s2, t2, eps) * Fg(a2, s2, t2, eps) *
    fg(a1, s1, t1, eps) * (1 - 2 * Fg(a1, s1, t1, eps))
}

dS2 <- function(alpha0, alpha1, beta0, beta1,
                gamma0, gamma1, lambda0, lambda1, phi,
                t1, t2, x, eps = 1e-10) {
  a1 <- alpha0 + alpha1 * x
  s1 <- exp(beta0 + beta1 * x)
  a2 <- gamma0 + gamma1 * x
  s2 <- exp(lambda0 + lambda1 * x)
  
  -fg(a2, s2, t2, eps) * Sg(a1, s1, t1, eps) +
    phi * Sg(a1, s1, t1, eps) * Fg(a1, s1, t1, eps) *
    fg(a2, s2, t2, eps) * (1 - 2 * Fg(a2, s2, t2, eps))
}

d2S12 <- function(alpha0, alpha1, beta0, beta1,
                  gamma0, gamma1, lambda0, lambda1, phi,
                  t1, t2, x, eps = 1e-10) {
  a1 <- alpha0 + alpha1 * x
  s1 <- exp(beta0 + beta1 * x)
  a2 <- gamma0 + gamma1 * x
  s2 <- exp(lambda0 + lambda1 * x)
  
  fg(a1, s1, t1, eps) * fg(a2, s2, t2, eps) +
    phi * fg(a1, s1, t1, eps) * fg(a2, s2, t2, eps) *
    (1 - 2 * Fg(a1, s1, t1, eps)) * (1 - 2 * Fg(a2, s2, t2, eps))
}

############################
## 3) Negative loglik     ##
############################
nll_txa <- function(par, t1, t2, cens1, cens2, x, eps = 1e-12) {
  alpha0  <- par[1]
  alpha1  <- par[2]
  beta0   <- par[3]
  beta1   <- par[4]
  gamma0  <- par[5]
  gamma1  <- par[6]
  lambda0 <- par[7]
  lambda1 <- par[8]
  phi     <- tanh(par[9])
  
  ll <- sum(cens1 * cens2 *
              safe_log(d2S12(alpha0, alpha1, beta0, beta1,
                             gamma0, gamma1, lambda0, lambda1, phi,
                             t1, t2, x, eps), eps)) +
    sum(cens1 * (1 - cens2) *
          safe_log(-dS1(alpha0, alpha1, beta0, beta1,
                        gamma0, gamma1, lambda0, lambda1, phi,
                        t1, t2, x, eps), eps)) +
    sum(cens2 * (1 - cens1) *
          safe_log(-dS2(alpha0, alpha1, beta0, beta1,
                        gamma0, gamma1, lambda0, lambda1, phi,
                        t1, t2, x, eps), eps)) +
    sum((1 - cens1) * (1 - cens2) *
          safe_log(S12(alpha0, alpha1, beta0, beta1,
                       gamma0, gamma1, lambda0, lambda1, phi,
                       t1, t2, x, eps), eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}

nll_txa_phi0 <- function(par, t1, t2, cens1, cens2, x, eps = 1e-12) {
  alpha0  <- par[1]
  alpha1  <- par[2]
  beta0   <- par[3]
  beta1   <- par[4]
  gamma0  <- par[5]
  gamma1  <- par[6]
  lambda0 <- par[7]
  lambda1 <- par[8]
  phi     <- 0
  
  ll <- sum(cens1 * cens2 *
              safe_log(d2S12(alpha0, alpha1, beta0, beta1,
                             gamma0, gamma1, lambda0, lambda1, phi,
                             t1, t2, x, eps), eps)) +
    sum(cens1 * (1 - cens2) *
          safe_log(-dS1(alpha0, alpha1, beta0, beta1,
                        gamma0, gamma1, lambda0, lambda1, phi,
                        t1, t2, x, eps), eps)) +
    sum(cens2 * (1 - cens1) *
          safe_log(-dS2(alpha0, alpha1, beta0, beta1,
                        gamma0, gamma1, lambda0, lambda1, phi,
                        t1, t2, x, eps), eps)) +
    sum((1 - cens1) * (1 - cens2) *
          safe_log(S12(alpha0, alpha1, beta0, beta1,
                       gamma0, gamma1, lambda0, lambda1, phi,
                       t1, t2, x, eps), eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}

############################
## 4) Simulation          ##
############################
qSg <- function(u, a, scale, eps = 1e-10) {
  u <- clip01(u)
  if (abs(a) < eps) return(-log(u) / scale)
  
  p_inf <- pcure(a, scale, eps)
  if (a < -eps && u < p_inf) return(Inf)
  
  val <- 1 - (a / scale) * log(u)
  if (!is.finite(val) || val <= 0) return(Inf)
  
  t <- log(val) / a
  if (!is.finite(t) || t < 0) t <- Inf
  t
}

rfgm <- function(n, phi, eps = 1e-10) {
  u1 <- runif(n)
  w  <- runif(n)
  u2 <- numeric(n)
  
  for (i in seq_len(n)) {
    theta <- phi * (1 - 2 * u1[i])
    if (abs(theta) < eps) {
      u2[i] <- w[i]
    } else {
      disc <- max((1 + theta)^2 - 4 * theta * w[i], 0)
      r1 <- ((1 + theta) - sqrt(disc)) / (2 * theta)
      r2 <- ((1 + theta) + sqrt(disc)) / (2 * theta)
      cand <- c(r1, r2)
      cand <- cand[cand >= 0 & cand <= 1]
      u2[i] <- if (length(cand)) cand[1] else w[i]
    }
  }
  
  cbind(u1 = u1, u2 = u2)
}

sim_data <- function(n,
                     alpha0, alpha1, beta0, beta1,
                     gamma0, gamma1, lambda0, lambda1, phi,
                     x = NULL, cmax1 = 25, cmax2 = 25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(x)) x <- rbinom(n, 1, 0.406)
  if (length(x) == 1) x <- rep(x, n)
  
  uv <- rfgm(n, phi)
  
  a1 <- alpha0 + alpha1 * x
  s1 <- exp(beta0 + beta1 * x)
  a2 <- gamma0 + gamma1 * x
  s2 <- exp(lambda0 + lambda1 * x)
  
  t1_true <- mapply(function(u, a, sc) qSg(u, a, sc), uv[, 1], a1, s1)
  t2_true <- mapply(function(u, a, sc) qSg(u, a, sc), uv[, 2], a2, s2)
  
  c1_obs <- runif(n, 0, cmax1)
  c2_obs <- runif(n, 0, cmax2)
  
  t1 <- ifelse(is.infinite(t1_true), c1_obs, pmin(t1_true, c1_obs))
  t2 <- ifelse(is.infinite(t2_true), c2_obs, pmin(t2_true, c2_obs))
  
  cens1 <- ifelse(is.infinite(t1_true), 0L, as.integer(t1_true <= c1_obs))
  cens2 <- ifelse(is.infinite(t2_true), 0L, as.integer(t2_true <= c2_obs))
  cure1 <- as.integer(is.infinite(t1_true))
  cure2 <- as.integer(is.infinite(t2_true))
  
  data.frame(
    x = x,
    t1 = t1,
    t2 = t2,
    cens1 = cens1,
    cens2 = cens2,
    cure1 = cure1,
    cure2 = cure2
  )
}

############################
## 5) Estimation          ##
############################
fit_model_free <- function(
    dat,
    init = c(alpha0 = -0.70, alpha1 = -0.1,
             beta0 = log(0.3), beta1 = 0.2,
             gamma0 = -0.1, gamma1 = 0.05,
             lambda0 = log(0.25), lambda1 = 0.05,
             phi_z = atanh(0.10))
) {
  fit <- try(
    optim(
      par = init,
      fn = nll_txa,
      t1 = dat$t1, t2 = dat$t2,
      cens1 = dat$cens1, cens2 = dat$cens2,
      x = dat$x,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = 2000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      converged = FALSE,
      logLik = NA_real_,
      fit = NULL,
      est = setNames(rep(NA_real_, 9),
                     c("alpha0", "alpha1", "beta0", "beta1",
                       "gamma0", "gamma1", "lambda0", "lambda1", "phi"))
    ))
  }
  
  est_nat <- c(fit$par[1:8], tanh(fit$par[9]))
  names(est_nat) <- c("alpha0", "alpha1", "beta0", "beta1",
                      "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  
  list(
    converged = fit$convergence == 0,
    logLik = -fit$value,
    fit = fit,
    est = est_nat
  )
}

fit_model_phi0 <- function(
    dat,
    init = c(alpha0 = -0.70, alpha1 = -0.1,
             beta0 = log(0.3), beta1 = 0.2,
             gamma0 = -0.1, gamma1 = 0.05,
             lambda0 = log(0.25), lambda1 = 0.05)
) {
  fit <- try(
    optim(
      par = init,
      fn = nll_txa_phi0,
      t1 = dat$t1, t2 = dat$t2,
      cens1 = dat$cens1, cens2 = dat$cens2,
      x = dat$x,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = 2000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      converged = FALSE,
      logLik = NA_real_,
      fit = NULL,
      est = setNames(rep(NA_real_, 9),
                     c("alpha0", "alpha1", "beta0", "beta1",
                       "gamma0", "gamma1", "lambda0", "lambda1", "phi"))
    ))
  }
  
  est_nat <- c(fit$par[1:8], 0)
  names(est_nat) <- c("alpha0", "alpha1", "beta0", "beta1",
                      "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  
  list(
    converged = fit$convergence == 0,
    logLik = -fit$value,
    fit = fit,
    est = est_nat
  )
}

############################
## 6) Likelihood ratio    ##
############################
lr_test_phi <- function(fit_free, fit_phi0) {
  if (is.null(fit_free) || is.null(fit_phi0) ||
      !isTRUE(fit_free$converged) || !isTRUE(fit_phi0$converged) ||
      !is.finite(fit_free$logLik) || !is.finite(fit_phi0$logLik)) {
    return(c(LR = NA_real_, p_value = NA_real_))
  }
  
  LR <- 2 * (fit_free$logLik - fit_phi0$logLik)
  LR <- pmax(LR, 0)
  p_value <- 1 - pchisq(LR, df = 1)
  c(LR = LR, p_value = p_value)
}

############################
## 7) One replication     ##
############################
one_rep_lrt <- function(n, true, cmax1 = 132.4, cmax2 = 132.4,
                        seed = NULL, alpha_test = 0.05) {
  dat <- sim_data(
    n = n,
    alpha0 = true["alpha0"], alpha1 = true["alpha1"],
    beta0 = true["beta0"], beta1 = true["beta1"],
    gamma0 = true["gamma0"], gamma1 = true["gamma1"],
    lambda0 = true["lambda0"], lambda1 = true["lambda1"],
    phi = true["phi"],
    cmax1 = cmax1, cmax2 = cmax2,
    seed = seed
  )
  
  fit_free <- fit_model_free(dat)
  init_phi0 <- if (isTRUE(fit_free$converged) && !is.null(fit_free$fit)) {
    fit_free$fit$par[1:8]
  } else {
    NULL
  }
  
  fit_phi0 <- fit_model_phi0(
    dat,
    init = if (is.null(init_phi0)) {
      c(alpha0 = -0.70, alpha1 = -0.1,
        beta0 = log(0.3), beta1 = 0.2,
        gamma0 = -0.1, gamma1 = 0.05,
        lambda0 = log(0.25), lambda1 = 0.05)
    } else {
      init_phi0
    }
  )
  
  lrt <- lr_test_phi(fit_free, fit_phi0)
  
  data.frame(
    converged_free = as.integer(isTRUE(fit_free$converged)),
    converged_phi0 = as.integer(isTRUE(fit_phi0$converged)),
    phi_true = unname(true["phi"]),
    phi_hat = unname(fit_free$est["phi"]),
    logLik_free = fit_free$logLik,
    logLik_phi0 = fit_phi0$logLik,
    LR = unname(lrt["LR"]),
    p_value = unname(lrt["p_value"]),
    reject_5pct = as.integer(!is.na(lrt["p_value"]) && lrt["p_value"] < alpha_test)
    #cure_rate_T1 = mean(dat$cure1),
    #cure_rate_T2 = mean(dat$cure2),
    #event_rate_T1 = mean(dat$cens1),
   # event_rate_T2 = mean(dat$cens2)
  )
}

############################
## 8) Summary functions   ##
############################
summarize_lrt_scenario <- function(mc_lrt, scenario_label, n_value, phi_value,
                                   alpha_test = 0.05) {
  data.frame(
    scenario = scenario_label,
    n = n_value,
    phi_true = phi_value,
    replications = nrow(mc_lrt),
    valid_p = sum(!is.na(mc_lrt$p_value)),
    convergence_rate_free = mean(mc_lrt$converged_free, na.rm = TRUE),
    convergence_rate_phi0 = mean(mc_lrt$converged_phi0, na.rm = TRUE),
    mean_phi_hat = mean(mc_lrt$phi_hat, na.rm = TRUE),
    sd_phi_hat = sd(mc_lrt$phi_hat, na.rm = TRUE),
  #  mean_LR = mean(mc_lrt$LR, na.rm = TRUE),
  #  mean_p_value = mean(mc_lrt$p_value, na.rm = TRUE),
  #  median_p_value = median(mc_lrt$p_value, na.rm = TRUE),
    rejection_rate = mean(mc_lrt$reject_5pct, na.rm = TRUE),
   # mean_event_rate_T1 = mean(mc_lrt$event_rate_T1, na.rm = TRUE),
   #mean_event_rate_T2 = mean(mc_lrt$event_rate_T2, na.rm = TRUE),
   # mean_cure_rate_T1 = mean(mc_lrt$cure_rate_T1, na.rm = TRUE),
   # mean_cure_rate_T2 = mean(mc_lrt$cure_rate_T2, na.rm = TRUE),
    alpha_test = alpha_test
  )
}

make_scenario_label <- function(n_value, phi_value, digits = 3) {
  phi_txt <- sprintf(paste0("%.", digits, "f"), phi_value)
  phi_txt <- gsub("-", "m", phi_txt)
  phi_txt <- gsub("\\.", "_", phi_txt)
  paste0("n", n_value, "_phi_", phi_txt)
}

############################
## 9) Parallel engine     ##
############################
run_one_phi_parallel <- function(R, n, phi_value, true_base,
                                 cmax1 = 132.4, cmax2 = 132.4,
                                 seed = 123, alpha_test = 0.05,
                                 n_cores = max(1, parallel::detectCores() - 1),
                                 type = "PSOCK") {
  true_j <- true_base
  true_j["phi"] <- phi_value
  
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, R)
  
  n_cores <- max(1, min(n_cores, R))
  cl <- parallel::makeCluster(n_cores, type = type)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "safe_log", "clip01",
      "Sg", "fg", "Fg", "pcure",
      "S12", "dS1", "dS2", "d2S12",
      "nll_txa", "nll_txa_phi0",
      "qSg", "rfgm", "sim_data",
      "fit_model_free", "fit_model_phi0",
      "lr_test_phi", "one_rep_lrt"
    ),
    envir = environment()
  )
  
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  res_j <- parallel::parLapply(
    cl,
    X = seq_len(R),
    fun = function(r, n, true_j, cmax1, cmax2, seeds, alpha_test) {
      one_rep_lrt(
        n = n,
        true = true_j,
        cmax1 = cmax1,
        cmax2 = cmax2,
        seed = seeds[r],
        alpha_test = alpha_test
      )
    },
    n = n,
    true_j = true_j,
    cmax1 = cmax1,
    cmax2 = cmax2,
    seeds = seeds,
    alpha_test = alpha_test
  )
  
  out <- do.call(rbind, res_j)
  rownames(out) <- NULL
  out
}

run_lrt_scenarios_parallel <- function(R,
                                       n_values,
                                       phi_values,
                                       true_base,
                                       cmax1 = 132.4,
                                       cmax2 = 132.4,
                                       seed = 123,
                                       alpha_test = 0.05,
                                       n_cores = max(1, parallel::detectCores() - 1),
                                       type = if (.Platform$OS.type == "windows") "PSOCK" else "PSOCK",
                                       output_dir = "simulation_results_lrt_parallel") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  scenarios <- expand.grid(
    n = n_values,
    phi = phi_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  scenario_results <- vector("list", nrow(scenarios))
  scenario_summaries <- vector("list", nrow(scenarios))
  
  for (j in seq_len(nrow(scenarios))) {
    n_j <- scenarios$n[j]
    phi_j <- scenarios$phi[j]
    scenario_label <- make_scenario_label(n_j, phi_j)
    
    cat("\n------------------------------------------\n")
    cat("Running scenario:", scenario_label, "\n")
    cat("------------------------------------------\n")
    
    res_j <- run_one_phi_parallel(
      R = R,
      n = n_j,
      phi_value = phi_j,
      true_base = true_base,
      cmax1 = cmax1,
      cmax2 = cmax2,
      seed = seed + j,
      alpha_test = alpha_test,
      n_cores = n_cores,
      type = type
    )
    
    res_j$scenario <- scenario_label
    res_j$n <- n_j
    res_j$alpha_test <- alpha_test
    
    summary_j <- summarize_lrt_scenario(
      mc_lrt = res_j,
      scenario_label = scenario_label,
      n_value = n_j,
      phi_value = phi_j,
      alpha_test = alpha_test
    )
    
    scenario_results[[j]] <- res_j
    scenario_summaries[[j]] <- summary_j
    
    utils::write.csv(
      res_j,
      file = file.path(output_dir, paste0("results_", scenario_label, ".csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      summary_j,
      file = file.path(output_dir, paste0("summary_", scenario_label, ".csv")),
      row.names = FALSE
    )
  }
  
  final_results <- do.call(rbind, scenario_results)
  final_summary <- do.call(rbind, scenario_summaries)
  rownames(final_results) <- NULL
  rownames(final_summary) <- NULL
  
  final_summary <- final_summary[order(final_summary$n, final_summary$phi_true), ]
  
  utils::write.csv(
    final_results,
    file = file.path(output_dir, "results_all_scenarios.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    final_summary,
    file = file.path(output_dir, "summary_all_scenarios.csv"),
    row.names = FALSE
  )
  
  end.time <- Sys.time()
  elapsed <- difftime(end.time, start.time, units = "mins")
  cat("Total elapsed time:", round(as.numeric(elapsed), 2), "minutes\n")
  
  list(
    results_all = final_results,
    summary_all = final_summary
  )
}

############################
## 10) User parameters    ##
############################
true_par <- c(
  alpha0 = -0.7077,
  alpha1 = -0.1669,
  beta0 = -1.2404,
  beta1 = 0.4399,
  gamma0 = -0.0560,
  gamma1 = 0.0025,
  lambda0 = -2.3630,
  lambda1 = -0.0504,
  phi = -0.8
)

############################
## 11) Example            ##
############################
## Run this block in R:
##
out_lrt <- run_lrt_scenarios_parallel(
  R = 1000,
  n_values = c(100, 200, 400, 600, 800, 1000),
  phi_values = c(-0.8, -0.6, -0.4,-0.2, 0, 0.2, 0.4, 0.6,0.8),
  true_base = true_par,
  cmax1 = 132.4,
  cmax2 = 132.4,
  seed = 20260403,
  alpha_test = 0.05,
  n_cores = 20,
  output_dir = "simulation_results_lrt"
)

out_lrt$summary_all


library(ggplot2)
library(scales)

plot_data <- out_lrt$summary_all
plot_data$n <- factor(plot_data$n)

p <- ggplot(
  plot_data,
  aes(
    x = phi_true,
    y = rejection_rate,
    group = n,
    color = n,
    shape = n,
    linetype = n
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = sort(unique(plot_data$phi_true))) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  ) +
  scale_shape_manual(values = c(16, 17, 15, 18, 22,21)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
  labs(
    x = expression(phi),
    y = "Rejection rate (%)",
    color = "Sample size (n)",
    shape = "Sample size (n)",
    linetype = "Sample size (n)",
    title = ""
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top"
  )

p
ggsave(filename = "plot_simulation_results_lrt.pdf",plot = p,width = 9,height = 7,dpi = 300)

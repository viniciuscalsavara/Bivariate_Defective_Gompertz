############################################################
## Bivariate Gompertz-FGM model with multiple covariates  ##
## Complete single-file script with Monte Carlo summary    ##
## Incremental saving in .csv by scenario and by block     ##
## Example covariates: z1 ~ Bernoulli(0.5), z2 ~ N(0,1)    ##
############################################################
rm(list = ls())
start.time <- Sys.time()

############################
## 1) Utilities           ##
############################
library(truncnorm)
safe_log <- function(z, eps = 1e-12) log(pmax(z, eps))
clip01   <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

append_csv_safely <- function(df, file) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    df,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file),
    quote = TRUE,
    na = "NA"
  )
}

############################
## 2) Gompertz functions  ##
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

############################
## 3) FGM copula          ##
############################

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

############################
## 4) Design matrix       ##
############################

build_X <- function(z1, z2) {
  cbind(1, z1, z2)
}

make_par_names <- function(p) {
  c(
    paste0("alpha", 0:(p - 1)),
    paste0("b",     0:(p - 1)),
    paste0("gamma", 0:(p - 1)),
    paste0("d",     0:(p - 1)),
    "phi"
  )
}

split_par <- function(par, p) {
  list(
    alpha = par[1:p],
    beta  = par[(p + 1):(2 * p)],
    gamma = par[(2 * p + 1):(3 * p)],
    delta = par[(3 * p + 1):(4 * p)],
    phi_z = par[4 * p + 1]
  )
}

linpred_components <- function(X, alpha, beta, gamma, delta) {
  a1   <- as.vector(X %*% alpha)
  eta1 <- as.vector(X %*% beta)
  a2   <- as.vector(X %*% gamma)
  eta2 <- as.vector(X %*% delta)

  list(
    a1 = a1,
    s1 = exp(eta1),
    a2 = a2,
    s2 = exp(eta2)
  )
}

############################
## 5) Joint survival      ##
############################

S12_general <- function(par, t1, t2, X, eps = 1e-10) {
  p <- ncol(X)
  pa <- split_par(par, p)
  co <- linpred_components(X, pa$alpha, pa$beta, pa$gamma, pa$delta)
  phi <- tanh(pa$phi_z)

  S1 <- Sg(co$a1, co$s1, t1, eps)
  S2 <- Sg(co$a2, co$s2, t2, eps)

  S1 * S2 * (1 + phi * (1 - S1) * (1 - S2))
}

dS1_general <- function(par, t1, t2, X, eps = 1e-10) {
  p <- ncol(X)
  pa <- split_par(par, p)
  co <- linpred_components(X, pa$alpha, pa$beta, pa$gamma, pa$delta)
  phi <- tanh(pa$phi_z)

  -fg(co$a1, co$s1, t1, eps) * Sg(co$a2, co$s2, t2, eps) +
    phi * Sg(co$a2, co$s2, t2, eps) * Fg(co$a2, co$s2, t2, eps) *
    fg(co$a1, co$s1, t1, eps) * (1 - 2 * Fg(co$a1, co$s1, t1, eps))
}

dS2_general <- function(par, t1, t2, X, eps = 1e-10) {
  p <- ncol(X)
  pa <- split_par(par, p)
  co <- linpred_components(X, pa$alpha, pa$beta, pa$gamma, pa$delta)
  phi <- tanh(pa$phi_z)

  -fg(co$a2, co$s2, t2, eps) * Sg(co$a1, co$s1, t1, eps) +
    phi * Sg(co$a1, co$s1, t1, eps) * Fg(co$a1, co$s1, t1, eps) *
    fg(co$a2, co$s2, t2, eps) * (1 - 2 * Fg(co$a2, co$s2, t2, eps))
}

d2S12_general <- function(par, t1, t2, X, eps = 1e-10) {
  p <- ncol(X)
  pa <- split_par(par, p)
  co <- linpred_components(X, pa$alpha, pa$beta, pa$gamma, pa$delta)
  phi <- tanh(pa$phi_z)

  fg(co$a1, co$s1, t1, eps) * fg(co$a2, co$s2, t2, eps) +
    phi * fg(co$a1, co$s1, t1, eps) * fg(co$a2, co$s2, t2, eps) *
    (1 - 2 * Fg(co$a1, co$s1, t1, eps)) * (1 - 2 * Fg(co$a2, co$s2, t2, eps))
}

############################
## 6) Negative loglik     ##
############################

nll_general <- function(par, t1, t2, cens1, cens2, X, eps = 1e-12) {
  ll <- sum(cens1 * cens2 * safe_log(d2S12_general(par, t1, t2, X, eps), eps)) +
    sum(cens1 * (1 - cens2) * safe_log(-dS1_general(par, t1, t2, X, eps), eps)) +
    sum(cens2 * (1 - cens1) * safe_log(-dS2_general(par, t1, t2, X, eps), eps)) +
    sum((1 - cens1) * (1 - cens2) * safe_log(S12_general(par, t1, t2, X, eps), eps))

  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}

############################
## 7) Data generation     ##
############################

sim_data_general <- function(n,
                             alpha, beta, gamma, delta, phi,
                             z1 = NULL, z2 = NULL,
                             cmax1 = 132, cmax2 = 132,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(z1)) z1 <- rbinom(n, 1, 0.406)
  #if (is.null(z2)) z2 <- rbinom(n, 1, 0.5)
  #if (is.null(z2)) z2 <- rtruncnorm(n, a=5, b=53, mean=32, sd=11)
  if (is.null(z2)) z2 <- rnorm(n, mean=0, sd=1)

  X <- build_X(z1, z2)
  uv <- rfgm(n, phi)
  co <- linpred_components(X, alpha, beta, gamma, delta)

  t1_true <- mapply(function(u, a, sc) qSg(u, a, sc), uv[, 1], co$a1, co$s1)
  t2_true <- mapply(function(u, a, sc) qSg(u, a, sc), uv[, 2], co$a2, co$s2)

  c1_obs <- runif(n, 0, cmax1)
  c2_obs <- runif(n, 0, cmax2)

  t1 <- ifelse(is.infinite(t1_true), c1_obs, pmin(t1_true, c1_obs))
  t2 <- ifelse(is.infinite(t2_true), c2_obs, pmin(t2_true, c2_obs))

  cens1 <- ifelse(is.infinite(t1_true), 0L, as.integer(t1_true <= c1_obs))
  cens2 <- ifelse(is.infinite(t2_true), 0L, as.integer(t2_true <= c2_obs))

  data.frame(
    z1 = z1,
    z2 = z2,
    t1 = t1,
    t2 = t2,
    cens1 = cens1,
    cens2 = cens2
  )
}

############################
## 8) Estimation          ##
############################

fit_model_general <- function(dat,
                              init_alpha = c(-1.4, 0.2, 0),
                              init_beta  = c(0.7, 0.7, -0.1),
                              init_gamma = c(-0.5, 0.1, 0),
                              init_delta = c(-1, -0.6, 0),
                              init_phi   = atanh(0.10)) {
  X <- build_X(dat$z1, dat$z2)
  p <- ncol(X)

  init <- c(init_alpha, init_beta, init_gamma, init_delta, init_phi)
  par_names <- make_par_names(p)

  fit <- try(
    optim(
      par = init,
      fn = nll_general,
      t1 = dat$t1,
      t2 = dat$t2,
      cens1 = dat$cens1,
      cens2 = dat$cens2,
      X = X,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = 3000, reltol = 1e-10)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    est <- setNames(rep(NA_real_, length(par_names)), par_names)
    se  <- setNames(rep(NA_real_, length(par_names)), paste0("SE_", par_names))
    return(list(converged = FALSE, est = est, se = se, vc = NULL, logLik = NA_real_, fit = NULL))
  }

  est_int <- fit$par
  est_nat <- est_int
  est_nat[length(est_nat)] <- tanh(est_int[length(est_int)])
  names(est_nat) <- par_names

  vc_int <- try(solve(fit$hessian), silent = TRUE)

  if (inherits(vc_int, "try-error") || any(!is.finite(vc_int))) {
    vc_nat <- matrix(NA_real_, length(par_names), length(par_names),
                     dimnames = list(par_names, par_names))
    se_nat <- setNames(rep(NA_real_, length(par_names)), paste0("SE_", par_names))
  } else {
    J <- diag(length(par_names))
    J[length(par_names), length(par_names)] <- 1 - tanh(est_int[length(est_int)])^2
    vc_nat <- J %*% vc_int %*% t(J)
    colnames(vc_nat) <- rownames(vc_nat) <- par_names
    se_nat <- setNames(sqrt(diag(vc_nat)), paste0("SE_", par_names))
  }

  list(
    converged = isTRUE(fit$convergence == 0),
    est = est_nat,
    se = se_nat,
    vc = vc_nat,
    logLik = -fit$value,
    fit = fit
  )
}

############################
## 9) One MC replication  ##
############################

one_rep_general <- function(n, true, cmax1 = 132, cmax2 = 132, seed = NULL, rep_id = NA_integer_) {
  dat <- sim_data_general(
    n = n,
    alpha = true$alpha,
    beta  = true$beta,
    gamma = true$gamma,
    delta = true$delta,
    phi   = true$phi,
    cmax1 = cmax1,
    cmax2 = cmax2,
    seed = seed
  )

  fit <- fit_model_general(dat)

  out <- data.frame(
    rep = rep_id,
    converged = as.integer(fit$converged),
    event_rate_T1 = mean(dat$cens1),
    event_rate_T2 = mean(dat$cens2),
    check.names = FALSE
  )

  for (nm in names(fit$est)) out[[nm]] <- unname(fit$est[[nm]])
  for (nm in names(fit$se))  out[[nm]] <- unname(fit$se[[nm]])

  out <- out[, c("rep", "converged",
                 names(fit$est), names(fit$se),
                 "event_rate_T1", "event_rate_T2"), drop = FALSE]
  rownames(out) <- NULL
  out
}

############################
## 10) Parallel MC study  ##
############################

mc_study_parallel_general <- function(R, n, true,
                                      cmax1 = 132, cmax2 = 132,
                                      seed = 123,
                                      n_cores = max(1, parallel::detectCores() - 1),
                                      type = if (.Platform$OS.type == "windows") "PSOCK" else "PSOCK") {
  if (n_cores < 1) n_cores <- 1
  n_cores <- min(n_cores, R)

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, R)

  cl <- parallel::makeCluster(n_cores, type = type)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(truncnorm)
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "safe_log", "clip01", "Sg", "fg", "Fg", "pcure", "qSg", "rfgm",
      "build_X", "make_par_names", "split_par", "linpred_components",
      "S12_general", "dS1_general", "dS2_general", "d2S12_general", "nll_general",
      "sim_data_general", "fit_model_general", "one_rep_general"
    ),
    envir = environment()
  )

  res <- parallel::parLapply(
    cl,
    X = seq_len(R),
    fun = function(r, n, true, cmax1, cmax2, seeds) {
      one_rep_general(
        n = n, true = true,
        cmax1 = cmax1, cmax2 = cmax2,
        seed = seeds[r], rep_id = r
      )
    },
    n = n, true = true, cmax1 = cmax1, cmax2 = cmax2, seeds = seeds
  )

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

############################
## 11) Sequential blocks  ##
############################

mc_study_save_by_blocks_general <- function(R, n, true,
                                            cmax1 = 132, cmax2 = 132,
                                            seed = 123,
                                            block_size = 50,
                                            scenario_id = "scenario",
                                            output_dir = ".",
                                            save_csv = TRUE) {
  if (block_size < 1) block_size <- 1
  block_starts <- seq(1, R, by = block_size)
  mc_file <- file.path(output_dir, paste0("mc_results_", scenario_id, ".csv"))

  if (isTRUE(save_csv) && file.exists(mc_file)) file.remove(mc_file)

  all_blocks <- vector("list", length(block_starts))

  for (b in seq_along(block_starts)) {
    start_rep <- block_starts[b]
    end_rep <- min(start_rep + block_size - 1, R)
    reps_here <- start_rep:end_rep

    cat(sprintf("  Block %d/%d: reps %d-%d\n", b, length(block_starts), start_rep, end_rep))

    block_res <- mc_study_parallel_general(
      R = length(reps_here),
      n = n,
      true = true,
      cmax1 = cmax1,
      cmax2 = cmax2,
      seed = seed + b - 1,
      n_cores = min(max(1, parallel::detectCores() - 1), length(reps_here))
    )

    block_res$rep <- reps_here
    block_res$Scenario <- scenario_id
    block_res$n <- n
    block_res$phi_true <- true$phi

    block_res <- block_res[, c("Scenario", "n", "phi_true", "rep",
                               setdiff(names(block_res), c("Scenario", "n", "phi_true", "rep"))),
                           drop = FALSE]

    if (isTRUE(save_csv)) append_csv_safely(block_res, mc_file)
    all_blocks[[b]] <- block_res
  }

  out <- do.call(rbind, all_blocks)
  rownames(out) <- NULL
  out
}

############################
## 12) MC summary         ##
############################

true_vector_general <- function(true) {
  tv <- c(true$alpha, true$beta, true$gamma, true$delta, true$phi)
  names(tv) <- make_par_names(length(true$alpha))
  tv
}

summarize_mc_general <- function(mc_res, true,
                                 scenario_id = NA_character_,
                                 n = NA_integer_,
                                 phi_true = NA_real_) {
  true_vec <- true_vector_general(true)
  pars <- names(true_vec)
  conv <- mc_res$converged == 1

  tab_par <- do.call(rbind, lapply(pars, function(p) {
    est <- mc_res[conv, p]
    se  <- mc_res[conv, paste0("SE_", p)]
    tru <- unname(true_vec[p])

    data.frame(
      Scenario = scenario_id,
      n = n,
      phi_true = phi_true,
      Parameter = p,
      Truth = tru,
      Mean_MLE = if (all(is.na(est))) NA_real_ else mean(est, na.rm = TRUE),
      Bias = if (all(is.na(est))) NA_real_ else mean(est - tru, na.rm = TRUE),
      SD_MLE = if (all(is.na(est))) NA_real_ else stats::sd(est, na.rm = TRUE),
      Mean_SE = if (all(is.na(se))) NA_real_ else mean(se, na.rm = TRUE),
      RMSE = if (all(is.na(est))) NA_real_ else sqrt(mean((est - tru)^2, na.rm = TRUE)),
      Coverage_95 = if (all(is.na(est)) || all(is.na(se))) {
        NA_real_
      } else {
        mean(est - 1.96 * se <= tru & tru <= est + 1.96 * se, na.rm = TRUE)
      },
      Convergence_Rate = mean(mc_res$converged == 1, na.rm = TRUE),
      Mean_Event_Rate_T1 = mean(mc_res$event_rate_T1, na.rm = TRUE),
      Mean_Event_Rate_T2 = mean(mc_res$event_rate_T2, na.rm = TRUE),
      N_Converged = sum(conv, na.rm = TRUE),
      R = nrow(mc_res)
    )
  }))

  rownames(tab_par) <- NULL
  tab_par
}

############################
## 13) Multiple scenarios ##
############################

run_mc_scenarios_general <- function(scenarios,
                                     n_values,
                                     R = 500,
                                     cmax1 = 132,
                                     cmax2 = 132,
                                     seed = 123,
                                     n_cores = max(1, parallel::detectCores() - 1),
                                     output_dir = NULL,
                                     save_csv = TRUE,
                                     save_incremental = TRUE,
                                     block_size = 50) {
  all_mc <- list()
  all_sum <- list()
  counter <- 1L

  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (s in seq_along(scenarios)) {
    true <- scenarios[[s]]
    scenario_name <- names(scenarios)[s]
    if (is.null(scenario_name) || identical(scenario_name, "")) {
      scenario_name <- paste0("scenario_", s)
    }

    for (n_i in n_values) {
      scenario_id <- paste0(scenario_name, "_n", n_i)
      cat("Running scenario:", scenario_id, "\n")

      if (isTRUE(save_incremental) && !is.null(output_dir)) {
        mc_res <- mc_study_save_by_blocks_general(
          R = R,
          n = n_i,
          true = true,
          cmax1 = cmax1,
          cmax2 = cmax2,
          seed = seed + counter,
          block_size = block_size,
          scenario_id = scenario_id,
          output_dir = output_dir,
          save_csv = save_csv
        )
      } else {
        mc_res <- mc_study_parallel_general(
          R = R,
          n = n_i,
          true = true,
          cmax1 = cmax1,
          cmax2 = cmax2,
          seed = seed + counter,
          n_cores = n_cores
        )
        mc_res$Scenario <- scenario_id
        mc_res$n <- n_i
        mc_res$phi_true <- true$phi
        mc_res <- mc_res[, c("Scenario", "n", "phi_true", "rep",
                             setdiff(names(mc_res), c("Scenario", "n", "phi_true", "rep"))),
                         drop = FALSE]

        if (!is.null(output_dir) && isTRUE(save_csv)) {
          utils::write.csv(mc_res,
                           file.path(output_dir, paste0("mc_results_", scenario_id, ".csv")),
                           row.names = FALSE)
        }
      }

      sum_res <- summarize_mc_general(
        mc_res = mc_res,
        true = true,
        scenario_id = scenario_id,
        n = n_i,
        phi_true = true$phi
      )

      all_mc[[counter]] <- mc_res
      all_sum[[counter]] <- sum_res

      if (!is.null(output_dir) && isTRUE(save_csv)) {
        utils::write.csv(sum_res,
                         file.path(output_dir, paste0("summary_", scenario_id, ".csv")),
                         row.names = FALSE)
      }

      counter <- counter + 1L
    }
  }

  all_mc_df <- do.call(rbind, all_mc)
  all_sum_df <- do.call(rbind, all_sum)
  rownames(all_mc_df) <- NULL
  rownames(all_sum_df) <- NULL

  if (!is.null(output_dir) && isTRUE(save_csv)) {
    utils::write.csv(all_mc_df,
                     file.path(output_dir, "all_mc_results_general.csv"),
                     row.names = FALSE)
    utils::write.csv(all_sum_df,
                     file.path(output_dir, "overall_summary_general.csv"),
                     row.names = FALSE)
  }

  list(all_mc_results = all_mc_df, overall_summary = all_sum_df)
}

############################
## 14) Optional ggplot    ##
############################

plot_all_parameters_by_n_gg <- function(summary_tab,
                                        parameters = unique(summary_tab$Parameter),
                                        phi_values = NULL,
                                        metrics = c("Bias", "RMSE", "Coverage_95")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggh4x", quietly = TRUE)) stop("Package 'ggh4x' is required.")

  needed_cols <- c("n", "phi_true", "Parameter", metrics)
  miss <- setdiff(needed_cols, names(summary_tab))
  if (length(miss) > 0) stop("summary_tab is missing: ", paste(miss, collapse = ", "))

  dat <- summary_tab[summary_tab$Parameter %in% parameters, , drop = FALSE]
  if (!is.null(phi_values)) dat <- dat[dat$phi_true %in% phi_values, , drop = FALSE]

  long_dat <- do.call(
    rbind,
    lapply(metrics, function(m) {
      data.frame(
        n = dat$n,
        phi_true = dat$phi_true,
        Parameter = dat$Parameter,
        Metric = m,
        Value = dat[[m]]
      )
    })
  )

  long_dat$Parameter <- factor(long_dat$Parameter, levels = parameters)
  long_dat$Metric <- factor(long_dat$Metric,
                            levels = c("Bias", "RMSE", "Coverage_95"),
                            labels = c("Bias", "RMSE", "CP95%"))
  long_dat$phi_true <- factor(
    long_dat$phi_true,
    levels = sort(unique(long_dat$phi_true)),
    labels = paste0("phi = ", sort(unique(long_dat$phi_true)))
  )

  ggplot2::ggplot(long_dat, ggplot2::aes(x = n, y = Value, color = phi_true, group = phi_true)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggh4x::facet_grid2(Parameter ~ Metric, scales = "free_y", independent = "y") +
    ggplot2::labs(x = "Sample size (n)", y = NULL, color = "Scenario") +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )    +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("CP95%", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0.93649),
      inherit.aes = FALSE,
      linetype = 2
    ) +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("CP95%", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0.9635),
      inherit.aes = FALSE,
      linetype = 2
    ) +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("Bias", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0),
      inherit.aes = FALSE,
      linetype = 2
    )
}

############################
## 15) Example settings   ##
############################

#age in the raw data
# true_example <- list(
#   alpha = c( -1.3752,0.2217, 0.0163),
#   beta  = c(0.754,0.7015,-0.0719),
#   gamma = c(-0.4063,0.0872,0.0094),
#   delta = c(-0.9672,-0.5959,-0.0322),
#   phi   = -0.8
# )

#coefficients after standardized age

true_example <- list(
  alpha = c(-0.8536, 0.2217, 0.1793),
  beta  = c(-1.5468, 0.7015, -0.7909),
  gamma = c(-0.1055, 0.0872, 0.1034),
  delta = c(-1.9976, -0.5959, -0.3542),   # lambdas
  phi   = -0.8
)

scenarios_example <- list(
  scenario_A = true_example,
  scenario_B = modifyList(true_example, list(phi =  -0.7255)),
  scenario_C = modifyList(true_example, list(phi = -0.6)),
  scenario_D = modifyList(true_example, list(phi = -0.4)),
  scenario_E = modifyList(true_example, list(phi = -0.2)),
  scenario_F = modifyList(true_example, list(phi = 0)),
  scenario_G = modifyList(true_example, list(phi = 0.2)),
  scenario_H = modifyList(true_example, list(phi = 0.4)),
  scenario_I = modifyList(true_example, list(phi = 0.6)),
  scenario_J = modifyList(true_example, list(phi = 0.8))
)

# Example multiple scenarios with incremental CSV saving:
res_all <- run_mc_scenarios_general(
  scenarios = scenarios_example,
  n_values = c(100,200,400,600,800,1000,2000),
  R = 1000,
  cmax1 = 132.4,
  cmax2 = 132.4,
  seed = 2026,
  n_cores = 20,
  output_dir = "simulation_results_paper_two_covariates",
  save_csv = TRUE,
  save_incremental = TRUE,
  block_size = 50
)



p <- plot_all_parameters_by_n_gg(res_all$overall_summary)
print(p)
ggplot2::ggsave("simulation_study2_copula_two_covariate.pdf", plot = p, width = 18, height = 24)




end.time <- Sys.time()
print(end.time - start.time)





## Updated plot 


get_valid_mc_results <- function(mc_res) {
  se_cols <- grep("^SE_", names(mc_res), value = TRUE)
  
  if (length(se_cols) == 0) {
    stop("No columns starting with 'SE_' were found in mc_res.")
  }
  
  keep <- rowSums(is.na(mc_res[, se_cols, drop = FALSE])) == 0
  mc_res[keep, , drop = FALSE]
}


summarize_mc_general_valid <- function(mc_res, true,
                                       scenario_id = NA_character_,
                                       n = NA_integer_,
                                       phi_true = NA_real_) {
  true_vec <- true_vector_general(true)
  pars <- names(true_vec)
  
  
  conv <- mc_res$converged == 1
  
  tab_par <- do.call(rbind, lapply(pars, function(p) {
    est <- mc_res[conv, p]
    se  <- mc_res[conv, paste0("SE_", p)]
    tru <- unname(true_vec[p])
    
    data.frame(
      Scenario = scenario_id,
      n = n,
      phi_true = phi_true,
      Parameter = p,
      Truth = tru,
      Mean_MLE = if (length(est) == 0 || all(is.na(est))) NA_real_ else mean(est, na.rm = TRUE),
      Bias = if (length(est) == 0 || all(is.na(est))) NA_real_ else mean(est - tru, na.rm = TRUE),
      SD_MLE = if (length(est) == 0 || all(is.na(est))) NA_real_ else stats::sd(est, na.rm = TRUE),
      Mean_SE = if (length(se) == 0 || all(is.na(se))) NA_real_ else mean(se, na.rm = TRUE),
      RMSE = if (length(est) == 0 || all(is.na(est))) NA_real_ else sqrt(mean((est - tru)^2, na.rm = TRUE)),
      Coverage_95 = if (length(est) == 0 || all(is.na(est)) || all(is.na(se))) {
        NA_real_
      } else {
        mean(est - 1.96 * se <= tru & tru <= est + 1.96 * se, na.rm = TRUE)
      },
      Convergence_Rate = mean(conv, na.rm = TRUE),
      Mean_Event_Rate_T1 = mean(mc_res$event_rate_T1, na.rm = TRUE),
      Mean_Event_Rate_T2 = mean(mc_res$event_rate_T2, na.rm = TRUE),
      N_Converged = sum(conv, na.rm = TRUE),
      R = nrow(mc_res)
    )
  }))
  
  rownames(tab_par) <- NULL
  tab_par
}



mc_valid <- get_valid_mc_results(res_all$all_mc_results)

overall_summary_valid <- do.call(
  rbind,
  lapply(split(mc_valid, mc_valid$Scenario), function(df_sc) {
    sc_id   <- unique(df_sc$Scenario)
    n_i     <- unique(df_sc$n)
    phi_i   <- unique(df_sc$phi_true)
    sc_name <- sub("_n[0-9]+$", "", sc_id)
    
    summarize_mc_general_valid(
      mc_res = df_sc,
      true = scenarios_example[[sc_name]],
      scenario_id = sc_id,
      n = n_i,
      phi_true = phi_i
    )
  })
)

rownames(overall_summary_valid) <- NULL

plot_all_parameters_by_n_gg <- function(summary_tab,
                                        parameters = unique(summary_tab$Parameter),
                                        phi_values = NULL,
                                        metrics = c("Bias", "RMSE", "Coverage_95")) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggh4x", quietly = TRUE)) stop("Package 'ggh4x' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")
  
  needed_cols <- c("n", "phi_true", "Parameter", metrics)
  miss <- setdiff(needed_cols, names(summary_tab))
  if (length(miss) > 0) stop("summary_tab is missing: ", paste(miss, collapse = ", "))
  
  dat <- summary_tab[summary_tab$Parameter %in% parameters, , drop = FALSE]
  
  if (!is.null(phi_values)) {
    dat <- dat[dat$phi_true %in% phi_values, , drop = FALSE]
  }
  
  keep_metric_row <- rowSums(is.na(dat[, metrics, drop = FALSE])) == 0
  dat <- dat[keep_metric_row, , drop = FALSE]
  
  long_dat <- do.call(
    rbind,
    lapply(metrics, function(m) {
      data.frame(
        n = dat$n,
        phi_true = dat$phi_true,
        Parameter = dat$Parameter,
        Metric = m,
        Value = dat[[m]]
      )
    })
  )
  
  long_dat <- long_dat[!is.na(long_dat$Value), , drop = FALSE]
  
  long_dat$n <- factor(long_dat$n, levels = sort(unique(long_dat$n)))
  
  parameter_labels <- c(
    alpha0 = "alpha[0]",
    alpha1 = "alpha[1]",
    alpha2 = "alpha[2]",
    b0     = "beta[0]",
    b1     = "beta[1]",
    b2     = "beta[2]",
    gamma0 = "gamma[0]",
    gamma1 = "gamma[1]",
    gamma2 = "gamma[2]",
    d0     = "lambda[0]",
    d1     = "lambda[1]",
    d2     = "lambda[2]",
    phi    = "phi"
  )
  
  long_dat$Parameter_raw <- long_dat$Parameter
  
  long_dat$Parameter <- factor(
    long_dat$Parameter,
    levels = names(parameter_labels),
    labels = parameter_labels
  )
  
  long_dat$Metric <- factor(
    long_dat$Metric,
    levels = c("Bias", "RMSE", "Coverage_95"),
    labels = c("Bias", "RMSE", "CP95%")
  )
  
  phi_levels <- sort(unique(long_dat$phi_true))
  long_dat$phi_true <- factor(
    long_dat$phi_true,
    levels = phi_levels,
    labels = paste0("phi = ", phi_levels)
  )
  
  p <- ggplot2::ggplot(
    long_dat,
    ggplot2::aes(
      x = n,
      y = Value,
      color = phi_true,
      group = phi_true
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggh4x::facet_grid2(
      Parameter ~ Metric,
      scales = "free_y",
      independent = "y",
      labeller = ggplot2::labeller(Parameter = ggplot2::label_parsed)
    ) +
    ggplot2::labs(
      x = "Sample size (n)",
      y = NULL,
      color = "Scenario"
    ) +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("CP95%", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0.93649),
      inherit.aes = FALSE,
      linetype = 2
    ) +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("CP95%", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0.9635),
      inherit.aes = FALSE,
      linetype = 2
    ) +
    ggplot2::geom_hline(
      data = data.frame(Metric = factor("Bias", levels = levels(long_dat$Metric))),
      ggplot2::aes(yintercept = 0),
      inherit.aes = FALSE,
      linetype = 2
    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        Parameter == "phi" & Metric == "CP95%" ~
          ggplot2::scale_y_continuous(
            trans = "logit",
            breaks = c(0.6, 0.8, 0.90, 0.95, 0.98)
          )
      )
    )
  
  return(p)
}

p_valid <- plot_all_parameters_by_n_gg(overall_summary_valid)
print(p_valid)

ggplot2::ggsave(
  "simulation_study2_copula_two_covariate.pdf",
  plot = p_valid,
  width = 18,
  height = 24
)


library(dplyr)

se_cols <- grep("^SE_", names(res_all$all_mc_results), value = TRUE)

removal_summary <- res_all$all_mc_results %>%
  mutate(valid_rep = rowSums(is.na(across(all_of(se_cols)))) == 0) %>%
  group_by(Scenario, n, phi_true) %>%
  summarise(
    R_total = n(),
    R_valid = sum(valid_rep),
    R_removed = sum(!valid_rep),
    prop_removed = mean(!valid_rep),
    .groups = "drop"
  )

print(removal_summary)




plot_phi_cp95_gg <- function(summary_tab, phi_values = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")
  
  dat <- subset(summary_tab, Parameter == "phi")
  
  if (!is.null(phi_values)) {
    dat <- dat[dat$phi_true %in% phi_values, , drop = FALSE]
  }
  
  dat <- dat[!is.na(dat$Coverage_95), , drop = FALSE]
  
  dat$n <- factor(dat$n, levels = sort(unique(dat$n)))
  
  phi_levels <- sort(unique(dat$phi_true))
  dat$phi_true <- factor(
    dat$phi_true,
    levels = phi_levels,
    labels = paste0("phi = ", phi_levels)
  )
  
  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = n,
      y = Coverage_95,
      color = phi_true,
      group = phi_true
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = 0.93649, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0.9635, linetype = 2) +
    ggplot2::scale_y_continuous(
      trans = "logit",
      breaks = c(0.6,0.7,0.8,0.85,0.90, 0.93, 0.95, 0.97, 0.99)
    ) +
    ggplot2::labs(
      x = "Sample size (n)",
      y = "CP95%",
      color = "Scenario"
    ) +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

p_phi_cp95 <- plot_phi_cp95_gg(overall_summary_valid)
print(p_phi_cp95)

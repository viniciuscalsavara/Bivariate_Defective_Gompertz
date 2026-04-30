############################################################
## Bivariate Gompertz-FGM model with covariate in shape   ##
## and scale parameters + Uniform(0,25) censoring         ##
## Includes cure-fraction estimates by covariate level    ##
## and delta-method standard errors                       ##
############################################################

rm(list=ls(all=TRUE))
start.time <- Sys.time()


############################
## 1) Utilities           ##
############################

safe_log <- function(z, eps = 1e-12) log(pmax(z, eps))
clip01   <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

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

############################
## 3) Cure fraction       ##
############################
cure_gompertz <- function(a, scale, eps = 1e-10) {
  a <- unname(a)
  scale <- unname(scale)
  
  if (!is.finite(a) || !is.finite(scale)) return(NA_real_)
  if (a < -eps) return(exp(scale / a))
  return(0)
}

grad_cure_a_eta <- function(a, eta, eps = 1e-10) {
  a <- unname(a)
  eta <- unname(eta)
  
  if (!is.finite(a) || !is.finite(eta)) {
    return(c(da = NA_real_, deta = NA_real_))
  }
  
  scale <- exp(eta)
  
  if (!is.finite(scale)) {
    return(c(da = NA_real_, deta = NA_real_))
  }
  
  if (a >= -eps) {
    return(c(da = 0, deta = 0))
  }
  
  p <- exp(scale / a)
  
  if (!is.finite(p)) {
    return(c(da = NA_real_, deta = NA_real_))
  }
  
  c(
    da   = p * (-scale / (a^2)),
    deta = p * ( scale / a )
  )
}



get_true_cure <- function(true, eps = 1e-10) {
  a10 <- true["alpha0"]
  eta10 <- true["beta0"]
  sc10 <- exp(eta10)
  
  a11 <- true["alpha0"] + true["alpha1"]
  eta11 <- true["beta0"] + true["beta1"]
  sc11 <- exp(eta11)
  
  a20 <- true["gamma0"]
  eta20 <- true["lambda0"]
  sc20 <- exp(eta20)
  
  a21 <- true["gamma0"] + true["gamma1"]
  eta21 <- true["lambda0"] + true["lambda1"]
  sc21 <- exp(eta21)
  
  c(
    cure_T1_x0 = cure_gompertz(a10, sc10, eps),
    cure_T1_x1 = cure_gompertz(a11, sc11, eps),
    cure_T2_x0 = cure_gompertz(a20, sc20, eps),
    cure_T2_x1 = cure_gompertz(a21, sc21, eps)
  )
}



get_cure_estimates <- function(est, vc, eps = 1e-10) {
  out_names <- c(
    "cure_T1_x0", "cure_T1_x1", "cure_T2_x0", "cure_T2_x1",
    "SE_cure_T1_x0", "SE_cure_T1_x1", "SE_cure_T2_x0", "SE_cure_T2_x1"
  )
  
  if (is.null(est) || is.null(vc)) {
    out <- setNames(rep(NA_real_, length(out_names)), out_names)
    return(out)
  }
  
  ## garante nomes corretos
  est <- unname(est)
  names(est) <- c("alpha0", "alpha1", "beta0", "beta1",
                  "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  
  ## função auxiliar para delta method usando apenas o sub-bloco relevante
  se_delta <- function(L, idx_eta, vc_full, cure_value, a_value, eta_value) {
    cure_value <- unname(cure_value)
    a_value    <- unname(a_value)
    eta_value  <- unname(eta_value)
    
    if (!is.finite(cure_value) || !is.finite(a_value) || !is.finite(eta_value)) {
      return(NA_real_)
    }
    
    g_small <- grad_cure_a_eta(a_value, eta_value, eps)
    
    if (any(!is.finite(g_small))) {
      return(NA_real_)
    }
    
    ## gradiente em relação aos parâmetros originais:
    ## grad = L %*% (dp/da, dp/deta)
    grad_full <- as.vector(L %*% g_small)
    V_sub <- vc_full[idx_eta, idx_eta, drop = FALSE]
    
    if (any(!is.finite(V_sub))) {
      return(NA_real_)
    }
    
    val <- drop(t(grad_full) %*% V_sub %*% grad_full)
    
    if (!is.finite(val) || val < 0) {
      return(NA_real_)
    }
    
    sqrt(val)
  }
  
  ## ========= T1, x = 0 =========
  a10   <- unname(est["alpha0"])
  eta10 <- unname(est["beta0"])
  sc10  <- exp(eta10)
  cure10 <- unname(cure_gompertz(a10, sc10, eps))
  
  ## a10 = alpha0
  ## eta10 = beta0
  idx10 <- c("alpha0", "alpha1", "beta0", "beta1")
  L10 <- rbind(
    c(1, 0),
    c(0, 0),
    c(0, 1),
    c(0, 0)
  )
  se_cure10 <- se_delta(L10, idx10, vc, cure10, a10, eta10)
  
  ## ========= T1, x = 1 =========
  a11   <- unname(est["alpha0"] + est["alpha1"])
  eta11 <- unname(est["beta0"] + est["beta1"])
  sc11  <- exp(eta11)
  cure11 <- unname(cure_gompertz(a11, sc11, eps))
  
  ## a11 = alpha0 + alpha1
  ## eta11 = beta0 + beta1
  idx11 <- c("alpha0", "alpha1", "beta0", "beta1")
  L11 <- rbind(
    c(1, 0),
    c(1, 0),
    c(0, 1),
    c(0, 1)
  )
  se_cure11 <- se_delta(L11, idx11, vc, cure11, a11, eta11)
  
  ## ========= T2, x = 0 =========
  a20   <- unname(est["gamma0"])
  eta20 <- unname(est["lambda0"])
  sc20  <- exp(eta20)
  cure20 <- unname(cure_gompertz(a20, sc20, eps))
  
  idx20 <- c("gamma0", "gamma1", "lambda0", "lambda1")
  L20 <- rbind(
    c(1, 0),
    c(0, 0),
    c(0, 1),
    c(0, 0)
  )
  se_cure20 <- se_delta(L20, idx20, vc, cure20, a20, eta20)
  
  ## ========= T2, x = 1 =========
  a21   <- unname(est["gamma0"] + est["gamma1"])
  eta21 <- unname(est["lambda0"] + est["lambda1"])
  sc21  <- exp(eta21)
  cure21 <- unname(cure_gompertz(a21, sc21, eps))
  
  idx21 <- c("gamma0", "gamma1", "lambda0", "lambda1")
  L21 <- rbind(
    c(1, 0),
    c(1, 0),
    c(0, 1),
    c(0, 1)
  )
  se_cure21 <- se_delta(L21, idx21, vc, cure21, a21, eta21)
  
  out <- c(
    cure_T1_x0 = cure10,
    cure_T1_x1 = cure11,
    cure_T2_x0 = cure20,
    cure_T2_x1 = cure21,
    SE_cure_T1_x0 = se_cure10,
    SE_cure_T1_x1 = se_cure11,
    SE_cure_T2_x0 = se_cure20,
    SE_cure_T2_x1 = se_cure21
  )
  
  out
}


############################
## 4) Joint survival      ##
############################

S12 <- function(alpha0, alpha1,beta0, beta1,
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
                gamma0, gamma1,lambda0, lambda1, phi,
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
                gamma0, gamma1, lambda0, lambda1,phi,
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
## 5) Negative loglik     ##
############################

nll_txa <- function(par, t1, t2, cens1, cens2, x, eps = 1e-12) {
  alpha0 <- par[1]
  alpha1 <- par[2]
  beta0     <- par[3]
  beta1     <- par[4]
  gamma0 <- par[5]
  gamma1 <- par[6]
  lambda0     <- par[7]
  lambda1     <- par[8]
  phi    <- tanh(par[9])
  
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
## 6) Simulation          ##
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
  
  if (is.null(x)) x <- rbinom(n, 1, 0.406)  ##Female proportion
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
  
  data.frame(
    x = x,
    t1 = t1,
    t2 = t2,
    cens1 = cens1,
    cens2 = cens2
  )
}

############################
## 7) Estimation          ##
############################

fit_model <- function(
    dat,
    init = c(alpha0 = -0.70, alpha1 =  -0.1,
             beta0     = log(0.3), beta1 =  0.2,
             gamma0 = -0.1, gamma1 = 0.05,
             lambda0     = log(0.25), lambda1 =  0.05,
             phi_z  = atanh(0.10))
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
  
  par_names <- c("alpha0", "alpha1", "beta0", "beta1",
                 "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  cure_names <- c("cure_T1_x0", "cure_T1_x1", "cure_T2_x0", "cure_T2_x1",
                  "SE_cure_T1_x0", "SE_cure_T1_x1", "SE_cure_T2_x0", "SE_cure_T2_x1")
  
  if (inherits(fit, "try-error")) {
    est <- setNames(rep(NA_real_, 9), par_names)
    se  <- setNames(rep(NA_real_, 9), paste0("SE_", par_names))
    cure_est <- setNames(rep(NA_real_, 8), cure_names)
    return(list(converged = FALSE, est = est, se = se, vc = NULL,
                cure = cure_est, logLik = NA_real_, fit = NULL))
  }
  
  est_int <- fit$par
 
  est_nat <- c(
    unname(est_int[1]),
    unname(est_int[2]),
    unname(est_int[3]),
    unname(est_int[4]),
    unname(est_int[5]),
    unname(est_int[6]),
    unname(est_int[7]),
    unname(est_int[8]),
    unname(tanh(est_int[9]))
  )
  names(est_nat) <- c("alpha0", "alpha1", "beta0", "beta1",
                      "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  
  
  vc_int <- try(solve(fit$hessian), silent = TRUE)
  
  if (inherits(vc_int, "try-error") || any(!is.finite(vc_int))) {
    vc_nat <- matrix(NA_real_, 9, 9, dimnames = list(par_names, par_names))
    se_nat <- setNames(rep(NA_real_, 9), paste0("SE_", par_names))
    cure_est <- setNames(rep(NA_real_, 8), cure_names)
  } else {
    J <- diag(c(1, 1, 1, 1, 1, 1, 1, 1, 1 - tanh(est_int[9])^2))
    vc_nat <- J %*% vc_int %*% t(J)
    colnames(vc_nat) <- rownames(vc_nat) <- par_names
    se_nat <- setNames(sqrt(diag(vc_nat)), paste0("SE_", par_names))
    cure_est <- get_cure_estimates(est = est_nat, vc = vc_nat)
  }
  
  list(
    converged = fit$convergence == 0,
    est = est_nat,
    se = se_nat,
    vc = vc_nat,
    cure = cure_est,
    logLik = -fit$value,
    fit = fit
  )
}

############################
## 8) Monte Carlo         ##
############################

one_rep <- function(n, true, cmax1 = 130, cmax2 = 130, seed = NULL) {
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
  
  fit <- fit_model(dat)
  
  data.frame(
    converged = as.integer(fit$converged),
    
    alpha0 = fit$est["alpha0"],
    alpha1 = fit$est["alpha1"],
    beta0     = fit$est["beta0"],
    beta1     = fit$est["beta1"],
    gamma0 = fit$est["gamma0"],
    gamma1 = fit$est["gamma1"],
    lambda0     = fit$est["lambda0"],
    lambda1     = fit$est["lambda1"],
    phi    = fit$est["phi"],
    
    SE_alpha0 = fit$se["SE_alpha0"],
    SE_alpha1 = fit$se["SE_alpha1"],
    SE_beta0     = fit$se["SE_beta0"],
    SE_beta1     = fit$se["SE_beta1"],
    SE_gamma0 = fit$se["SE_gamma0"],
    SE_gamma1 = fit$se["SE_gamma1"],
    SE_lambda0     = fit$se["SE_lambda0"],
    SE_lambda1     = fit$se["SE_lambda1"],
    SE_phi    = fit$se["SE_phi"],
    
    cure_T1_x0 = fit$cure["cure_T1_x0"],
    cure_T1_x1 = fit$cure["cure_T1_x1"],
    cure_T2_x0 = fit$cure["cure_T2_x0"],
    cure_T2_x1 = fit$cure["cure_T2_x1"],
    
    SE_cure_T1_x0 = fit$cure["SE_cure_T1_x0"],
    SE_cure_T1_x1 = fit$cure["SE_cure_T1_x1"],
    SE_cure_T2_x0 = fit$cure["SE_cure_T2_x0"],
    SE_cure_T2_x1 = fit$cure["SE_cure_T2_x1"],
    
    cure_rate_T1 = mean(dat$cure1),
    cure_rate_T2 = mean(dat$cure2),
    event_rate_T1 = mean(dat$cens1),
    event_rate_T2 = mean(dat$cens2)
  )
}

mc_study_parallel <- function(R, n, true,
                              cmax1 = 130, cmax2 = 130,
                              seed = 123,
                              n_cores = max(1, parallel::detectCores() - 1),
                              type = if (.Platform$OS.type == "windows") "PSOCK" else "PSOCK") {
  if (n_cores < 1) n_cores <- 1
  n_cores <- min(n_cores, R)
  
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, R)
  
  cl <- parallel::makeCluster(n_cores, type = type)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c("safe_log", "clip01", "Sg", "fg", "Fg", "pcure",
                "cure_gompertz", "grad_cure_a_eta", "get_true_cure", "get_cure_estimates",
                "S12", "dS1", "dS2", "d2S12", "nll_txa",
                "qSg", "rfgm", "sim_data", "fit_model", "one_rep"),
    envir = environment()
  )
  
  res <- parallel::parLapply(
    cl,
    X = seq_len(R),
    fun = function(r, n, true, cmax1, cmax2, seeds) {
      one_rep(n = n, true = true, cmax1 = cmax1, cmax2 = cmax2, seed = seeds[r])
    },
    n = n, true = true, cmax1 = cmax1, cmax2 = cmax2, seeds = seeds
  )
  
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

############################
## 9) Monte Carlo summary ##
############################

summarize_mc <- function(mc_res, true, scenario_id = NA_character_,
                         n = NA_integer_, phi_true = NA_real_) {
  pars <- c("alpha0", "alpha1", "beta0", "beta1",
            "gamma0", "gamma1", "lambda0", "lambda1", "phi")
  cure_pars <- c("cure_T1_x0", "cure_T1_x1", "cure_T2_x0", "cure_T2_x1")
  cure_true <- get_true_cure(true)
  
  tab_par <- do.call(rbind, lapply(pars, function(p) {
    est <- mc_res[[p]][mc_res$converged == 1]
    se  <- mc_res[[paste0("SE_", p)]][mc_res$converged == 1]
    tru <- unname(true[p])
    
    data.frame(
      Scenario = scenario_id,
      n = n,
      phi_true = phi_true,
      Parameter = p,
      Truth = tru,
      Mean_MLE = mean(est, na.rm = TRUE),
      Bias = mean(est - tru, na.rm = TRUE),
      SD_MLE = sd(est, na.rm = TRUE),
      Mean_SE = mean(se, na.rm = TRUE),
      RMSE = sqrt(mean((est - tru)^2, na.rm = TRUE)),
      Coverage_95 = mean(est - 1.96 * se <= tru & tru <= est + 1.96 * se, na.rm = TRUE),
      Convergence_Rate = mean(mc_res$converged == 1, na.rm = TRUE)
   #   Mean_Cure_Rate_T1 = mean(mc_res$cure_rate_T1, na.rm = TRUE),
  #    Mean_Cure_Rate_T2 = mean(mc_res$cure_rate_T2, na.rm = TRUE),
   #   Mean_Event_Rate_T1 = mean(mc_res$event_rate_T1, na.rm = TRUE),
  #    Mean_Event_Rate_T2 = mean(mc_res$event_rate_T2, na.rm = TRUE)
    )
  }))
  
  tab_cure <- do.call(rbind, lapply(cure_pars, function(p) {
    est <- mc_res[[p]][mc_res$converged == 1]
    se  <- mc_res[[paste0("SE_", p)]][mc_res$converged == 1]
    tru <- unname(cure_true[p])
    
    data.frame(
      Scenario = scenario_id,
      n = n,
      phi_true = phi_true,
      Parameter = p,
      Truth = tru,
      Mean_MLE = mean(est, na.rm = TRUE),
      Bias = mean(est - tru, na.rm = TRUE),
      SD_MLE = sd(est, na.rm = TRUE),
      Mean_SE = mean(se, na.rm = TRUE),
      RMSE = sqrt(mean((est - tru)^2, na.rm = TRUE)),
      Coverage_95 = mean(est - 1.96 * se <= tru & tru <= est + 1.96 * se, na.rm = TRUE),
      Convergence_Rate = mean(mc_res$converged == 1, na.rm = TRUE)
     # Mean_Cure_Rate_T1 = mean(mc_res$cure_rate_T1, na.rm = TRUE),
    #  Mean_Cure_Rate_T2 = mean(mc_res$cure_rate_T2, na.rm = TRUE),
    #  Mean_Event_Rate_T1 = mean(mc_res$event_rate_T1, na.rm = TRUE),
    #  Mean_Event_Rate_T2 = mean(mc_res$event_rate_T2, na.rm = TRUE)
    )
  }))
  
  out <- rbind(tab_par, tab_cure)
  rownames(out) <- NULL
  out
}

############################
## 10) Plots              ##
############################

plot_parameter_histograms <- function(mc_res, true,
                                      params = c("alpha0", "alpha1", "beta0", "beta1",
                                                 "gamma0", "gamma1", "lambda0", "lambda1", "phi"),
                                      converged_only = TRUE,
                                      file = NULL) {
  if (converged_only) mc_res <- mc_res[mc_res$converged == 1, , drop = FALSE]
  if (!is.null(file)) grDevices::pdf(file, width = 10, height = 8)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); if (!is.null(file)) grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(3, 3))
  for (p in params) {
    hist(mc_res[[p]], main = paste("Histogram of", p), xlab = p, breaks = "FD")
    abline(v = true[p], lwd = 2, lty = 2)
  }
}

plot_scenario_summary <- function(summary_tab,
                                  metric = c("Bias", "RMSE", "Coverage_95", "SD_MLE", "Mean_SE"),
                                  parameter = "phi",
                                  file = NULL) {
  metric <- match.arg(metric)
  dat <- summary_tab[summary_tab$Parameter == parameter, , drop = FALSE]
  n_levels <- sort(unique(dat$n))
  phi_levels <- sort(unique(dat$phi_true))
  mat <- sapply(phi_levels, function(ph) {
    vals <- dat[match(n_levels, dat$n[dat$phi_true == ph]), metric]
    vals
  })
  if (!is.null(file)) grDevices::pdf(file, width = 9, height = 6)
  on.exit(if (!is.null(file)) grDevices::dev.off(), add = TRUE)
  matplot(n_levels, mat, type = "b", pch = 1:ncol(mat), lty = 1:ncol(mat),
          xlab = "Sample size (n)", ylab = metric,
          main = paste(metric, "for", parameter))
  legend("topright", legend = paste0("phi=", phi_levels),
         col = 1:ncol(mat), pch = 1:ncol(mat), lty = 1:ncol(mat), bty = "n")
}

plot_overall_summary_grid <- function(summary_tab,
                                      parameters = c("alpha0", "alpha1", "beta0", "beta1",
                                                     "gamma0", "gamma1", "lambda0", "lambda1", "phi",
                                                     "cure_T1_x0", "cure_T1_x1", "cure_T2_x0", "cure_T2_x1"),
                                      metric = c("RMSE", "Bias"), file = NULL) {
  metric <- match.arg(metric)
  if (!is.null(file)) grDevices::pdf(file, width = 14, height = 10)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); if (!is.null(file)) grDevices::dev.off()}, add = TRUE)
  nr <- ceiling(sqrt(length(parameters)))
  nc <- ceiling(length(parameters) / nr)
  par(mfrow = c(nr, nc))
  for (p in parameters) {
    dat <- summary_tab[summary_tab$Parameter == p, , drop = FALSE]
    n_levels <- sort(unique(dat$n))
    phi_levels <- sort(unique(dat$phi_true))
    mat <- sapply(phi_levels, function(ph) dat[match(n_levels, dat$n[dat$phi_true == ph]), metric])
    matplot(n_levels, mat, type = "b", pch = 1:ncol(mat), lty = 1:ncol(mat),
            xlab = "n", ylab = metric, main = p)
  }
}

############################
## 11) Scenario engine    ##
############################

run_scenarios <- function(R,
                          true_base,
                          n_values = c(100, 200, 400, 600, 800, 1000),
                          phi_values = c(0.2, 0.4, 0.6, 0.8),
                          cmax1 = 25,
                          cmax2 = 25,
                          seed = 123,
                          n_cores = max(1, parallel::detectCores() - 1),
                          output_dir = "simulation_results") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  scenarios <- expand.grid(
    n = n_values,
    phi = phi_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  all_mc <- vector("list", nrow(scenarios))
  all_sum <- vector("list", nrow(scenarios))
  
  for (i in seq_len(nrow(scenarios))) {
    n_i <- scenarios$n[i]
    phi_i <- scenarios$phi[i]
    scenario_id <- paste0("n", n_i, "_phi", gsub("\\.", "", sprintf("%.1f", phi_i)))
    
    true_i <- true_base
    true_i["phi"] <- phi_i
    
    cat("Running scenario:", scenario_id, "\n")
    
    mc_res_i <- mc_study_parallel(
      R = R, n = n_i, true = true_i,
      cmax1 = cmax1, cmax2 = cmax2,
      seed = seed + i, n_cores = n_cores
    )
    
    mc_res_i$Scenario <- scenario_id
    mc_res_i$n <- n_i
    mc_res_i$phi_true <- phi_i
    
    sum_i <- summarize_mc(
      mc_res_i, true = true_i,
      scenario_id = scenario_id, n = n_i, phi_true = phi_i
    )
    
    utils::write.csv(mc_res_i,
                     file.path(output_dir, paste0("mc_results_", scenario_id, ".csv")),
                     row.names = FALSE)
    utils::write.csv(sum_i,
                     file.path(output_dir, paste0("summary_", scenario_id, ".csv")),
                     row.names = FALSE)
    
    plot_parameter_histograms(
      mc_res_i, true = true_i,
      file = file.path(output_dir, paste0("histograms_", scenario_id, ".pdf"))
    )
    
    all_mc[[i]] <- mc_res_i
    all_sum[[i]] <- sum_i
  }
  
  all_mc_df <- do.call(rbind, all_mc)
  all_sum_df <- do.call(rbind, all_sum)
  rownames(all_mc_df) <- NULL
  rownames(all_sum_df) <- NULL
  
  utils::write.csv(all_mc_df,
                   file.path(output_dir, "all_mc_results.csv"),
                   row.names = FALSE)
  utils::write.csv(all_sum_df,
                   file.path(output_dir, "overall_summary_for_article.csv"),
                   row.names = FALSE)
  
  plot_scenario_summary(all_sum_df, metric = "RMSE", parameter = "phi",
                        file = file.path(output_dir, "plot_RMSE_phi.pdf"))
  plot_scenario_summary(all_sum_df, metric = "Bias", parameter = "phi",
                        file = file.path(output_dir, "plot_Bias_phi.pdf"))
  plot_scenario_summary(all_sum_df, metric = "Coverage_95", parameter = "phi",
                        file = file.path(output_dir, "plot_Coverage_phi.pdf"))
  plot_overall_summary_grid(all_sum_df, metric = "RMSE",
                            file = file.path(output_dir, "plot_overall_RMSE_grid.pdf"))
  plot_overall_summary_grid(all_sum_df, metric = "Bias",
                            file = file.path(output_dir, "plot_overall_Bias_grid.pdf"))
  
  list(all_mc_results = all_mc_df, overall_summary = all_sum_df)
}

############################
## 12) Optional ggplot    ##
############################
plot_all_parameters_by_n_gg <- function(summary_tab,
                                        parameters = c("alpha0", "alpha1", "beta0", "beta1",
                                                       "gamma0", "gamma1", "lambda0", "lambda1", "phi",
                                                       "cure_T1_x0", "cure_T1_x1",
                                                       "cure_T2_x0", "cure_T2_x1"),
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
  
  param_labels <- c(
    alpha0 = "alpha[0]",
    alpha1 = "alpha[1]",
    beta0 = "beta[0]",
    beta1 = "beta[1]",
    gamma0 = "gamma[0]",
    gamma1 = "gamma[1]",
    lambda0 = "lambda[0]",
    lambda1 = "lambda[1]",
    phi = "phi",
    cure_T1_x0 = "p[10]",
    cure_T1_x1 = "p[11]",
    cure_T2_x0 = "p[20]",
    cure_T2_x1 = "p[21]"
  )
  
  long_dat$Parameter <- factor(long_dat$Parameter, levels = parameters)
  
  long_dat$Metric <- factor(
    long_dat$Metric,
    levels = c("Bias", "RMSE", "Coverage_95"),
    labels = c("Bias", "RMSE", "CP95%")
  )
  
  long_dat$phi_true <- factor(
    long_dat$phi_true,
    levels = sort(unique(long_dat$phi_true)),
    labels = paste0("phi = ", sort(unique(long_dat$phi_true)))
  )
  
  long_dat$n <- factor(long_dat$n, levels = sort(unique(long_dat$n)))
  
  ggplot2::ggplot(
    long_dat,
    ggplot2::aes(x = n, y = Value, color = phi_true, group = phi_true)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggh4x::facet_grid2(
      Parameter ~ Metric,
      scales = "free_y",
      independent = "y",
      labeller = ggplot2::labeller(
        Parameter = ggplot2::as_labeller(param_labels, default = ggplot2::label_parsed)
      )
    ) +
    ggplot2::labs(x = "Sample size (n)", y = NULL, color = "Scenario") +
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 22),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 22),
      legend.text = ggplot2::element_text(size = 22),
      axis.title = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 18),
      panel.grid.minor = ggplot2::element_blank()
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
}

############################
## 13) Example            ##
############################

true_par <- c(
  alpha0 = -0.7077,
  alpha1 =  -0.1669,
  beta0     = -1.2404,
  beta1     =  0.4399,
  gamma0 = -0.056,
  gamma1 =  0.0025,
  lambda0     = -2.363,
  lambda1     =  -0.0504,
  phi    =  -0.436
  
)

## One fit example
set.seed(123)
dat <- sim_data(
  n = 2000,
  alpha0 = true_par["alpha0"],
  alpha1 = true_par["alpha1"],
  beta0     = true_par["beta0"],
  beta1     = true_par["beta1"],
  gamma0 = true_par["gamma0"],
  gamma1 = true_par["gamma1"],
  lambda0     = true_par["lambda0"],
  lambda1     = true_par["lambda1"],
  phi    = true_par["phi"],
  cmax1  = 132.4,
  cmax2  = 132.4
)

fit <- fit_model(dat)
print(round(fit$est, 4))
print(round(fit$se, 4))
print(round(fit$cure, 4))
print(fit$converged)



## Full scenario study
 scenario_out <- run_scenarios(
   R = 1000,
   true_base = true_par,
   n_values = c(100,200,400,600,800,1000,2000),
   phi_values = c(-0.8,-0.6,-0.436,-0.2,0,0.2,0.4,0.6,0.8),
   cmax1 = 132.4,
   cmax2 = 132.4,
   seed = 321,
   n_cores = 20,
   output_dir = "simulation_results_paper_one_covariate"
 )
 head(scenario_out$overall_summary)

 p <- plot_all_parameters_by_n_gg(scenario_out$overall_summary)
 print(p)
 ggplot2::ggsave("simulation_study1_copula_one_covariate.pdf", plot = p, width = 18, height = 24)

 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
 
  
## Plots ##

library(readr)
overall_summary <- read_csv("simulation_results_paper_one_covariate/overall_summary_for_article.csv")
p <- plot_all_parameters_by_n_gg(overall_summary)
print(p)
ggplot2::ggsave("simulation_study1_copula_one_covariate.pdf", plot = p, width = 18, height = 24)

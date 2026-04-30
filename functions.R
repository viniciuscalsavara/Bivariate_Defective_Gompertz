############################################################
## Bivariate Gompertz-FGM model ##
############################################################

## Parameterization
## T1: a1(x) = alpha0 + alpha'x, b1(x) = beta1'x
## T2: a2(x) = gamma0 + gamma'x, b2(x) = beta2'x
## Dependence: phi (FGM copula)

############################
## 1) Utilities           ##
############################

safe_log <- function(z, eps = 1e-12) log(pmax(z, eps))
clip01   <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

time_var = NULL
event_var = NULL

############################
## 2) Core model pieces   ##
###########################
surv_g <- function(a, b, t, eps = 1e-10) {
  ifelse(abs(a) < eps,
         exp(-b * t),
         exp(-(b / a) * (exp(a * t) - 1)))
}

f_g <- function(a, b, t, eps = 1e-10) {
  ifelse(abs(a) < eps,
         b * exp(-b * t),
         b * exp(a * t) * exp(-(b / a) * (exp(a * t) - 1)))
}

dist_g <- function(a, b, t, eps = 1e-10) 1 - surv_g(a, b, t, eps)

pcure <- function(a, b, eps = 1e-10) ifelse(a < -eps, exp(b / a), 0)

s12 <- function(a1, a2, b1, b2, phi, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  s1 * s2 * (1 + phi * (1 - s1) * (1 - s2))
}

ds1 <- function(a1, a2, b1, b2, phi, t1, t2, eps = 1e-10) {
  -f_g(a1, b1, t1, eps) * surv_g(a2, b2, t2, eps) +
    phi * surv_g(a2, b2, t2, eps) * dist_g(a2, b2, t2, eps) * f_g(a1, b1, t1, eps) * (1 - 2 * dist_g(a1, b1, t1, eps))
}

ds2 <- function(a1, a2, b1, b2, phi, t1, t2, eps = 1e-10) {
  -f_g(a2, b2, t2, eps) * surv_g(a1, b1, t1, eps) +
    phi * surv_g(a1, b1, t1, eps) * dist_g(a1, b1, t1, eps) * f_g(a2, b2, t2, eps) * (1 - 2 * dist_g(a2, b2, t2, eps))
}

ds12 <- function(a1, a2, b1, b2, phi, t1, t2, eps = 1e-10) {
  f_g(a1, b1, t1, eps) * f_g(a2, b2, t2, eps) +
    phi * f_g(a1, b1, t1, eps) * f_g(a2, b2, t2, eps) *
    (1 - 2 * dist_g(a1, b1, t1, eps)) * (1 - 2 * dist_g(a2, b2, t2, eps))
}


####################### Likelihood function #################################################
lik_dg_fgm_mv <- function(par, data, cens1, cens2, linear_pred_a1, linear_pred_a2,
                          linear_pred_b1, linear_pred_b2, eps = 1e-12) {
  # linear_pred_a1 for a1
  # linear_pred_a2 for a2
  # linear_pred_b1 for b1
  # linear_pred_b2 para b2
  
  Xa1 <- stats::model.matrix(linear_pred_a1, data)
  Xa2 <- stats::model.matrix(linear_pred_a2, data)
  
  Xb1 <- stats::model.matrix(linear_pred_b1, data)
  Xb2 <- stats::model.matrix(linear_pred_b2, data)
  
  j_a1 <- ncol(Xa1)
  j_a2 <- ncol(Xa2)
  j_b1 <- ncol(Xb1)
  j_b2 <- ncol(Xb2)
  
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  if (length(par) != npar) return(1e20)
  
  idx <- 0
  take_par <- function(k) {
    if (k == 0) return(numeric(0))
    out <- par[(idx + 1):(idx + k)]
    idx <<- idx + k
    out
  }
  
  reg_a1 <- take_par(j_a1)
  reg_b1 <- take_par(j_b1)
  reg_a2 <- take_par(j_a2)
  reg_b2 <- take_par(j_b2)
  phi <- tanh(take_par(1))
  
  cens1 <- as.numeric(data[[cens1]])
  cens2 <- as.numeric(data[[cens2]])
  
  t1 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a1, data)))
  t2 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a2, data)))
  
  a1 <- as.numeric(Xa1 %*% reg_a1)
  a2 <- as.numeric(Xa2 %*% reg_a2)
  
  b1 <- exp(as.numeric(Xb1 %*% reg_b1) )
  b2 <- exp(as.numeric(Xb2 %*% reg_b2))
  
  v_ds12 <- ds12(a1, a2, b1, b2, phi, t1, t2, eps)
  v_mds1 <- -ds1(a1, a2, b1, b2, phi, t1, t2, eps)
  v_mds2 <- -ds2(a1, a2, b1, b2, phi, t1, t2, eps)
  v_s12  <- s12(a1, a2, b1, b2, phi, t1, t2, eps)
  
  #  Needs non-negative terms; allow tiny negatives from
  # floating-point error and clip with safe_log.
  tol_neg <- -sqrt(eps)
  if (any(!is.finite(v_ds12)) || any(!is.finite(v_mds1)) ||
      any(!is.finite(v_mds2)) || any(!is.finite(v_s12)) ||
      any(v_ds12 < tol_neg) || any(v_mds1 < tol_neg) ||
      any(v_mds2 < tol_neg) || any(v_s12 < tol_neg)) {
    return(1e20)
  }
  
  ll <- sum(cens1 * cens2 * safe_log(v_ds12, eps)) +
    sum(cens1 * (1 - cens2) * safe_log(v_mds1, eps)) +
    sum(cens2 * (1 - cens1) * safe_log(v_mds2, eps)) +
    sum((1 - cens1) * (1 - cens2) * safe_log(v_s12, eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
  
}

####################### Likelihood function phi = 0 #####################################
lik_dg_phi0_mv <- function(par, data, cens1, cens2, linear_pred_a1, linear_pred_a2,
                           linear_pred_b1, linear_pred_b2, eps = 1e-12) {
  # linear_pred_a1 for a1
  # linear_pred_a2 for a2
  # linear_pred_b1 for b1
  # linear_pred_b2 para b2
  
  Xa1 <- stats::model.matrix(linear_pred_a1, data)
  Xa2 <- stats::model.matrix(linear_pred_a2, data)
  
  Xb1 <- stats::model.matrix(linear_pred_b1, data)
  Xb2 <- stats::model.matrix(linear_pred_b2, data)
  
  j_a1 <- ncol(Xa1)
  j_a2 <- ncol(Xa2)
  j_b1 <- ncol(Xb1)
  j_b2 <- ncol(Xb2)
  
  npar <- j_a1 + j_b1 + j_a2 + j_b2 
  if (length(par) != npar) return(1e20)
  
  idx <- 0
  take_par <- function(k) {
    if (k == 0) return(numeric(0))
    out <- par[(idx + 1):(idx + k)]
    idx <<- idx + k
    out
  }
  
  reg_a1 <- take_par(j_a1)
  reg_b1 <- take_par(j_b1)
  reg_a2 <- take_par(j_a2)
  reg_b2 <- take_par(j_b2)
  phi <- 0
  
  cens1 <- as.numeric(data[[cens1]])
  cens2 <- as.numeric(data[[cens2]])
  
  t1 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a1, data)))
  t2 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a2, data)))
  
  a1 <- as.numeric(Xa1 %*% reg_a1)
  a2 <- as.numeric(Xa2 %*% reg_a2)
  
  b1 <- exp(as.numeric(Xb1 %*% reg_b1) )
  b2 <- exp(as.numeric(Xb2 %*% reg_b2))
  
  v_ds12 <- ds12(a1, a2, b1, b2, phi, t1, t2, eps)
  v_mds1 <- -ds1(a1, a2, b1, b2, phi, t1, t2, eps)
  v_mds2 <- -ds2(a1, a2, b1, b2, phi, t1, t2, eps)
  v_s12  <- s12(a1, a2, b1, b2, phi, t1, t2, eps)
  
  #  Needs non-negative terms; allow tiny negatives from
  # floating-point error and clip with safe_log.
  tol_neg <- -sqrt(eps)
  if (any(!is.finite(v_ds12)) || any(!is.finite(v_mds1)) ||
      any(!is.finite(v_mds2)) || any(!is.finite(v_s12)) ||
      any(v_ds12 < tol_neg) || any(v_mds1 < tol_neg) ||
      any(v_mds2 < tol_neg) || any(v_s12 < tol_neg)) {
    return(1e20)
  }
  
  ll <- sum(cens1 * cens2 * safe_log(v_ds12, eps)) +
    sum(cens1 * (1 - cens2) * safe_log(v_mds1, eps)) +
    sum(cens2 * (1 - cens1) * safe_log(v_mds2, eps)) +
    sum((1 - cens1) * (1 - cens2) * safe_log(v_s12, eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
  
}

### independent likelihood
lik_dg_mv <- function(par, data, cens, linear_pred_a, 
                      linear_pred_b, eps = 1e-12) {
  # linear_pred_a for a
  # linear_pred_b for b 
  
  Xa <- stats::model.matrix(linear_pred_a, data) 
  Xb <- stats::model.matrix(linear_pred_b, data) 
  
  j_a <- ncol(Xa) 
  j_b <- ncol(Xb)
  
  npar <- j_a + j_b
  if (length(par) != npar) return(1e20)
  
  idx <- 0
  take_par <- function(k) {
    if (k == 0) return(numeric(0))
    out <- par[(idx + 1):(idx + k)]
    idx <<- idx + k
    out
  }
  
  reg_a <- take_par(j_a)
  reg_b <- take_par(j_b)  
  
  cens <- as.numeric(data[[cens]]) 
  
  t <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a, data)))
  
  a <- as.numeric(Xa %*% reg_a)
  b <- exp(as.numeric(Xb %*% reg_b)) 
  
  f_fail <- f_g(a, b, t, eps) 
  s_cens <- surv_g(a, b, t, eps)
  
  ll <- sum(cens*safe_log(f_fail, eps)) + 
    sum((1 - cens) * safe_log(s_cens, eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}



#function to calculate standard error and confidence interval
variances <- function(lik, est, ...) {
  hes_numDeriv <- numDeriv::hessian(func = lik, x = est, ...)
  hes_numDeriv <- 0.5 * (hes_numDeriv + t(hes_numDeriv))
  aux_numDeriv <- MASS::ginv(hes_numDeriv)
  var_numDeriv <- diag(aux_numDeriv)
  cont_neg <- sapply(var_numDeriv, function(x) ifelse(x < 0, 1, 0))
  if(any(var_numDeriv < 0)) {
    var_numDeriv <- abs(var_numDeriv)
    diag(aux_numDeriv) <- var_numDeriv
  }
  li <- est - 1.96*sqrt(var_numDeriv)
  ls <- est + 1.96*sqrt(var_numDeriv)
  out <- NULL
  out$covar <- aux_numDeriv
  out$ic <- data.frame(var = var_numDeriv, li = li, ls = ls, cont_neg)
  return(out)
}

# Univariate Delta Method
delta_uni <- function(f, est_mv, vari) {
  # f is the function
  # est_mv is the MLE of the argument of f
  # vari is the variance of the parameter which est_mv is the MLE
  var_delta <- vari*(numDeriv::grad(f, est_mv))^2
  li <- f(est_mv) - 1.96*sqrt(var_delta)
  ls <- f(est_mv) + 1.96*sqrt(var_delta)
  out <- data.frame(var = var_delta, li = li, ls = ls)
  return(out)
}

# Multivariate Delta Method
delta_multi <- function(f, est_mv, vari) {
  # f is the function
  # est_mv is the MLE of the argument of f
  # vari is the variance matrix of the parameter which est_mv is the MLE
  var_delta <- t(numDeriv::grad(f, est_mv))%*%vari%*%numDeriv::grad(f, est_mv)
  li <- f(est_mv) - 1.96*sqrt(var_delta)
  ls <- f(est_mv) + 1.96*sqrt(var_delta)
  out <- data.frame(var = var_delta, li = li, ls = ls)
  return(out)
}



ajuste <- function(linear_pred_a1, linear_pred_b1, 
                   linear_pred_a2, linear_pred_b2, 
                   cens1, cens2, dados) {
  
  j_a1 <- dim(stats::model.matrix(linear_pred_a1, dados))[2]
  j_a2 <- dim(stats::model.matrix(linear_pred_a2, dados))[2]
  j_b1 <- dim(stats::model.matrix(linear_pred_b1, dados))[2]
  j_b2 <- dim(stats::model.matrix(linear_pred_b2, dados))[2] 
  
  # nomes dos parâmetros
  nm_a1 <- paste0("a1_", colnames(stats::model.matrix(linear_pred_a1, dados)))
  nm_b1 <- paste0("b1_", colnames(stats::model.matrix(linear_pred_b1, dados))) 
  nm_a2 <- paste0("a2_", colnames(stats::model.matrix(linear_pred_a2, dados)))
  nm_b2 <- paste0("b2_", colnames(stats::model.matrix(linear_pred_b2, dados))) 
  nm_par <- c(nm_a1, nm_b1, nm_a2, nm_b2, "eta_phi")
  
  #### Model fit
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  theta0 <- c(rep(0, j_a1), rep(-6, j_b1), rep(0, j_a2), rep(-6, j_b2), 0)
  stopifnot(length(theta0) == npar)
  
  fit_nm <- optim(
    par = theta0, fn = lik_dg_fgm_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )
  
  fit <- try(
    optim(par = fit_nm$par, fn = lik_dg_fgm_mv, data = dados, 
          cens1 = cens1, cens2 = cens2, linear_pred_a1 = linear_pred_a1, 
          linear_pred_a2 = linear_pred_a2, linear_pred_b1 = linear_pred_b1, 
          linear_pred_b2 = linear_pred_b2, eps = 1e-12, 
          method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE)
  
  if(!inherits(fit, "try-error")){
    var_aux <- try(
      variances(
        lik = lik_dg_fgm_mv,
        est = fit$par,
        data = dados,
        cens1 = cens1,
        cens2 = cens2,
        linear_pred_a1 = linear_pred_a1,
        linear_pred_a2 = linear_pred_a2,
        linear_pred_b1 = linear_pred_b1,
        linear_pred_b2 = linear_pred_b2,
        eps = 1e-12
      ),
      silent = TRUE
    )
  }
  
  estimativas <- data.frame(
    variavel = c(rep("t1", j_a1 + j_b1), rep("t2", j_a2 + j_b2), "phi_reta"),
    parametro = nm_par,
    emv = fit$par,
    var = var_aux$ic$var,
    li = var_aux$ic$li,
    ls = var_aux$ic$ls 
  )
  
  l_mv <- length(fit$par)
  
  ## var e IC para phi #############
  efe_phi <- function(phi) tanh(phi)
  l_mv <- length(fit$par)
  out_phi <- delta_uni(f = efe_phi, est_mv = fit$par[l_mv], vari = var_aux$ic$var[l_mv])
  
  estimativas <- rbind(
    estimativas,
    data.frame(
      variavel = "dep",
      parametro = "phi",
      emv = efe_phi(fit$par[l_mv]),
      var = out_phi$var,
      li = out_phi$li,
      ls = out_phi$ls
    )
  )
  
  eta_hat <- estimativas$emv[estimativas$parametro == "eta_phi"]
  var_eta <- estimativas$var[estimativas$parametro == "eta_phi"]
  z <- qnorm(1 - 0.05/2)
  # phi estimado
  phi_hat <- tanh(eta_hat)
  
  # limites na escala de eta
  li_eta <- eta_hat - z * sqrt(var_eta)
  ls_eta <- eta_hat + z * sqrt(var_eta)
  
  # volta para escala de phi
  li_phi <- tanh(li_eta)
  ls_phi <- tanh(ls_eta)
  
  estimativas <- rbind(
    estimativas, 
    data.frame(
      variavel = "dep", 
      parametro = "phi_transf", 
      emv = phi_hat, 
      var = out_phi$var, 
      li = li_phi, 
      ls = ls_phi
    ) 
  )
  
  ## Cálculo e IC de p
  # Índices dos parâmetros no vetor total
  idx_a1 <- 1:j_a1
  idx_b1 <- (max(idx_a1) + 1):(max(idx_a1) + j_b1)
  idx_a2 <- (max(idx_b1) + 1):(max(idx_b1) + j_a2)
  idx_b2 <- (max(idx_a2) + 1):(max(idx_a2) + j_b2)
  idx_eta_phi <- max(idx_b2) + 1
  
  theta_hat <- fit$par            
  V_hat <- var_aux$covar
  
  # covariáveis usadas em cada submodelo (sem resposta)
  vars_a1 <- all.vars(delete.response(terms(linear_pred_a1)))
  vars_b1 <- all.vars(delete.response(terms(linear_pred_b1)))
  vars_a2 <- all.vars(delete.response(terms(linear_pred_a2)))
  vars_b2 <- all.vars(delete.response(terms(linear_pred_b2)))
  
  # união das covariáveis -> cenário geral
  vars_all <- unique(c(vars_a1, vars_b1, vars_a2, vars_b2))
  
  # combinações distintas de covariáveis no banco
  newdata <- unique(dados[, vars_all, drop = FALSE])
  
  # matrizes de desenho alinhadas na MESMA ordem de linhas
  X_a1 <- stats::model.matrix(delete.response(terms(linear_pred_a1)), newdata)
  X_b1 <- stats::model.matrix(delete.response(terms(linear_pred_b1)), newdata)
  X_a2 <- stats::model.matrix(delete.response(terms(linear_pred_a2)), newdata)
  X_b2 <- stats::model.matrix(delete.response(terms(linear_pred_b2)), newdata)
  
  eps_cure <- 1e-5
  clip01 <- function(x) pmin(pmax(x, 0), 1)
  
  out <- vector("list", nrow(newdata))
  
  for (i in seq_len(nrow(newdata))) {
    xa1 <- X_a1[i, , drop = FALSE]
    xb1 <- X_b1[i, , drop = FALSE]
    xa2 <- X_a2[i, , drop = FALSE]
    xb2 <- X_b2[i, , drop = FALSE]
    
    f_p1 <- function(pars) {
      a1 <- as.numeric(xa1 %*% pars[idx_a1])
      b1 <- exp(as.numeric(xb1 %*% pars[idx_b1])) 
      ifelse(a1 < -eps_cure, exp(b1 / a1), 0)
    }
    
    f_p2 <- function(pars) {
      a2 <- as.numeric(xa2 %*% pars[idx_a2])
      b2 <- exp(as.numeric(xb2 %*% pars[idx_b2])) 
      ifelse(a2 < -eps_cure, exp(b2 / a2), 0)
    }
    
    f_p12 <- function(pars) {
      p1 <- f_p1(pars)
      p2 <- f_p2(pars)
      phi <- tanh(pars[idx_eta_phi])
      p1 * p2 * (1 + phi * (1 - p1) * (1 - p2))
    }
    
    d1  <- delta_multi(f_p1,  theta_hat, V_hat)
    d2  <- delta_multi(f_p2,  theta_hat, V_hat)
    d12 <- delta_multi(f_p12, theta_hat, V_hat)
    
    out[[i]] <- cbind(
      newdata[i, , drop = FALSE],
      data.frame(
        p1_emv  = f_p1(theta_hat),  p1_var  = as.numeric(d1$var),  p1_li  = clip01(as.numeric(d1$li)),  p1_ls  = clip01(as.numeric(d1$ls)),
        p2_emv  = f_p2(theta_hat),  p2_var  = as.numeric(d2$var),  p2_li  = clip01(as.numeric(d2$li)),  p2_ls  = clip01(as.numeric(d2$ls)),
        p12_emv = f_p12(theta_hat), p12_var = as.numeric(d12$var), p12_li = clip01(as.numeric(d12$li)), p12_ls = clip01(as.numeric(d12$ls))
      )
    )
  }
  
  fracoes_cura_delta <- do.call(rbind, out)
  
  
  ##Gráficos
  reg_a1 <- fit$par[idx_a1]
  reg_b1 <-  fit$par[idx_b1] 
  reg_a2 <- fit$par[idx_a2]
  reg_b2 <- fit$par[idx_b2] 
  
  make_plot <- function(form_a, form_b, reg_a, reg_b, time_var, event_var, dados) {
    # união das covariáveis de a e b
    vars_a <- all.vars(delete.response(terms(form_a)))
    vars_b <- all.vars(delete.response(terms(form_b)))
    vars_all <- unique(c(vars_a, vars_b))
    
    # combinações observadas no banco  
    newdata <- unique(dados[, vars_all, drop = FALSE])
    
    # matrizes de desenho separadas para a e b, ambas no MESMO newdata
    Xa <- model.matrix(delete.response(terms(form_a)), newdata)
    Xb <- model.matrix(delete.response(terms(form_b)), newdata) 
    # rótulo do grupo com união das variáveis
    rotulo <- apply(newdata[, vars_all, drop = FALSE], 1, function(z) {
      paste(paste(vars_all, z, sep = "="), collapse = " | ")
    }) %>%  as.vector()
    
    # grade de tempo
    t_grid <- seq(0, max(dados[[time_var]], na.rm = TRUE), length.out = 300)
    
    # curva do modelo por grupo
    S_estimates <- do.call(rbind, lapply(seq_len(nrow(newdata)), function(i) {
      a_i <- as.numeric(Xa[i, , drop = FALSE] %*% reg_a)
      b_i <- exp(as.numeric(Xb[i, , drop = FALSE] %*% reg_b)) 
      surv_i <- vapply(t_grid, function(tt) surv_g(a = a_i, b = b_i, t = tt) * 100, numeric(1))
      data.frame(time = t_grid, surv = surv_i, grupo = rotulo[i], tipo = "Modelo")
    }))
    
    # KM com o mesmo agrupamento
    dados$grupo_plot <- apply(dados[, vars_all, drop = FALSE], 1, function(z) {
      paste(paste(vars_all, z, sep = "="), collapse = " | ")
    })
    
    sf <- survival::survfit(
      as.formula(paste0("survival::Surv(", time_var, ",", event_var, ") ~ grupo_plot")),
      data = dados
    )
    
    km_levels <- sub("^grupo_plot=", "", names(sf$strata))
    
    ekm <- ggsurvplot(
      sf, data = dados,
      ylab = "Survival probability",
      xlab = "Time",
      linetype = "solid",
      palette = "lancet",
      fun = "pct",
      censor.shape = " ",
      size = 1.1,
      legend.title = "",
      legend.labs = km_levels,
      legend = "bottom"
    )
    
    # KM: remove prefixo e fixa níveis
    ekm$plot$data$strata <- factor(
      sub("^grupo_plot=", "", as.character(ekm$plot$data$strata)),
      levels = km_levels
    )
    
    # Modelo: usa a MESMA variável e MESMOS níveis
    S_estimates$strata <- factor(S_estimates$grupo, levels = km_levels)
    
    # Overlay (sem criar nova escala)
    ekm_group <- ekm$plot +
      geom_line(
        data = S_estimates,
        aes(x = time, y = surv, colour = strata, group = strata),
        linetype = "dashed",
        linewidth = 1.1
      ) +
      theme_minimal() +
      theme(legend.position = "bottom") 
    # +
    #   coord_cartesian(xlim = c(0, 40)) 
    
    ekm_group
  }
  
  plot_t1 <- make_plot(
    form_a = linear_pred_a1, form_b = linear_pred_b1,
    reg_a = reg_a1, reg_b = reg_b1,
    time_var = all.vars(linear_pred_a1)[1],
    event_var = cens1,
    dados = dados
  )
  
  plot_t2 <- make_plot(
    form_a = linear_pred_a2, form_b = linear_pred_b2,
    reg_a = reg_a2, reg_b = reg_b2,
    time_var = all.vars(linear_pred_a2)[1],
    event_var = cens2,
    dados = dados
  )
  
  out_final <- list(estimates = estimativas, 
                    fracoes_cura = fracoes_cura_delta,
                    plot_t1 = plot_t1,
                    plot_t2 = plot_t2)
  return(out_final)
  
}




ajuste_covar_cont <- function(linear_pred_a1, linear_pred_b1, 
                              linear_pred_a2, linear_pred_b2, 
                              cens1, cens2, dados) {
  
  j_a1 <- dim(stats::model.matrix(linear_pred_a1, dados))[2]
  j_a2 <- dim(stats::model.matrix(linear_pred_a2, dados))[2]
  j_b1 <- dim(stats::model.matrix(linear_pred_b1, dados))[2]
  j_b2 <- dim(stats::model.matrix(linear_pred_b2, dados))[2] 
  
  # nomes dos parâmetros
  nm_a1 <- paste0("a1_", colnames(stats::model.matrix(linear_pred_a1, dados)))
  nm_b1 <- paste0("b1_", colnames(stats::model.matrix(linear_pred_b1, dados))) 
  nm_a2 <- paste0("a2_", colnames(stats::model.matrix(linear_pred_a2, dados)))
  nm_b2 <- paste0("b2_", colnames(stats::model.matrix(linear_pred_b2, dados))) 
  nm_par <- c(nm_a1, nm_b1, nm_a2, nm_b2, "eta_phi")
  
  #### Model fit
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  theta0 <- c(rep(0, j_a1), rep(0, j_b1), rep(0, j_a2), rep(0, j_b2), 0)
  stopifnot(length(theta0) == npar)
  
  fit_nm <- optim(
    par = theta0, fn = lik_dg_fgm_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )
  
  fit <- try(
    optim(par = fit_nm$par, fn = lik_dg_fgm_mv, data = dados, 
          cens1 = cens1, cens2 = cens2, linear_pred_a1 = linear_pred_a1, 
          linear_pred_a2 = linear_pred_a2, linear_pred_b1 = linear_pred_b1, 
          linear_pred_b2 = linear_pred_b2, eps = 1e-12, 
          method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE)
  
  if(!inherits(fit, "try-error")){
    var_aux <- try(
      variances(
        lik = lik_dg_fgm_mv,
        est = fit$par,
        data = dados,
        cens1 = cens1,
        cens2 = cens2,
        linear_pred_a1 = linear_pred_a1,
        linear_pred_a2 = linear_pred_a2,
        linear_pred_b1 = linear_pred_b1,
        linear_pred_b2 = linear_pred_b2,
        eps = 1e-12
      ),
      silent = TRUE
    )
  }
  
  estimativas <- data.frame(
    variavel = c(rep("t1", j_a1 + j_b1), rep("t2", j_a2 + j_b2), "phi_reta"),
    parametro = nm_par,
    emv = fit$par,
    var = var_aux$ic$var,
    li = var_aux$ic$li,
    ls = var_aux$ic$ls 
  )
  
  l_mv <- length(fit$par)
  
  ## var e IC para phi #############
  efe_phi <- function(phi) tanh(phi)
  l_mv <- length(fit$par)
  out_phi <- delta_uni(f = efe_phi, est_mv = fit$par[l_mv], vari = var_aux$ic$var[l_mv])
  
  estimativas <- rbind(
    estimativas,
    data.frame(
      variavel = "dep",
      parametro = "phi",
      emv = efe_phi(fit$par[l_mv]),
      var = out_phi$var,
      li = out_phi$li,
      ls = out_phi$ls
    )
  )
  
  eta_hat <- estimativas$emv[estimativas$parametro == "eta_phi"]
  var_eta <- estimativas$var[estimativas$parametro == "eta_phi"]
  z <- qnorm(1 - 0.05/2)
  # phi estimado
  phi_hat <- tanh(eta_hat)
  
  # limites na escala de eta
  li_eta <- eta_hat - z * sqrt(var_eta)
  ls_eta <- eta_hat + z * sqrt(var_eta)
  
  # volta para escala de phi
  li_phi <- tanh(li_eta)
  ls_phi <- tanh(ls_eta)
  
  estimativas <- rbind(
    estimativas, 
    data.frame(
      variavel = "dep", 
      parametro = "phi_transf", 
      emv = phi_hat, 
      var = out_phi$var, 
      li = li_phi, 
      ls = ls_phi
    ) 
  )
  out_final <- NULL
  out_final$estimates <- estimativas
  out_final$var_covar <- var_aux$covar
  return(out_final)
  
}

fracoes_covar_cont <- function(linear_pred_a1, linear_pred_b1,
                                       linear_pred_a2, linear_pred_b2,
                                       cens1, cens2, dados,
                                       var_cont,   # nome da variável contínua (string)
                                       var_cat,    # nome da variável categórica (string)
                                       fit_out,    # saída de ajuste_covar_cont (estimativas)
                                       fit_par,    # fit$par (vetor de parâmetros)
                                       var_aux_covar # var_aux$covar (matriz de covariâncias)
                                       ) {
  
  # ---------- índices ----------
  j_a1 <- dim(stats::model.matrix(linear_pred_a1, dados))[2]
  j_b1 <- dim(stats::model.matrix(linear_pred_b1, dados))[2]
  j_a2 <- dim(stats::model.matrix(linear_pred_a2, dados))[2]
  j_b2 <- dim(stats::model.matrix(linear_pred_b2, dados))[2]
  
  idx_a1      <- 1:j_a1
  idx_b1      <- (max(idx_a1) + 1):(max(idx_a1) + j_b1)
  idx_a2      <- (max(idx_b1) + 1):(max(idx_b1) + j_a2)
  idx_b2      <- (max(idx_a2) + 1):(max(idx_a2) + j_b2)
  idx_eta_phi <- max(idx_b2) + 1
  
  theta_hat <- fit_par
  V_hat     <- var_aux_covar
  
  # ---------- grid de predição ----------
  q_cont <- quantile(dados[[var_cont]], probs = c(0.25, 0.5, 0.75))
  cats   <- sort(unique(dados[[var_cat]]))
  
  newdata <- expand.grid(
    cat  = cats,
    cont = q_cont
  )
  names(newdata) <- c(var_cat, var_cont)
  newdata$quartil <- rep(c("Q1", "Mediana", "Q3"), each = length(cats))
  
  # matrizes de desenho no newdata
  X_a1 <- stats::model.matrix(delete.response(terms(linear_pred_a1)), newdata)
  X_b1 <- stats::model.matrix(delete.response(terms(linear_pred_b1)), newdata)
  X_a2 <- stats::model.matrix(delete.response(terms(linear_pred_a2)), newdata)
  X_b2 <- stats::model.matrix(delete.response(terms(linear_pred_b2)), newdata)
  
  eps_cure <- 1e-5
  clip01   <- function(x) pmin(pmax(x, 0), 1)
  
  # ---------- frações de cura ----------
  out <- vector("list", nrow(newdata))
  
  for (i in seq_len(nrow(newdata))) {
    xa1 <- X_a1[i, , drop = FALSE]
    xb1 <- X_b1[i, , drop = FALSE]
    xa2 <- X_a2[i, , drop = FALSE]
    xb2 <- X_b2[i, , drop = FALSE]
    
    f_p1 <- function(pars) {
      a1 <- as.numeric(xa1 %*% pars[idx_a1])
      b1 <- exp(as.numeric(xb1 %*% pars[idx_b1]))
      ifelse(a1 < -eps_cure, exp(b1 / a1), 0)
    }
    
    f_p2 <- function(pars) {
      a2 <- as.numeric(xa2 %*% pars[idx_a2])
      b2 <- exp(as.numeric(xb2 %*% pars[idx_b2]))
      ifelse(a2 < -eps_cure, exp(b2 / a2), 0)
    }
    
    f_p12 <- function(pars) {
      p1  <- f_p1(pars)
      p2  <- f_p2(pars)
      phi <- tanh(pars[idx_eta_phi])
      p1 * p2 * (1 + phi * (1 - p1) * (1 - p2))
    }
    
    d1  <- delta_multi(f_p1,  theta_hat, V_hat)
    d2  <- delta_multi(f_p2,  theta_hat, V_hat)
    d12 <- delta_multi(f_p12, theta_hat, V_hat)
    
    out[[i]] <- data.frame(
      newdata[i, , drop = FALSE],
      p1_emv  = f_p1(theta_hat),
      p1_se   = sqrt(as.numeric(d1$var)),
      p1_li   = clip01(as.numeric(d1$li)),
      p1_ls   = clip01(as.numeric(d1$ls)),
      p2_emv  = f_p2(theta_hat),
      p2_se   = sqrt(as.numeric(d2$var)),
      p2_li   = clip01(as.numeric(d2$li)),
      p2_ls   = clip01(as.numeric(d2$ls)),
      p12_emv = f_p12(theta_hat),
      p12_se  = sqrt(as.numeric(d12$var)),
      p12_li  = clip01(as.numeric(d12$li)),
      p12_ls  = clip01(as.numeric(d12$ls))
    )
  }
  
  fracoes_cura <- do.call(rbind, out)
  

  list(
    fracoes_cura = fracoes_cura
  )
}


ajuste_covar_cont_phi0 <- function(linear_pred_a1, linear_pred_b1, 
                                   linear_pred_a2, linear_pred_b2, 
                                   cens1, cens2, dados) {
  
  j_a1 <- dim(stats::model.matrix(linear_pred_a1, dados))[2]
  j_a2 <- dim(stats::model.matrix(linear_pred_a2, dados))[2]
  j_b1 <- dim(stats::model.matrix(linear_pred_b1, dados))[2]
  j_b2 <- dim(stats::model.matrix(linear_pred_b2, dados))[2] 
  
  # nomes dos parâmetros
  nm_a1 <- paste0("a1_", colnames(stats::model.matrix(linear_pred_a1, dados)))
  nm_b1 <- paste0("b1_", colnames(stats::model.matrix(linear_pred_b1, dados))) 
  nm_a2 <- paste0("a2_", colnames(stats::model.matrix(linear_pred_a2, dados)))
  nm_b2 <- paste0("b2_", colnames(stats::model.matrix(linear_pred_b2, dados))) 
  nm_par <- c(nm_a1, nm_b1, nm_a2, nm_b2)
  
  #### Model fit
  npar <- j_a1 + j_b1 + j_a2 + j_b2 
  theta0 <- c(rep(0, j_a1), rep(0, j_b1), rep(0, j_a2), rep(0, j_b2))
  stopifnot(length(theta0) == npar)
  
  fit_nm <- optim(
    par = theta0, fn = lik_dg_phi0_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )
  
  fit <- try(
    optim(par = fit_nm$par, fn = lik_dg_phi0_mv, data = dados, 
          cens1 = cens1, cens2 = cens2, linear_pred_a1 = linear_pred_a1, 
          linear_pred_a2 = linear_pred_a2, linear_pred_b1 = linear_pred_b1, 
          linear_pred_b2 = linear_pred_b2, eps = 1e-12, 
          method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE)
  
  if(!inherits(fit, "try-error")){
    var_aux <- try(
      variances(
        lik = lik_dg_phi0_mv,
        est = fit$par,
        data = dados,
        cens1 = cens1,
        cens2 = cens2,
        linear_pred_a1 = linear_pred_a1,
        linear_pred_a2 = linear_pred_a2,
        linear_pred_b1 = linear_pred_b1,
        linear_pred_b2 = linear_pred_b2,
        eps = 1e-12
      ),
      silent = TRUE
    )

  }
  
  estimativas <- data.frame(
    variavel = c(rep("t1", j_a1 + j_b1), rep("t2", j_a2 + j_b2)),
    parametro = nm_par,
    emv = fit$par,
    var = var_aux$ic$var,
    li = var_aux$ic$li,
    ls = var_aux$ic$ls
  )
  
  out_final <- NULL
  out_final$estimates <- estimativas
  out_final$var_covar <- var_aux$covar
  return(out_final)
  
}

AIC_calc <- function(loglik, k) {
  aic <- -2 * loglik + 2 * k
  return(aic)
}

AICc_calc <- function(loglik, k, n) {
  aic <- AIC_calc(loglik = loglik, k = k)
  den <- n - k - 1
  if (den <= 0) return(Inf)
  aic + (2 * k * (k + 1)) / den
}

######################## Frank copula ############################
############################
##    Frank survival      ##
##    copula pieces       ##
############################

## --- Standard Frank copula C_F(u, v; theta), u,v in [0,1] ---
frank_cdf_copula <- function(u, v, theta, eps = 1e-12) {
  eps_th <- 1e-8
  ifelse(
    abs(theta) < eps_th,
    u * v,
    {
      num <- (exp(-theta * u) - 1) * (exp(-theta * v) - 1)
      den <- exp(-theta) - 1
      arg <- pmax(1 + num / den, eps)
      -(1 / theta) * log(arg)
    }
  )
}

## --- dC_F/du ---
frank_dcdu <- function(u, v, theta, eps = 1e-12) {
  eps_th <- 1e-8
  ifelse(
    abs(theta) < eps_th,
    v,
    {
      A <- exp(-theta * u) - 1
      B <- exp(-theta * v) - 1
      D <- exp(-theta) - 1
      exp(-theta * u) * B / pmax(D + A * B, eps)
    }
  )
}

## --- Copula density c_F(u,v) = d2C_F/(du dv) ---
frank_density <- function(u, v, theta, eps = 1e-12) {
  eps_th <- 1e-8
  ifelse(
    abs(theta) < eps_th,
    1,
    {
      A <- exp(-theta * u) - 1
      B <- exp(-theta * v) - 1
      D <- exp(-theta) - 1
      -theta * D * exp(-theta * (u + v)) / pmax((D + A * B)^2, eps)
    }
  )
}

## -------------------------------------------------------------------
## SURVIVAL COPULA of Frank:
##   Chat_F(s1,s2) = s1 + s2 - 1 + C_F(1-s1, 1-s2; theta)
##
## Used as S(t1,t2) = Chat_F(S1(t1), S2(t2)).
## -------------------------------------------------------------------

frank_surv_copula <- function(s1, s2, theta, eps = 1e-12) {
  s1 + s2 - 1 + frank_cdf_copula(1 - s1, 1 - s2, theta, eps)
}

## dChat/ds1 = 1 - dC_F/du |_{u=1-s1, v=1-s2}
frank_surv_ds1 <- function(s1, s2, theta, eps = 1e-12) {
  1 - frank_dcdu(1 - s1, 1 - s2, theta, eps)
}

## dChat/ds2 = 1 - dC_F/dv |_{u=1-s1, v=1-s2}
## By symmetry of Frank: dC_F/dv(u,v) = dC_F/du(v,u)
frank_surv_ds2 <- function(s1, s2, theta, eps = 1e-12) {
  1 - frank_dcdu(1 - s2, 1 - s1, theta, eps)
}

## d2Chat/(ds1 ds2) = c_F(1-s1, 1-s2; theta)
frank_surv_d2 <- function(s1, s2, theta, eps = 1e-12) {
  frank_density(1 - s1, 1 - s2, theta, eps)
}

############################
## Joint survival and  ##
##    likelihood pieces   ##
############################

s12_frank <- function(a1, a2, b1, b2, theta, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  frank_surv_copula(s1, s2, theta, eps)
}

## dS/dt1 = dChat/ds1 * (-f1)
ds1_frank <- function(a1, a2, b1, b2, theta, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f1 <- f_g(a1, b1, t1, eps)
  frank_surv_ds1(s1, s2, theta, eps) * (-f1)
}

## dS/dt2 = dChat/ds2 * (-f2)
ds2_frank <- function(a1, a2, b1, b2, theta, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f2 <- f_g(a2, b2, t2, eps)
  frank_surv_ds2(s1, s2, theta, eps) * (-f2)
}

## d2S/(dt1 dt2) = c_F(1-s1, 1-s2) * f1 * f2
ds12_frank <- function(a1, a2, b1, b2, theta, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f1 <- f_g(a1, b1, t1, eps)
  f2 <- f_g(a2, b2, t2, eps)
  frank_surv_d2(s1, s2, theta, eps) * f1 * f2
}

## Joint cure fraction
p12_frank <- function(p1, p2, theta, eps = 1e-12) {
  frank_surv_copula(p1, p2, theta, eps)
}

############################
##  Kendall's tau       ##
############################

debye1 <- function(theta, eps_th = 1e-8) {
  if (abs(theta) < eps_th) return(1)
  integrate(function(t) t / (exp(t) - 1),
            lower = 0, upper = theta, rel.tol = 1e-10)$value / theta
}

tau_frank <- function(theta, eps_th = 1e-8) {
  if (abs(theta) < eps_th) return(0)
  1 - 4 / theta * (1 - debye1(theta))
}

############################
## Likelihood          ##
############################

lik_dg_frank_mv <- function(par, data, cens1, cens2,
                            linear_pred_a1, linear_pred_a2,
                            linear_pred_b1, linear_pred_b2,
                            eps = 1e-12) {
  Xa1 <- stats::model.matrix(linear_pred_a1, data)
  Xa2 <- stats::model.matrix(linear_pred_a2, data)
  Xb1 <- stats::model.matrix(linear_pred_b1, data)
  Xb2 <- stats::model.matrix(linear_pred_b2, data)
  
  j_a1 <- ncol(Xa1); j_a2 <- ncol(Xa2)
  j_b1 <- ncol(Xb1); j_b2 <- ncol(Xb2)
  
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  if (length(par) != npar) return(1e20)
  
  idx <- 0
  take_par <- function(k) {
    if (k == 0) return(numeric(0))
    out <- par[(idx + 1):(idx + k)]; idx <<- idx + k; out
  }
  
  reg_a1 <- take_par(j_a1); reg_b1 <- take_par(j_b1)
  reg_a2 <- take_par(j_a2); reg_b2 <- take_par(j_b2)
  theta  <- take_par(1)
  
  d1 <- as.numeric(data[[cens1]])
  d2 <- as.numeric(data[[cens2]])
  t1 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a1, data)))
  t2 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a2, data)))
  
  a1 <- as.numeric(Xa1 %*% reg_a1)
  a2 <- as.numeric(Xa2 %*% reg_a2)
  b1 <- exp(as.numeric(Xb1 %*% reg_b1))
  b2 <- exp(as.numeric(Xb2 %*% reg_b2))
  
  v_ds12 <-  ds12_frank(a1, a2, b1, b2, theta, t1, t2, eps)
  v_mds1 <- -ds1_frank (a1, a2, b1, b2, theta, t1, t2, eps)
  v_mds2 <- -ds2_frank (a1, a2, b1, b2, theta, t1, t2, eps)
  v_s12  <-  s12_frank (a1, a2, b1, b2, theta, t1, t2, eps)
  
  tol_neg <- -sqrt(eps)
  if (any(!is.finite(v_ds12)) || any(!is.finite(v_mds1)) ||
      any(!is.finite(v_mds2)) || any(!is.finite(v_s12))  ||
      any(v_ds12 < tol_neg)   || any(v_mds1 < tol_neg)  ||
      any(v_mds2 < tol_neg)   || any(v_s12  < tol_neg)) return(1e20)
  
  ll <- sum( d1      *  d2      * safe_log(v_ds12, eps)) +
    sum( d1      * (1 - d2) * safe_log(v_mds1, eps)) +
    sum((1 - d1) *  d2      * safe_log(v_mds2, eps)) +
    sum((1 - d1) * (1 - d2) * safe_log(v_s12,  eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}

####################################
## ajuste_frank_covar_cont()   ##
####################################

ajuste_frank_covar_cont <- function(linear_pred_a1, linear_pred_b1,
                                    linear_pred_a2, linear_pred_b2,
                                    cens1, cens2, dados) {
  
  j_a1 <- ncol(stats::model.matrix(linear_pred_a1, dados))
  j_a2 <- ncol(stats::model.matrix(linear_pred_a2, dados))
  j_b1 <- ncol(stats::model.matrix(linear_pred_b1, dados))
  j_b2 <- ncol(stats::model.matrix(linear_pred_b2, dados))
  
  nm_a1  <- paste0("a1_", colnames(stats::model.matrix(linear_pred_a1, dados)))
  nm_b1  <- paste0("b1_", colnames(stats::model.matrix(linear_pred_b1, dados)))
  nm_a2  <- paste0("a2_", colnames(stats::model.matrix(linear_pred_a2, dados)))
  nm_b2  <- paste0("b2_", colnames(stats::model.matrix(linear_pred_b2, dados)))
  nm_par <- c(nm_a1, nm_b1, nm_a2, nm_b2, "theta_frank")
  
  npar   <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  theta0 <- c(rep(0, j_a1), rep(0, j_b1), rep(0, j_a2), rep(0, j_b2), -1)
  stopifnot(length(theta0) == npar)
  
  fit_nm <- optim(
    par = theta0, fn = lik_dg_frank_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10))
  
  fit <- try(optim(
    par = fit_nm$par, fn = lik_dg_frank_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "BFGS", hessian = TRUE,
    control = list(maxit = 2000, reltol = 1e-10)), silent = TRUE)
  
  var_aux <- try(variances(
    lik = lik_dg_frank_mv, est = fit$par,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12), silent = TRUE)
  
  idx_a1    <- 1:j_a1
  idx_b1    <- (max(idx_a1)+1):(max(idx_a1)+j_b1)
  idx_a2    <- (max(idx_b1)+1):(max(idx_b1)+j_a2)
  idx_b2    <- (max(idx_a2)+1):(max(idx_a2)+j_b2)
  idx_theta <- max(idx_b2)+1
  theta_hat <- fit$par
  V_hat     <- var_aux$covar
  
  estimativas <- data.frame(
    variavel  = c(rep("t1",j_a1+j_b1), rep("t2",j_a2+j_b2), "theta_frank"),
    parametro = nm_par,
    emv       = fit$par,
    var       = var_aux$ic$var,
    li        = var_aux$ic$li,
    ls        = var_aux$ic$ls)
  
  f_tau   <- function(pars) tau_frank(pars[idx_theta])
  out_tau <- delta_multi(f_tau, theta_hat, V_hat)
  estimativas <- rbind(estimativas, data.frame(
    variavel="dep", parametro="tau_kendall",
    emv=tau_frank(theta_hat[idx_theta]),
    var=as.numeric(out_tau$var),
    li=as.numeric(out_tau$li), ls=as.numeric(out_tau$ls)))
  
  list(estimates = estimativas, var_covar = var_aux$covar)
}


####################################################
#### Gaussian fit #################################
###################################################
gaussian_cdf_copula <- function(u, v, rho, eps = 1e-12) {
  u <- ifelse(is.finite(u), u, eps)
  v <- ifelse(is.finite(v), v, eps)
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)
  z1 <- stats::qnorm(u)
  z2 <- stats::qnorm(v)
  corr <- matrix(c(1, rho, rho, 1), nrow = 2)
  out <- mapply(function(x, y) {
    as.numeric(mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(x, y), sigma = corr))
  }, z1, z2)
  as.numeric(out)
}

gaussian_dcdu <- function(u, v, rho, eps = 1e-12) {
  u <- ifelse(is.finite(u), u, eps)
  v <- ifelse(is.finite(v), v, eps)
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)
  z1 <- stats::qnorm(u)
  z2 <- stats::qnorm(v)
  stats::pnorm((z2 - rho * z1) / sqrt(1 - rho^2))
}

gaussian_density <- function(u, v, rho, eps = 1e-12) {
  u <- ifelse(is.finite(u), u, eps)
  v <- ifelse(is.finite(v), v, eps)
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)
  z1 <- stats::qnorm(u)
  z2 <- stats::qnorm(v)
  (1 / sqrt(1 - rho^2)) * exp((2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) / (2 * (1 - rho^2)))
}

gaussian_surv_copula <- function(s1, s2, rho, eps = 1e-12) {
  s1 + s2 - 1 + gaussian_cdf_copula(1 - s1, 1 - s2, rho, eps)
}

gaussian_surv_ds1 <- function(s1, s2, rho, eps = 1e-12) {
  1 - gaussian_dcdu(1 - s1, 1 - s2, rho, eps)
}

gaussian_surv_ds2 <- function(s1, s2, rho, eps = 1e-12) {
  1 - gaussian_dcdu(1 - s2, 1 - s1, rho, eps)
}

gaussian_surv_d2 <- function(s1, s2, rho, eps = 1e-12) {
  gaussian_density(1 - s1, 1 - s2, rho, eps)
}

p12_gaussian <- function(p1, p2, rho, eps = 1e-12) {
  gaussian_surv_copula(p1, p2, rho, eps)
}

tau_gaussian <- function(rho) {
  (2 / pi) * asin(rho)
}

s12_gaussian <- function(a1, a2, b1, b2, rho, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  gaussian_surv_copula(s1, s2, rho, eps)
}

ds1_gaussian <- function(a1, a2, b1, b2, rho, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f1 <- f_g(a1, b1, t1, eps)
  gaussian_surv_ds1(s1, s2, rho, eps) * (-f1)
}

ds2_gaussian <- function(a1, a2, b1, b2, rho, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f2 <- f_g(a2, b2, t2, eps)
  gaussian_surv_ds2(s1, s2, rho, eps) * (-f2)
}

ds12_gaussian <- function(a1, a2, b1, b2, rho, t1, t2, eps = 1e-10) {
  s1 <- surv_g(a1, b1, t1, eps)
  s2 <- surv_g(a2, b2, t2, eps)
  f1 <- f_g(a1, b1, t1, eps)
  f2 <- f_g(a2, b2, t2, eps)
  gaussian_surv_d2(s1, s2, rho, eps) * f1 * f2
}

lik_dg_gaussian_mv <- function(par, data, cens1, cens2,
                               linear_pred_a1, linear_pred_a2,
                               linear_pred_b1, linear_pred_b2,
                               eps = 1e-12) {
  Xa1 <- stats::model.matrix(linear_pred_a1, data)
  Xa2 <- stats::model.matrix(linear_pred_a2, data)
  Xb1 <- stats::model.matrix(linear_pred_b1, data)
  Xb2 <- stats::model.matrix(linear_pred_b2, data)
  
  j_a1 <- ncol(Xa1); j_a2 <- ncol(Xa2)
  j_b1 <- ncol(Xb1); j_b2 <- ncol(Xb2)
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  if (length(par) != npar) return(1e20)
  
  idx <- 0
  take <- function(k) { out <- par[(idx + 1):(idx + k)]; idx <<- idx + k; out }
  
  reg_a1 <- take(j_a1); reg_b1 <- take(j_b1)
  reg_a2 <- take(j_a2); reg_b2 <- take(j_b2)
  rho <- tanh(take(1))
  
  d1 <- as.numeric(data[[cens1]])
  d2 <- as.numeric(data[[cens2]])
  t1 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a1, data)))
  t2 <- as.numeric(stats::model.response(stats::model.frame(linear_pred_a2, data)))
  
  a1 <- as.numeric(Xa1 %*% reg_a1)
  a2 <- as.numeric(Xa2 %*% reg_a2)
  b1 <- exp(as.numeric(Xb1 %*% reg_b1))
  b2 <- exp(as.numeric(Xb2 %*% reg_b2))
  
  v_ds12 <- ds12_gaussian(a1, a2, b1, b2, rho, t1, t2, eps)
  v_mds1 <- -ds1_gaussian(a1, a2, b1, b2, rho, t1, t2, eps)
  v_mds2 <- -ds2_gaussian(a1, a2, b1, b2, rho, t1, t2, eps)
  v_s12  <- s12_gaussian(a1, a2, b1, b2, rho, t1, t2, eps)
  
  tol_neg <- -sqrt(eps)
  if (any(!is.finite(v_ds12)) || any(!is.finite(v_mds1)) || any(!is.finite(v_mds2)) || any(!is.finite(v_s12)) ||
      any(v_ds12 < tol_neg) || any(v_mds1 < tol_neg) || any(v_mds2 < tol_neg) || any(v_s12 < tol_neg)) return(1e20)
  
  ll <- sum(d1 * d2 * safe_log(v_ds12, eps)) +
    sum(d1 * (1 - d2) * safe_log(v_mds1, eps)) +
    sum((1 - d1) * d2 * safe_log(v_mds2, eps)) +
    sum((1 - d1) * (1 - d2) * safe_log(v_s12, eps))
  
  out <- -ll
  if (!is.finite(out)) out <- 1e20
  out
}

ajuste_gaussian_covar_cont <- function(linear_pred_a1, linear_pred_b1,
                                       linear_pred_a2, linear_pred_b2,
                                       cens1, cens2, dados) {
  j_a1 <- ncol(stats::model.matrix(linear_pred_a1, dados))
  j_a2 <- ncol(stats::model.matrix(linear_pred_a2, dados))
  j_b1 <- ncol(stats::model.matrix(linear_pred_b1, dados))
  j_b2 <- ncol(stats::model.matrix(linear_pred_b2, dados))
  
  nm_a1 <- paste0("a1_", colnames(stats::model.matrix(linear_pred_a1, dados)))
  nm_b1 <- paste0("b1_", colnames(stats::model.matrix(linear_pred_b1, dados)))
  nm_a2 <- paste0("a2_", colnames(stats::model.matrix(linear_pred_a2, dados)))
  nm_b2 <- paste0("b2_", colnames(stats::model.matrix(linear_pred_b2, dados)))
  nm_par <- c(nm_a1, nm_b1, nm_a2, nm_b2, "eta_rho")
  
  npar <- j_a1 + j_b1 + j_a2 + j_b2 + 1
  theta0 <- c(rep(0, j_a1), rep(0, j_b1), rep(0, j_a2), rep(0, j_b2), 0)
  
  fit_nm <- optim(
    par = theta0, fn = lik_dg_gaussian_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )
  
  fit <- optim(
    par = fit_nm$par, fn = lik_dg_gaussian_mv,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12, method = "BFGS", hessian = TRUE,
    control = list(maxit = 2000, reltol = 1e-10)
  )
  
  var_aux <- variances(
    lik = lik_dg_gaussian_mv, est = fit$par,
    data = dados, cens1 = cens1, cens2 = cens2,
    linear_pred_a1 = linear_pred_a1, linear_pred_a2 = linear_pred_a2,
    linear_pred_b1 = linear_pred_b1, linear_pred_b2 = linear_pred_b2,
    eps = 1e-12
  )
  
  idx_a1 <- 1:j_a1
  idx_b1 <- (max(idx_a1) + 1):(max(idx_a1) + j_b1)
  idx_a2 <- (max(idx_b1) + 1):(max(idx_b1) + j_a2)
  idx_b2 <- (max(idx_a2) + 1):(max(idx_a2) + j_b2)
  idx_eta <- max(idx_b2) + 1
  
  theta_hat <- fit$par
  V_hat <- var_aux$covar
  
  estimativas <- data.frame(
    variavel = c(rep("t1", j_a1 + j_b1), rep("t2", j_a2 + j_b2), "rho_reta"),
    parametro = nm_par,
    emv = theta_hat,
    var = var_aux$ic$var,
    li = var_aux$ic$li,
    ls = var_aux$ic$ls 
  )
  
  f_rho <- function(pars) tanh(pars[idx_eta])
  out_rho <- delta_multi(f_rho, theta_hat, V_hat)
  estimativas <- rbind(estimativas, data.frame(
    variavel = "dep", parametro = "rho",
    emv = tanh(theta_hat[idx_eta]),
    var = as.numeric(out_rho$var),
    li = as.numeric(out_rho$li),
    ls = as.numeric(out_rho$ls) 
  ))
  
  out_tau <- delta_multi(function(pars) tau_gaussian(tanh(pars[idx_eta])), theta_hat, V_hat)
  estimativas <- rbind(estimativas, data.frame(
    variavel = "dep", parametro = "tau_kendall",
    emv = tau_gaussian(tanh(theta_hat[idx_eta])),
    var = as.numeric(out_tau$var),
    li = as.numeric(out_tau$li),
    ls = as.numeric(out_tau$ls) 
  ))
  
  
  
  list(
    estimates = estimativas,
    var_covar = var_aux$covar
  )
} 

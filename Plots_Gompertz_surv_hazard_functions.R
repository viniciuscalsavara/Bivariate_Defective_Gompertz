#===========================================================
# Hazard, Survival, Cumulative Hazard and Hazard Ratio
# under defective Gompertz regression model
#===========================================================

rm(list = ls())

library(ggplot2)

#-----------------------------------------------------------
# Functions
#-----------------------------------------------------------

# Hazard
haz_gomp <- function(t, x, beta0, beta1, alpha0, alpha1) {
  exp(beta0 + beta1*x + (alpha0 + alpha1*x)*t)
}

# Survival
surv_gomp <- function(t, x, beta0, beta1, alpha0, alpha1, eps = 1e-10) {
  a <- alpha0 + alpha1*x
  b <- exp(beta0 + beta1*x)
  
  if (abs(a) < eps) {
    return(exp(-b*t))
  } else {
    return(exp(-(b/a)*(exp(a*t)-1)))
  }
}

# Cumulative hazard
cumhaz_gomp <- function(t, x, beta0, beta1, alpha0, alpha1, eps = 1e-10) {
  a <- alpha0 + alpha1*x
  b <- exp(beta0 + beta1*x)
  
  if (abs(a) < eps) {
    return(b*t)
  } else {
    return((b/a)*(exp(a*t)-1))
  }
}

# Hazard ratio
hr_gomp <- function(t, beta1, alpha1) {
  exp(beta1 + alpha1*t)
}

#-----------------------------------------------------------
# Time grid
#-----------------------------------------------------------
t_grid <- seq(0, 20, length.out = 500)

#-----------------------------------------------------------
# Scenarios
#-----------------------------------------------------------
scenarios <- list(
  
  list(panel="Improper distribution", scenario="Pattern 1",
       beta0=-2.0, beta1=1, alpha0=-0.25, alpha1=-0.8),
  
  list(panel="Improper distribution", scenario="Pattern 2",
       beta0=-3.0, beta1=1.5, alpha0=-0.05, alpha1=-0.10),
  
  list(panel="Proper distribution", scenario="Pattern 3",
       beta0=-2, beta1=-0.6, alpha0=0.05, alpha1=0.05),
  
  list(panel="Proper distribution", scenario="Pattern 4",
       beta0=-1, beta1=0.5, alpha0=0, alpha1=0)
)

#-----------------------------------------------------------
# Build data
#-----------------------------------------------------------
build_data <- function(fun, yname) {
  
  do.call(
    rbind,
    lapply(scenarios, function(s) {
      
      d0 <- data.frame(
        t = t_grid,
        y = fun(t_grid, 0,
                s$beta0, s$beta1,
                s$alpha0, s$alpha1),
        Group = "x = 0",
        Scenario = s$scenario,
        Panel = s$panel
      )
      
      d1 <- data.frame(
        t = t_grid,
        y = fun(t_grid, 1,
                s$beta0, s$beta1,
                s$alpha0, s$alpha1),
        Group = "x = 1",
        Scenario = s$scenario,
        Panel = s$panel
      )
      
      rbind(d0, d1)
    })
  )
}

haz_dat  <- build_data(haz_gomp, "hazard")
surv_dat <- build_data(surv_gomp, "survival")
cum_dat  <- build_data(cumhaz_gomp, "cumhaz")

# Hazard ratio data
hr_dat <- do.call(
  rbind,
  lapply(scenarios, function(s) {
    data.frame(
      t = t_grid,
      y = hr_gomp(t_grid, s$beta1, s$alpha1),
      Scenario = s$scenario,
      Panel = s$panel
    )
  })
)

#-----------------------------------------------------------
# Generic plotting function
#-----------------------------------------------------------
plot_fun <- function(dat, ylab, filename, use_group=TRUE) {
  
  aes_base <- if (use_group) {
    aes(x=t, y=y, color=Scenario, linetype=Group)
  } else {
    aes(x=t, y=y, color=Scenario)
  }
  
  p <- ggplot(dat, aes_base) +
    geom_line(linewidth=1) +
    facet_wrap(~Panel, ncol=2, scales="free_y") +
    scale_y_continuous(
      limits = c(0, 1)
    ) +
    labs(
      x = "Time",
      y = ylab
    ) +
    theme_bw(base_size=18) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size=18),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  print(p)
  ggsave(filename, plot = p, width = 12, height = 8)
}


plot_fun_hazard <- function(dat, ylab, filename, use_group=TRUE) {
  
  aes_base <- if (use_group) {
    aes(x=t, y=y, color=Scenario, linetype=Group)
  } else {
    aes(x=t, y=y, color=Scenario)
  }
  
  p <- ggplot(dat, aes_base) +
    geom_line(linewidth=1) +
    facet_wrap(~Panel, ncol=2, scales="free_y") +
 #   scale_y_continuous(
#      limits = c(0, 1)
#    ) +
    labs(
      x = "Time",
      y = ylab
    ) +
    theme_bw(base_size=18) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size=18),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  print(p)
  ggsave(filename, plot = p, width = 12, height = 8)
}
#-----------------------------------------------------------
# Plot builder
#-----------------------------------------------------------
build_plot <- function(dat, ylab, use_group=TRUE) {
  
  aes_base <- if (use_group) {
    aes(x = t, y = y, color = Scenario, linetype = Group)
  } else {
    aes(x = t, y = y, color = Scenario)
  }
  
  ggplot(dat, aes_base) +
    geom_line(linewidth = 1) +
    facet_wrap(~Panel, ncol = 2, scales = "free_y") +
    labs(
      x = "Time",
      y = ylab
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    )+
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

build_plot_hazard <- function(dat, ylab, use_group=TRUE) {
  
  aes_base <- if (use_group) {
    aes(x = t, y = y, color = Scenario, linetype = Group)
  } else {
    aes(x = t, y = y, color = Scenario)
  }
  
  ggplot(dat, aes_base) +
    geom_line(linewidth = 1) +
    facet_wrap(~Panel, ncol = 2, scales = "free_y") +
    labs(
      x = "Time",
      y = ylab
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

#-----------------------------------------------------------
# Generate plots
#-----------------------------------------------------------
p1<-plot_fun_hazard(haz_dat,  "Hazard function",            "hazard_function_gompertz_plot.pdf")
p2<-plot_fun(surv_dat, "Survival function",          "survival_function_gompertz_plot.pdf")
p3<-plot_fun_hazard(hr_dat,   "Hazard ratio",               "hazard_ratio_plot_gompertz.pdf", use_group=FALSE)


library(patchwork)

p1 <- build_plot_hazard(haz_dat,  "Hazard function")
p2 <- build_plot(surv_dat, "Survival function")
p3 <- build_plot_hazard(hr_dat,   "Hazard ratio", use_group = FALSE)

combined_plot <- (p2 / p1 ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

ggplot2::ggsave("Example_Gompertz.pdf", plot = combined_plot, width = 12, height = 9)

















##Aplication



rm(list = ls())

library(ggplot2)
library(patchwork)

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

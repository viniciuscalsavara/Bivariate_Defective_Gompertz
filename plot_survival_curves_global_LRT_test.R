library(ggplot2)
library(dplyr)

## Survival function
Sg <- function(a, scale, t, eps = 1e-10) {
  if (abs(a) < eps) {
    exp(-scale * t)
  } else {
    exp(-(scale / a) * (exp(a * t) - 1))
  }
}

## Function to create scenario
make_scenario <- function(effect) {
  c(
    alpha0  = -0.7077,
    alpha1  = -0.7077 * effect,
    beta0   = -1.2404,
    beta1   = -1.2404 * effect,
    gamma0  = -0.0560,
    gamma1  = -0.0560 * effect,
    lambda0 = -2.3630,
    lambda1 = -2.3630 * effect,
    phi     = -0.4
  )
}

## Effects considered
effects <- c(0, 0.05, 0.10, 0.15, 0.20)

## Time grid
t.grid <- seq(0, 132, by = 0.1)

## Build data for plotting
build_surv_df <- function(effect, t.grid) {
  
  pars <- make_scenario(effect)
  
  ## T1
  a1_g0  <- pars["alpha0"]
  sc1_g0 <- exp(pars["beta0"])
  
  a1_g1  <- pars["alpha0"] + pars["alpha1"]
  sc1_g1 <- exp(pars["beta0"] + pars["beta1"])
  
  ## T2
  a2_g0  <- pars["gamma0"]
  sc2_g0 <- exp(pars["lambda0"])
  
  a2_g1  <- pars["gamma0"] + pars["gamma1"]
  sc2_g1 <- exp(pars["lambda0"] + pars["lambda1"])
  
  data.frame(
    time = rep(t.grid, 4),
    survival = c(
      Sg(a1_g0, sc1_g0, t.grid),
      Sg(a1_g1, sc1_g1, t.grid),
      Sg(a2_g0, sc2_g0, t.grid),
      Sg(a2_g1, sc2_g1, t.grid)
    ),
    outcome = rep(c("T[1]", "T[1]", "T[2]", "T[2]"), each = length(t.grid)),
    group = rep(c("Group 0", "Group 1", "Group 0", "Group 1"), each = length(t.grid)),
    effect = factor(
      sprintf("%.2f", effect),
      levels = sprintf("%.2f", effects)
    )
  )
}

plot_data <- bind_rows(
  lapply(effects, build_surv_df, t.grid = t.grid)
)

p<-ggplot(plot_data, aes(x = time, y = survival, color = group, linetype = group)) +
  geom_line(linewidth = 2) +
  facet_grid(
    #cols = vars(outcome),
    rows = vars(outcome),
    #rows = vars(effect),
    cols = vars(effect),
    labeller = labeller(
      outcome = label_parsed,
      effect = function(x) paste("Relative effect =", x)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = "Time",
    y = "Survival function",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.spacing.y = unit(0.6, "cm"),
    panel.spacing.x = unit(0.5, "cm"),
    panel.grid.minor = element_blank()
  )
p

ggsave(
  filename = "plot_survival_curves_global_lrt_group.pdf",
  plot = p,
  width =20 , height = 12, dpi = 300
)


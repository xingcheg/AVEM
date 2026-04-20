library(ggplot2)
library(gridExtra)
library(grid)

source("main.R")

# ------------------------------------------------------------
# 1. Accuracy experiment
#    Fix K and d, vary n, T_len, tau2
# ------------------------------------------------------------
run_accuracy_experiment <- function(
    K_fixed = 3,
    d_fixed = 2,
    n_grid = c(20, 35, 50, 75, 100),
    T_grid = c(20, 40, 60, 80, 100),
    tau2_grid = c(0.25, 0.5, 1, 1.5, 2),
    n_rep = 20,
    max_iter = 100,
    tol = 1e-4,
    base_seed = 2026,
    verbose = TRUE,
    fix_tau2 = FALSE,
    fix_sigma2 = FALSE
) {
  scenarios <- expand.grid(
    K = K_fixed,
    n = n_grid,
    T_len = T_grid,
    d = d_fixed,
    tau2 = tau2_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  out <- vector("list", nrow(scenarios) * n_rep)
  idx <- 1

  for (s in seq_len(nrow(scenarios))) {
    sc <- scenarios[s, ]
    for (r in seq_len(n_rep)) {
      if (verbose) {
        cat(sprintf(
          "ACC | Scenario %d/%d | rep %d/%d | K=%d d=%d n=%d T=%d tau2=%.2f\n",
          s, nrow(scenarios), r, n_rep,
          sc$K, sc$d, sc$n, sc$T_len, sc$tau2
        ))
      }

      out[[idx]] <- run_one_replication(
        K = sc$K,
        n = sc$n,
        T_len = sc$T_len,
        d = sc$d,
        tau2 = sc$tau2,
        max_iter = max_iter,
        tol = tol,
        seed = base_seed + 1000 * s + r,
        verbose = FALSE,
        fix_tau2 = fix_tau2,
        fix_sigma2 = fix_sigma2
      )
      idx <- idx + 1
    }
  }

  res <- do.call(rbind, out)
  return(list(raw = res, summary = aggregate_results(res)))
}



# ------------------------------------------------------------
# 7. Example workflow
# ------------------------------------------------------------
# ---- accuracy experiment ----
acc_out <- run_accuracy_experiment(
  K_fixed = 3,
  d_fixed = 2,
  n_grid = c(20, 40, 60, 80, 100),
  T_grid = c(20, 40, 60, 80),
  tau2_grid = c(0.25, 0.5, 1, 2),
  n_rep = 100,
  max_iter = 80,
  tol = 1e-4,
  base_seed = 2026,
  verbose = TRUE
)
acc_out$summary
#saveRDS(acc_out, file = "accuracy_experiment.rds")




acc_out <- readRDS("accuracy_experiment.rds")
# ------------------------------------------------------------
# Add time per iteration
# ------------------------------------------------------------
acc_out$summary$time_per_iter <- acc_out$summary$time_sec / pmax(acc_out$summary$iterations, 1)

# ------------------------------------------------------------
# Common theme
# ------------------------------------------------------------
base_acc_theme <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.2, "cm"),
    legend.spacing.x = unit(0.4, "cm"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
# ------------------------------------------------------------
# Helper to extract legend
# ------------------------------------------------------------
get_legend <- function(myplot) {
  g <- ggplotGrob(myplot)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx]]
}

# ------------------------------------------------------------
# Common x scale
# ------------------------------------------------------------
x_breaks <- sort(unique(acc_out$summary$n))

# ------------------------------------------------------------
# Plot 1: mu
# ------------------------------------------------------------
p_mu <- ggplot(
  acc_out$summary,
  aes(x = n, y = mu_rmse, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = expression("RMSE of " * mu),
    colour = expression(tau^2),
    title = expression(mu)
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Plot 2: sigma^2
# ------------------------------------------------------------
p_sigma <- ggplot(
  acc_out$summary,
  aes(x = n, y = sigma2_rmse, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = expression("RMSE of " * sigma^2),
    colour = expression(tau^2),
    title = expression(sigma^2)
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Plot 3: random effects
# ------------------------------------------------------------
p_f <- ggplot(
  acc_out$summary,
  aes(x = n, y = f_mse, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = expression("MSE of " * f),
    colour = expression(tau^2),
    title = expression(f)
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Plot 4: Gamma
# ------------------------------------------------------------
p_Gamma <- ggplot(
  acc_out$summary,
  aes(x = n, y = gamma_abs, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = expression("MAE of " * Gamma),
    colour = expression(tau^2),
    title = expression(Gamma)
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Plot 5: total time
# ------------------------------------------------------------
p_time <- ggplot(
  acc_out$summary,
  aes(x = n, y = time_sec, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = "Total time (sec)",
    colour = expression(tau^2),
    title = "Total time"
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Plot 6: time per iteration
# ------------------------------------------------------------
p_time_iter <- ggplot(
  acc_out$summary,
  aes(x = n, y = time_per_iter, colour = factor(tau2), group = factor(tau2))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ T_len, nrow = 1,
    labeller = as_labeller(function(x) paste("T =", x))
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Sample size n",
    y = "Time per iteration (sec)",
    colour = expression(tau^2),
    title = "Time per iteration"
  ) +
  base_acc_theme

# ------------------------------------------------------------
# Extract one common legend
# ------------------------------------------------------------
legend_shared <- get_legend(p_mu)

# remove legends from all panels
p_mu_noleg        <- p_mu + theme(legend.position = "none")
p_sigma_noleg     <- p_sigma + theme(legend.position = "none")
p_f_noleg         <- p_f + theme(legend.position = "none")
p_Gamma_noleg     <- p_Gamma + theme(legend.position = "none")
p_time_noleg      <- p_time + theme(legend.position = "none")
p_time_iter_noleg <- p_time_iter + theme(legend.position = "none")

# ------------------------------------------------------------
# Arrange 6 plots + 1 legend
# ------------------------------------------------------------
grid.arrange(
  arrangeGrob(
    p_mu_noleg, p_sigma_noleg,
    p_f_noleg, p_Gamma_noleg,
    p_time_noleg, p_time_iter_noleg,
    ncol = 2
  ),
  legend_shared,
  ncol = 1,
  heights = c(18, 1)
)


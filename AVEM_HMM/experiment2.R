library(ggplot2)
library(gridExtra)
library(grid)

source("main.R")

# ------------------------------------------------------------
# Helper: summarize speed experiment
# ------------------------------------------------------------
aggregate_speed_results <- function(res) {
  res$time_per_iter <- res$time_sec / pmax(res$iterations, 1)

  out <- aggregate(
    cbind(time_sec, time_per_iter, iterations) ~ K + n + T_len + d + tau2,
    data = res,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  out <- out[order(out$K, out$d), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Speed experiment
# Fix n, T_len, tau2; vary K and d
# ------------------------------------------------------------
run_speed_experiment <- function(
    n_fixed = 50,
    T_fixed = 60,
    tau2_fixed = 1,
    K_grid = c(2, 3, 4, 5, 6),
    d_grid = c(1, 2, 3, 5, 8, 10),
    n_rep = 20,
    max_iter = 100,
    tol = 1e-4,
    base_seed = 4040,
    verbose = TRUE,
    fix_tau2 = FALSE,
    fix_sigma2 = FALSE
) {
  scenarios <- expand.grid(
    K = K_grid,
    n = n_fixed,
    T_len = T_fixed,
    d = d_grid,
    tau2 = tau2_fixed,
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
          "SPD | Scenario %d/%d | rep %d/%d | K=%d d=%d n=%d T=%d tau2=%.2f\n",
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
  res$time_per_iter <- res$time_sec / pmax(res$iterations, 1)

  list(
    raw = res,
    summary = aggregate_speed_results(res)
  )
}



# ------------------------------------------------------------
# Example workflow
# ------------------------------------------------------------
speed_out <- run_speed_experiment(
  n_fixed = 60,
  T_fixed = 60,
  tau2_fixed = 1,
  K_grid = c(2, 3, 4, 5, 6, 7),
  d_grid = c(2, 3, 4, 5),
  n_rep = 100,
  max_iter = 80,
  tol = 1e-4,
  base_seed = 4040,
  verbose = TRUE
)

# summarized table
speed_out$summary
#saveRDS(speed_out, file = "speed_experiment.rds")




speed_out <- readRDS("speed_experiment.rds")

# ------------------------------------------------------------
# Rebuild summary if needed
# ------------------------------------------------------------
speed_out$raw$time_per_iter <- speed_out$raw$time_sec / pmax(speed_out$raw$iterations, 1)

speed_out$summary <- aggregate(
  cbind(mu_rmse, sigma2_rmse, f_mse, gamma_abs,
        time_sec, time_per_iter, iterations) ~ K + n + T_len + d + tau2,
  data = speed_out$raw,
  FUN = function(x) mean(x, na.rm = TRUE)
)

speed_out$summary <- speed_out$summary[order(speed_out$summary$K, speed_out$summary$d), ]
rownames(speed_out$summary) <- NULL

# ------------------------------------------------------------
# Common theme
# ------------------------------------------------------------
base_spd_theme <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
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
K_breaks <- sort(unique(speed_out$summary$K))

# ------------------------------------------------------------
# Plot 1: mu
# ------------------------------------------------------------
p_mu <- ggplot(
  speed_out$summary,
  aes(x = K, y = mu_rmse, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = expression("RMSE of " * mu),
    colour = "Dimension d",
    title = expression(mu)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Plot 2: sigma^2
# ------------------------------------------------------------
p_sigma <- ggplot(
  speed_out$summary,
  aes(x = K, y = sigma2_rmse, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = expression("RMSE of " * sigma^2),
    colour = "Dimension d",
    title = expression(sigma^2)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Plot 3: random effects
# ------------------------------------------------------------
p_f <- ggplot(
  speed_out$summary,
  aes(x = K, y = f_mse, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = expression("MSE of " * f),
    colour = "Dimension d",
    title = expression(f)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Plot 4: Gamma
# ------------------------------------------------------------
p_Gamma <- ggplot(
  speed_out$summary,
  aes(x = K, y = gamma_abs, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = expression("MAE of " * Gamma),
    colour = "Dimension d",
    title = expression(Gamma)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Plot 5: total time
# ------------------------------------------------------------
p_time <- ggplot(
  speed_out$summary,
  aes(x = K, y = time_sec, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = "Total time (sec)",
    colour = "Dimension d",
    title = "Total time"
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Plot 6: time per iteration
# ------------------------------------------------------------
p_time_iter <- ggplot(
  speed_out$summary,
  aes(x = K, y = time_per_iter, colour = factor(d), group = factor(d))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = K_breaks) +
  labs(
    x = "Number of latent states K",
    y = "Time per iteration (sec)",
    colour = "Dimension d",
    title = "Time per iteration"
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8, size = 4))) +
  base_spd_theme

# ------------------------------------------------------------
# Shared legend
# ------------------------------------------------------------
legend_shared <- get_legend(p_mu)

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
    ncol = 3
  ),
  legend_shared,
  ncol = 1,
  heights = c(18, 1.2)
)

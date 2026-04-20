source("avem_messm_simulation.R")
source("avem_messm_mc.R")


base_config <- config
base_config$seed <- 2026
base_config$n <- 50
base_config$T <- 100
base_config$p <- 4
base_config$q <- 2
#base_config$tau_g2 <- 0.05
#base_config$tau_h2 <- 0.05
#base_config$r_diag <- rep(0.25, base_config$p)
base_config$max_iter <- 120
base_config$rel_tol <- 1e-4
base_config$verbose <- TRUE

mc_out <- run_monte_carlo_avem_all_curves(
  base_config = base_config,
  n_rep = 100,
  seed_start = 2026,
  save_each_rep = FALSE,
  verbose = TRUE
)


saveRDS(mc_out, "results/case6_out.rds")


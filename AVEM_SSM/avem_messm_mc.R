# ============================================================
# Monte Carlo wrapper for AVEM MESSM:
# store all ELBO curves + all metrics, and plot all normalized
# ============================================================

run_monte_carlo_avem_all_curves <- function(base_config,
                                            n_rep = 50,
                                            seed_start = NULL,
                                            save_each_rep = FALSE,
                                            save_prefix = "avem_messm_mc",
                                            verbose = TRUE) {
  required_fns <- c(
    "make_lower_maps", "make_default_truth", "simulate_messm_data",
    "run_avem_messm", "compute_metrics"
  )
  miss <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function")]
  if (length(miss) > 0) {
    stop("Missing required functions: ", paste(miss, collapse = ", "),
         ". Please source `avem_messm_simulation.R` first.")
  }
  
  if (is.null(seed_start)) {
    seed_start <- if (!is.null(base_config$seed)) base_config$seed else 1L
  }
  
  metrics_list <- vector("list", n_rep)
  elbo_long_list <- vector("list", n_rep)
  elapsed_vec <- numeric(n_rep)
  
  if (save_each_rep) {
    dir.create(paste0(save_prefix, "_rep_outputs"), showWarnings = FALSE, recursive = TRUE)
  }
  
  mc_start <- proc.time()[3]
  
  for (r in seq_len(n_rep)) {
    cfg <- base_config
    cfg$seed <- seed_start + r - 1L
    cfg$verbose <- FALSE
    cfg$save_prefix <- file.path(paste0(save_prefix, "_rep_outputs"),
                                 paste0("rep", sprintf("%03d", r)))
    
    maps <- make_lower_maps(cfg$p, cfg$q)
    truth <- make_default_truth(
      p = cfg$p,
      q = cfg$q,
      tau_g2 = cfg$tau_g2,
      tau_h2 = cfg$tau_h2,
      r_diag = cfg$r_diag,
      maps = maps
    )
    
    rep_start <- proc.time()[3]
    sim_data <- simulate_messm_data(cfg, truth, maps)
    fit <- run_avem_messm(sim_data, cfg, maps)
    rep_elapsed <- proc.time()[3] - rep_start
    met <- compute_metrics(fit, sim_data, cfg, maps)
    
    met$replicate_id <- r
    met$seed <- cfg$seed
    met$n <- cfg$n
    met$T <- cfg$T
    met$p <- cfg$p
    met$q <- cfg$q
    met$tau_g2 <- cfg$tau_g2
    met$tau_h2 <- cfg$tau_h2
    met$wall_time_rep_sec <- rep_elapsed
    
    metrics_list[[r]] <- met
    elapsed_vec[r] <- rep_elapsed
    
    elbo_vec <- fit$history$elbo
    rel_vec <- fit$history$rel_change
    n_iter_r <- length(elbo_vec)
    
    if (n_iter_r >= 1L) {
      iter_idx <- seq_len(n_iter_r)
      
      elbo_long_list[[r]] <- data.frame(
        replicate_id = r,
        seed = cfg$seed,
        iteration = iter_idx,
        elbo = as.numeric(elbo_vec),
        rel_change = as.numeric(rel_vec),
        n = cfg$n,
        T = cfg$T,
        p = cfg$p,
        q = cfg$q,
        tau_g2 = cfg$tau_g2,
        tau_h2 = cfg$tau_h2
      )
    } else {
      elbo_long_list[[r]] <- data.frame(
        replicate_id = integer(0),
        seed = integer(0),
        iteration = integer(0),
        elbo = numeric(0),
        rel_change = numeric(0),
        n = integer(0),
        T = integer(0),
        p = integer(0),
        q = integer(0),
        tau_g2 = numeric(0),
        tau_h2 = numeric(0)
      )
    }
    
    if (save_each_rep) {
      saveRDS(
        list(config = cfg, truth = truth, sim_data = sim_data, fit = fit,
             metrics = met, elbo_history = elbo_long_list[[r]]),
        paste0(cfg$save_prefix, "_full_output.rds")
      )
    }
    
    if (verbose) {
      cat(sprintf(
        "Rep %3d/%3d done | seed = %d | iter = %d | rep time = %.3f sec | latent MSE = %.6f\n",
        r, n_rep, cfg$seed, fit$n_iter, rep_elapsed, met$latent_mse
      ))
    }
  }
  
  metrics_df <- do.call(rbind, metrics_list)
  elbo_df <- do.call(rbind, elbo_long_list)
  
  metric_names <- setdiff(names(metrics_df), c(
    "replicate_id", "seed", "n", "T", "p", "q", "tau_g2", "tau_h2"
  ))
  
  summary_df <- rbind(
    cbind(summary_type = "mean", as.data.frame(as.list(colMeans(metrics_df[, metric_names, drop = FALSE], na.rm = TRUE)))),
    cbind(summary_type = "sd", as.data.frame(as.list(vapply(metrics_df[, metric_names, drop = FALSE], sd, numeric(1), na.rm = TRUE)))),
    cbind(summary_type = "q25", as.data.frame(as.list(vapply(metrics_df[, metric_names, drop = FALSE], quantile, numeric(1), probs = 0.25, na.rm = TRUE)))),
    cbind(summary_type = "median", as.data.frame(as.list(vapply(metrics_df[, metric_names, drop = FALSE], median, numeric(1), na.rm = TRUE)))),
    cbind(summary_type = "q75", as.data.frame(as.list(vapply(metrics_df[, metric_names, drop = FALSE], quantile, numeric(1), probs = 0.75, na.rm = TRUE))))
  )
  
  mc_total_time <- proc.time()[3] - mc_start
  meta_df <- data.frame(
    n_rep = n_rep,
    seed_start = seed_start,
    n = base_config$n,
    T = base_config$T,
    p = base_config$p,
    q = base_config$q,
    tau_g2 = base_config$tau_g2,
    tau_h2 = base_config$tau_h2,
    mc_total_time_sec = mc_total_time,
    avg_rep_time_sec = mean(elapsed_vec),
    avg_iter = mean(metrics_df$n_iter),
    max_iter_observed = max(metrics_df$n_iter),
    save_each_rep = save_each_rep
  )
  
  if (verbose) {
    cat("\n================ Monte Carlo Summary ================\n")
    print(summary_df)
    cat(sprintf("\nMonte Carlo total time: %.3f sec\n", mc_total_time))
  }
  
  out <- list(
    metrics = metrics_df,
    summary_metrics = summary_df,
    elbo_histories = elbo_df,
    meta = meta_df,
    save_prefix = save_prefix
  )
  class(out) <- c("avem_messm_mc_elbo_curves", class(out))
  invisible(out)
}



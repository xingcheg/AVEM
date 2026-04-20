# ================================================================
# Gaussian mixed HMM: mean-anchor AVEM only, generalized to K >= 2
#
# Model:
#   U_it in {1,...,K} follows a common K-state Markov chain
#   f_i ~ N_d(0, tau2 * I_d)
#   Y_it | (U_it = k, f_i) ~ N_d(mu_k + f_i, sigma2_k * I_d)
# ================================================================

# -----------------------------
# Utilities
# -----------------------------
softmax_log <- function(logw) {
  a <- max(logw)
  w <- exp(logw - a)
  w / sum(w)
}

safe_cor <- function(x, y) {
  if (length(x) <= 1 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}

# General stationary distribution for K-state Markov chain
stationary_dist <- function(Gamma) {
  K <- nrow(Gamma)
  A <- t(Gamma) - diag(K)
  A[K, ] <- 1
  b <- c(rep(0, K - 1), 1)
  delta <- as.numeric(solve(A, b))
  delta <- pmax(delta, 0)
  delta / sum(delta)
}

make_true_pars <- function(K = 3, d = 1, tau2 = 1,
                           mu_centers = NULL,
                           sigma2 = NULL,
                           Gamma = NULL) {
  if (is.null(mu_centers)) {
    # default: spread means roughly symmetrically
    centers <- seq(from = 1.5, to = -1.5, length.out = K)
    mu <- lapply(centers, function(a) rep(a, d))
  } else {
    stopifnot(length(mu_centers) == K)
    mu <- lapply(mu_centers, function(a) rep(a, d))
  }
  
  if (is.null(sigma2)) {
    sigma2 <- rep(1, K)
  } else {
    stopifnot(length(sigma2) == K)
  }
  
  if (is.null(Gamma)) {
    offdiag <- 0.08 / max(K - 1, 1)
    Gamma <- matrix(offdiag, K, K)
    diag(Gamma) <- 0.92
    Gamma <- Gamma / rowSums(Gamma)
  } else {
    stopifnot(all(dim(Gamma) == c(K, K)))
    Gamma <- Gamma / rowSums(Gamma)
  }
  
  list(
    mu = mu,
    sigma2 = sigma2,
    tau2 = tau2,
    Gamma = Gamma,
    delta = stationary_dist(Gamma),
    d = d,
    K = K
  )
}

copy_pars <- function(pars) {
  K <- length(pars$mu)
  list(
    mu = lapply(pars$mu, function(x) as.numeric(x)),
    sigma2 = as.numeric(pars$sigma2),
    tau2 = as.numeric(pars$tau2),
    Gamma = matrix(pars$Gamma, nrow = K, ncol = K),
    delta = as.numeric(pars$delta),
    d = as.integer(pars$d),
    K = as.integer(K)
  )
}

# -----------------------------
# Relabeling for general K
# Match estimated labels to truth by nearest means (greedy / exhaustive for small K)
# -----------------------------
all_permutations <- function(x) {
  if (length(x) == 1) return(list(x))
  out <- list()
  idx <- 1
  for (i in seq_along(x)) {
    rest <- x[-i]
    perms_rest <- all_permutations(rest)
    for (pr in perms_rest) {
      out[[idx]] <- c(x[i], pr)
      idx <- idx + 1
    }
  }
  out
}

find_best_permutation_by_mu <- function(est_mu, true_mu) {
  K <- length(est_mu)
  stopifnot(length(true_mu) == K)
  
  perms <- all_permutations(seq_len(K))
  best_perm <- seq_len(K)
  best_obj <- Inf
  
  for (perm in perms) {
    obj <- 0
    for (k in seq_len(K)) {
      obj <- obj + sum((est_mu[[perm[k]]] - true_mu[[k]])^2)
    }
    if (obj < best_obj) {
      best_obj <- obj
      best_perm <- perm
    }
  }
  best_perm
}

relabel_fit <- function(fit, perm) {
  out <- fit
  K <- length(perm)
  
  out$pars$mu <- out$pars$mu[perm]
  out$pars$sigma2 <- out$pars$sigma2[perm]
  out$pars$Gamma <- out$pars$Gamma[perm, perm, drop = FALSE]
  out$pars$delta <- out$pars$delta[perm]
  
  inv_perm <- integer(K)
  inv_perm[perm] <- seq_len(K)
  
  if (!is.null(out$gamma_list)) {
    out$gamma_list <- lapply(out$gamma_list, function(g) g[, perm, drop = FALSE])
  }
  if (!is.null(out$xi_list)) {
    out$xi_list <- lapply(out$xi_list, function(xi) xi[, perm, perm, drop = FALSE])
  }
  if (!is.null(out$state_hat)) {
    out$state_hat <- lapply(out$state_hat, function(u) inv_perm[u])
  }
  
  out
}

match_labels_general <- function(fit, truth) {
  perm <- find_best_permutation_by_mu(fit$pars$mu, truth$mu)
  relabel_fit(fit, perm)
}

# -----------------------------
# Data generation
# -----------------------------
simulate_gaussian_mhmm <- function(n, T_len, pars, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  K <- length(pars$mu)
  d <- pars$d
  Gamma <- pars$Gamma
  delta <- pars$delta
  mu <- pars$mu
  sigma2 <- pars$sigma2
  tau2 <- pars$tau2
  
  Y_list <- vector("list", n)
  U_list <- vector("list", n)
  F_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    f_i <- rnorm(d, mean = 0, sd = sqrt(tau2))
    U <- integer(T_len)
    Y <- matrix(0, nrow = T_len, ncol = d)
    
    U[1] <- sample.int(K, size = 1, prob = delta)
    Y[1, ] <- rnorm(d, mean = mu[[U[1]]] + f_i, sd = sqrt(sigma2[U[1]]))
    
    if (T_len >= 2) {
      for (t in 2:T_len) {
        U[t] <- sample.int(K, size = 1, prob = Gamma[U[t - 1], ])
        Y[t, ] <- rnorm(d, mean = mu[[U[t]]] + f_i, sd = sqrt(sigma2[U[t]]))
      }
    }
    
    Y_list[[i]] <- Y
    U_list[[i]] <- U
    F_list[[i]] <- f_i
  }
  
  list(Y = Y_list, U = U_list, F = F_list, pars = pars, n = n, T_len = T_len, d = d, K = K)
}

# -----------------------------
# HMM routines conditional on fixed f
# -----------------------------
log_emission_matrix <- function(y, f, pars) {
  T_len <- nrow(y)
  d <- ncol(y)
  K <- length(pars$mu)
  out <- matrix(0, nrow = T_len, ncol = K)
  
  for (k in seq_len(K)) {
    mean_mat <- matrix(pars$mu[[k]] + f, nrow = T_len, ncol = d, byrow = TRUE)
    sq <- rowSums((y - mean_mat)^2)
    out[, k] <- -0.5 * d * log(2 * pi * pars$sigma2[k]) - 0.5 * sq / pars$sigma2[k]
  }
  out
}

fb_conditional <- function(y, f, pars) {
  Gamma <- pars$Gamma
  delta <- pars$delta
  T_len <- nrow(y)
  K <- length(delta)
  
  log_emit <- log_emission_matrix(y, f, pars)
  emit <- exp(log_emit)
  emit[emit < 1e-300] <- 1e-300
  
  alpha <- matrix(0, nrow = T_len, ncol = K)
  beta <- matrix(0, nrow = T_len, ncol = K)
  cscale <- numeric(T_len)
  
  alpha[1, ] <- delta * emit[1, ]
  cscale[1] <- sum(alpha[1, ])
  alpha[1, ] <- alpha[1, ] / cscale[1]
  
  if (T_len >= 2) {
    for (t in 2:T_len) {
      alpha[t, ] <- as.numeric(alpha[t - 1, ] %*% Gamma) * emit[t, ]
      cscale[t] <- sum(alpha[t, ])
      alpha[t, ] <- alpha[t, ] / cscale[t]
    }
  }
  
  beta[T_len, ] <- 1
  if (T_len >= 2) {
    for (t in (T_len - 1):1) {
      beta[t, ] <- as.numeric(Gamma %*% (emit[t + 1, ] * beta[t + 1, ])) / cscale[t + 1]
    }
  }
  
  gamma <- alpha * beta
  gamma <- gamma / rowSums(gamma)
  
  xi <- array(0, dim = c(max(T_len - 1, 1), K, K))
  if (T_len >= 2) {
    for (t in 1:(T_len - 1)) {
      numer <- outer(alpha[t, ], emit[t + 1, ] * beta[t + 1, ]) * Gamma
      xi[t, , ] <- numer / sum(numer)
    }
  }
  
  list(logLik = sum(log(cscale)), gamma = gamma, xi = xi)
}

viterbi_conditional <- function(y, f, pars) {
  T_len <- nrow(y)
  K <- length(pars$delta)
  log_emit <- log_emission_matrix(y, f, pars)
  log_delta <- log(pmax(pars$delta, 1e-300))
  log_Gamma <- log(pmax(pars$Gamma, 1e-300))
  
  delta_mat <- matrix(-Inf, nrow = T_len, ncol = K)
  psi <- matrix(0L, nrow = T_len, ncol = K)
  delta_mat[1, ] <- log_delta + log_emit[1, ]
  
  if (T_len >= 2) {
    for (t in 2:T_len) {
      for (k in seq_len(K)) {
        vals <- delta_mat[t - 1, ] + log_Gamma[, k]
        psi[t, k] <- which.max(vals)
        delta_mat[t, k] <- max(vals) + log_emit[t, k]
      }
    }
  }
  
  path <- integer(T_len)
  path[T_len] <- which.max(delta_mat[T_len, ])
  if (T_len >= 2) {
    for (t in (T_len - 1):1) path[t] <- psi[t + 1, path[t + 1]]
  }
  path
}

# -----------------------------
# Parameter updates
# -----------------------------
update_delta <- function(gamma_list) {
  delta <- Reduce(`+`, lapply(gamma_list, function(g) g[1, ])) / length(gamma_list)
  delta / sum(delta)
}

update_Gamma <- function(gamma_list, xi_list) {
  K <- ncol(gamma_list[[1]])
  numer <- matrix(0, K, K)
  denom <- numeric(K)
  
  for (i in seq_along(gamma_list)) {
    g <- gamma_list[[i]]
    x <- xi_list[[i]]
    if (nrow(g) >= 2) {
      numer <- numer + apply(x, c(2, 3), sum)
      denom <- denom + colSums(g[1:(nrow(g) - 1), , drop = FALSE])
    }
  }
  
  G <- numer / pmax(denom, 1e-12)
  G <- G / rowSums(G)
  G
}

update_mu_sigma2_avem <- function(Y_list, gamma_list, m_list, V_list, pars_old) {
  K <- length(pars_old$mu)
  d <- pars_old$d
  mu_new <- vector("list", K)
  sigma2_new <- numeric(K)
  
  for (k in seq_len(K)) {
    num_mu <- rep(0, d)
    den_mu <- 0
    for (i in seq_along(Y_list)) {
      y <- Y_list[[i]]
      gk <- gamma_list[[i]][, k]
      num_mu <- num_mu + colSums(gk * (y - matrix(m_list[[i]], nrow = nrow(y), ncol = d, byrow = TRUE)))
      den_mu <- den_mu + sum(gk)
    }
    mu_new[[k]] <- num_mu / max(den_mu, 1e-12)
  }
  
  for (k in seq_len(K)) {
    numer <- 0
    denom <- 0
    for (i in seq_along(Y_list)) {
      y <- Y_list[[i]]
      T_len <- nrow(y)
      gk <- gamma_list[[i]][, k]
      m_i <- m_list[[i]]
      V_i <- V_list[[i]]
      mu_k <- mu_new[[k]]
      resid_mean <- y - matrix(mu_k + m_i, nrow = T_len, ncol = d, byrow = TRUE)
      sq <- rowSums(resid_mean^2)
      numer <- numer + sum(gk * (sq + sum(diag(V_i))))
      denom <- denom + d * sum(gk)
    }
    sigma2_new[k] <- numer / max(denom, 1e-12)
    sigma2_new[k] <- max(sigma2_new[k], 1e-6)
  }
  
  list(mu = mu_new, sigma2 = sigma2_new)
}

update_tau2_avem <- function(m_list, V_list, d) {
  val <- 0
  n <- length(m_list)
  for (i in seq_len(n)) val <- val + sum(m_list[[i]]^2) + sum(diag(V_list[[i]]))
  max(val / (n * d), 1e-6)
}

# -----------------------------
# Mean-anchor AVEM
# -----------------------------
fit_mean_anchor_avem <- function(dat, init_pars = NULL, max_iter = 100, tol = 1e-4, verbose = FALSE,
                                 fix_tau2 = FALSE, tau2_value = NULL,
                                 fix_sigma2 = FALSE, sigma2_value = NULL) {
  t0 <- proc.time()[3]
  n <- dat$n
  d <- dat$d
  K <- dat$K
  Y_list <- dat$Y
  
  if (is.null(init_pars)) {
    pars <- make_true_pars(K = K, d = d, tau2 = 0.5)
    centers <- seq(from = 0.8, to = -0.8, length.out = K)
    pars$mu <- lapply(centers, function(a) rep(a, d))
    pars$sigma2 <- rep(1.2, K)
  } else {
    pars <- copy_pars(init_pars)
  }
  
  if (fix_tau2) {
    if (is.null(tau2_value)) tau2_value <- dat$pars$tau2
    pars$tau2 <- tau2_value
  }
  if (fix_sigma2) {
    if (is.null(sigma2_value)) sigma2_value <- dat$pars$sigma2
    pars$sigma2 <- sigma2_value
  }
  
  m_list <- replicate(n, rep(0, d), simplify = FALSE)
  V_list <- replicate(n, diag(pars$tau2, d), simplify = FALSE)
  prev_obj <- -Inf
  iter_used <- 0
  
  for (iter in seq_len(max_iter)) {
    gamma_list <- vector("list", n)
    xi_list <- vector("list", n)
    logLik_sum <- 0
    
    # E-step for latent states at anchor m_i
    for (i in seq_len(n)) {
      fb <- fb_conditional(Y_list[[i]], m_list[[i]], pars)
      gamma_list[[i]] <- fb$gamma
      xi_list[[i]] <- fb$xi
      logLik_sum <- logLik_sum + fb$logLik
    }
    
    # Update q_i(f_i) approximately using gamma at anchor
    inv_tau2 <- 1 / pars$tau2
    for (i in seq_len(n)) {
      y <- Y_list[[i]]
      T_len <- nrow(y)
      
      weight_sum <- 0
      rhs <- rep(0, d)
      
      for (k in seq_len(K)) {
        wk <- gamma_list[[i]][, k] / pars$sigma2[k]
        weight_sum <- weight_sum + sum(wk)
        rhs <- rhs + colSums(wk * (y - matrix(pars$mu[[k]], nrow = T_len, ncol = d, byrow = TRUE)))
      }
      
      prec_scalar <- inv_tau2 + weight_sum
      v_scalar <- 1 / prec_scalar
      m_list[[i]] <- v_scalar * rhs
      V_list[[i]] <- diag(v_scalar, d)
    }
    
    # M-step
    pars$delta <- update_delta(gamma_list)
    pars$Gamma <- update_Gamma(gamma_list, xi_list)
    
    upd <- update_mu_sigma2_avem(Y_list, gamma_list, m_list, V_list, pars)
    pars$mu <- upd$mu
    
    if (fix_sigma2) {
      pars$sigma2 <- sigma2_value
    } else {
      pars$sigma2 <- upd$sigma2
    }
    
    if (fix_tau2) {
      pars$tau2 <- tau2_value
    } else {
      pars$tau2 <- update_tau2_avem(m_list, V_list, d)
    }
    
    obj <- logLik_sum
    iter_used <- iter
    if (verbose) cat(sprintf("AVEM iter %d, obj = %.4f, tau2 = %.4f\n", iter, obj, pars$tau2))
    if (iter > 1 && abs(obj - prev_obj) / (abs(prev_obj) + 1e-8) < tol) break
    prev_obj <- obj
  }
  
  state_hat <- vector("list", n)
  for (i in seq_len(n)) {
    state_hat[[i]] <- viterbi_conditional(Y_list[[i]], m_list[[i]], pars)
  }
  
  list(
    method = "AVEM",
    pars = pars,
    m_list = m_list,
    V_list = V_list,
    gamma_list = gamma_list,
    xi_list = xi_list,
    state_hat = state_hat,
    iterations = iter_used,
    total_nodes = 1,
    time_sec = proc.time()[3] - t0
  )
}

# -----------------------------
# Evaluation
# -----------------------------
evaluate_fit <- function(fit, truth, dat) {
  est <- fit$pars
  K <- truth$K
  d <- truth$d
  
  mu_rmse <- sqrt(mean(unlist(lapply(seq_len(K), function(k) {
    (est$mu[[k]] - truth$mu[[k]])^2
  }))))
  
  sigma2_rmse <- sqrt(mean((est$sigma2 - truth$sigma2)^2))
  tau2_abs <- abs(est$tau2 - truth$tau2)
  gamma_abs <- mean(abs(est$Gamma - truth$Gamma))
  
  f_mse <- mean(sapply(seq_along(dat$F), function(i) {
    mean((fit$m_list[[i]] - dat$F[[i]])^2)
  }))
  f_cor_mean <- mean(sapply(seq_along(dat$F), function(i) {
    safe_cor(fit$m_list[[i]], dat$F[[i]])
  }), na.rm = TRUE)
  
  state_acc <- mean(unlist(lapply(seq_along(dat$U), function(i) {
    fit$state_hat[[i]] == dat$U[[i]]
  })))
  
  data.frame(
    method = fit$method,
    K = K,
    n = dat$n,
    T_len = dat$T_len,
    d = d,
    tau2 = truth$tau2,
    time_sec = fit$time_sec,
    iterations = fit$iterations,
    total_nodes = fit$total_nodes,
    mu_rmse = mu_rmse,
    sigma2_rmse = sigma2_rmse,
    tau2_abs = tau2_abs,
    gamma_abs = gamma_abs,
    f_mse = f_mse,
    f_cor_mean = f_cor_mean,
    state_acc = state_acc,
    stringsAsFactors = FALSE
  )
}

print_fit_summary <- function(fit) {
  K <- length(fit$pars$mu)
  mu_means <- sapply(fit$pars$mu, mean)
  cat(sprintf("  %-6s | K=%d | iter=%3d | time=%8.3f sec | nodes=%6d | tau2=%.4f\n",
              fit$method, K, fit$iterations, fit$time_sec, fit$total_nodes, fit$pars$tau2))
  cat("    mean(mu_k):", paste(round(mu_means, 3), collapse = ", "), "\n")
}

# -----------------------------
# One replication
# -----------------------------
run_one_replication <- function(K = 3, n, T_len, d, tau2,
                                max_iter = 100, tol = 1e-4, seed = 1,
                                verbose = FALSE,
                                fix_tau2 = FALSE,
                                fix_sigma2 = FALSE) {
  truth <- make_true_pars(K = K, d = d, tau2 = tau2)
  dat <- simulate_gaussian_mhmm(n = n, T_len = T_len, pars = truth, seed = seed)
  
  init_pars <- make_true_pars(K = K, d = d, tau2 = max(0.4, tau2 * 0.7))
  init_centers <- seq(from = 0.8, to = -0.8, length.out = K)
  init_pars$mu <- lapply(init_centers, function(a) rep(a, d))
  init_pars$sigma2 <- rep(1.2, K)
  
  # mildly sticky init Gamma
  offdiag <- 0.15 / max(K - 1, 1)
  init_pars$Gamma <- matrix(offdiag, K, K)
  diag(init_pars$Gamma) <- 0.85
  init_pars$Gamma <- init_pars$Gamma / rowSums(init_pars$Gamma)
  init_pars$delta <- stationary_dist(init_pars$Gamma)
  
  fit <- fit_mean_anchor_avem(
    dat, init_pars = init_pars, max_iter = max_iter, tol = tol, verbose = verbose,
    fix_tau2 = fix_tau2, tau2_value = truth$tau2,
    fix_sigma2 = fix_sigma2, sigma2_value = truth$sigma2
  )
  
  fit <- match_labels_general(fit, truth)
  
  if (verbose) {
    cat("  Finished fitting AVEM:\n")
    print_fit_summary(fit)
  }
  
  res <- evaluate_fit(fit, truth = truth, dat = dat)
  
  if (verbose) {
    cat("  Accuracy snapshot:\n")
    print(res[, c("method", "K", "time_sec", "mu_rmse", "sigma2_rmse", "tau2_abs", "state_acc")])
    cat("\n")
  }
  
  res
}

# -----------------------------
# Experiment runner
# -----------------------------
run_experiment <- function(n_rep = 10,
                           scenarios = expand.grid(K = c(2, 3, 4),
                                                   n = 80,
                                                   T_len = 60,
                                                   d = c(1, 2, 3),
                                                   tau2 = 1,
                                                   KEEP.OUT.ATTRS = FALSE,
                                                   stringsAsFactors = FALSE),
                           max_iter = 100,
                           tol = 1e-4,
                           base_seed = 2026,
                           verbose = TRUE,
                           fix_tau2 = FALSE,
                           fix_sigma2 = FALSE) {
  out <- vector("list", nrow(scenarios) * n_rep)
  idx <- 1
  
  for (s in seq_len(nrow(scenarios))) {
    sc <- scenarios[s, ]
    for (r in seq_len(n_rep)) {
      if (verbose) {
        cat(sprintf("Scenario %d/%d | rep %d/%d | K=%d n=%d T=%d d=%d tau2=%.2f\n",
                    s, nrow(scenarios), r, n_rep, sc$K, sc$n, sc$T_len, sc$d, sc$tau2))
      }
      
      out[[idx]] <- run_one_replication(
        K = sc$K, n = sc$n, T_len = sc$T_len, d = sc$d, tau2 = sc$tau2,
        max_iter = max_iter, tol = tol,
        seed = base_seed + 1000 * s + r, verbose = verbose,
        fix_tau2 = fix_tau2, fix_sigma2 = fix_sigma2
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, out)
}

aggregate_results <- function(res) {
  wanted <- c("time_sec", "iterations", "total_nodes",
              "mu_rmse", "sigma2_rmse", "tau2_abs",
              "gamma_abs", "f_mse", "f_cor_mean", "state_acc")
  have <- intersect(wanted, names(res))
  if (length(have) == 0) stop("No expected metric columns were found in `res`.")
  
  out <- aggregate(
    res[, have, drop = FALSE],
    by = list(method = res$method, K = res$K, n = res$n, T_len = res$T_len, d = res$d, tau2 = res$tau2),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  out[order(out$K, out$d, out$tau2, out$method), ]
}

# # -----------------------------
# # Example usage
# # -----------------------------
# results_dim <- run_experiment(
#   n_rep = 20,
#   scenarios = expand.grid(
#     K = c(2, 3, 4),
#     n = 50,
#     T_len = c(50, 100),
#     d = c(1, 2, 3),
#     tau2 = 1,
#     KEEP.OUT.ATTRS = FALSE,
#     stringsAsFactors = FALSE
#   ),
#   max_iter = 80,
#   tol = 1e-4,
#   base_seed = 2026,
#   verbose = TRUE
# )
# 
# aggregate_results(results_dim)

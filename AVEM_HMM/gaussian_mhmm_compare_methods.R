# ================================================================
# Gaussian mixed HMM: AVEM vs Gaussian quadrature EM vs MCEM
#
# Model:
#   U_it in {1,...,K} follows a common K-state Markov chain
#   f_i ~ N_d(0, tau2 * I_d)
#   Y_it | (U_it = k, f_i) ~ N_d(mu_k + f_i, sigma2_k * I_d)
#
# This script extends the user's AVEM code to include:
#   1) Mean-anchor AVEM
#   2) Gaussian quadrature EM (QEM)
#   3) Monte Carlo EM (MCEM)
#
# The intended comparison setting is:
#   K = 3, n = 60, T = 60, d in {2, 3}
#   QEM single-dimension quadrature points = 3 or 5
#   MCEM Monte Carlo points = 50 or 100
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

parameter_distance <- function(pars_new, pars_old) {
  K <- pars_new$K
  out <- 0
  for (k in seq_len(K)) out <- max(out, max(abs(pars_new$mu[[k]] - pars_old$mu[[k]])))
  out <- max(out, max(abs(pars_new$sigma2 - pars_old$sigma2)))
  out <- max(out, abs(pars_new$tau2 - pars_old$tau2))
  out <- max(out, max(abs(pars_new$Gamma - pars_old$Gamma)))
  out <- max(out, max(abs(pars_new$delta - pars_old$delta)))
  out
}

# -----------------------------
# Relabeling for general K
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
  perms <- all_permutations(seq_len(K))
  best_perm <- seq_len(K)
  best_obj <- Inf
  for (perm in perms) {
    obj <- 0
    for (k in seq_len(K)) obj <- obj + sum((est_mu[[perm[k]]] - true_mu[[k]])^2)
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

  if (!is.null(out$gamma_list)) out$gamma_list <- lapply(out$gamma_list, function(g) g[, perm, drop = FALSE])
  if (!is.null(out$xi_list)) out$xi_list <- lapply(out$xi_list, function(xi) xi[, perm, perm, drop = FALSE])
  if (!is.null(out$state_hat)) out$state_hat <- lapply(out$state_hat, function(u) inv_perm[u])
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
# Parameter updates from posterior summaries
# -----------------------------
update_delta <- function(gamma_list) {
  delta <- Reduce(`+`, lapply(gamma_list, function(g) g[1, ])) / length(gamma_list)
  delta / sum(delta)
}

update_Gamma <- function(gamma_list, xi_list) {
  K <- ncol(gamma_list[[1]])
  numer <- matrix(0, K, K)
  
  for (i in seq_along(gamma_list)) {
    g <- gamma_list[[i]]
    x <- xi_list[[i]]
    if (nrow(g) >= 2) {
      numer <- numer + apply(x, c(2, 3), sum)
    }
  }
  
  denom <- rowSums(numer)
  
  G <- matrix(0, K, K)
  for (k in seq_len(K)) {
    if (denom[k] > 1e-12) {
      G[k, ] <- numer[k, ] / denom[k]
    } else {
      G[k, ] <- rep(1 / K, K)
    }
  }
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
    sigma2_new[k] <- max(numer / max(denom, 1e-12), 1e-6)
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

    for (i in seq_len(n)) {
      fb <- fb_conditional(Y_list[[i]], m_list[[i]], pars)
      gamma_list[[i]] <- fb$gamma
      xi_list[[i]] <- fb$xi
      logLik_sum <- logLik_sum + fb$logLik
    }

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

    pars$delta <- update_delta(gamma_list)
    pars$Gamma <- update_Gamma(gamma_list, xi_list)

    upd <- update_mu_sigma2_avem(Y_list, gamma_list, m_list, V_list, pars)
    pars$mu <- upd$mu
    pars$sigma2 <- if (fix_sigma2) sigma2_value else upd$sigma2
    pars$tau2 <- if (fix_tau2) tau2_value else update_tau2_avem(m_list, V_list, d)

    obj <- logLik_sum
    iter_used <- iter
    if (verbose) cat(sprintf("AVEM iter %d, obj = %.4f, tau2 = %.4f\n", iter, obj, pars$tau2))
    if (iter > 1 && abs(obj - prev_obj) / (abs(prev_obj) + 1e-8) < tol) break
    prev_obj <- obj
  }

  state_hat <- vector("list", n)
  for (i in seq_len(n)) state_hat[[i]] <- viterbi_conditional(Y_list[[i]], m_list[[i]], pars)

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
# Integration rules for QEM / MCEM
# -----------------------------
std_normal_gh_rule <- function(n_nodes) {
  
  if (length(n_nodes) != 1 || !is.numeric(n_nodes) || n_nodes < 1 || n_nodes != as.integer(n_nodes)) {
    stop("n_nodes must be a positive integer.")
  }
  
  gh <- statmod::gauss.quad.prob(n_nodes, dist = "normal")
  
  list(
    nodes = gh$nodes,
    weights = gh$weights
  )
}

tensor_rule_std_normal <- function(d, n_nodes_1d) {
  rule1 <- std_normal_gh_rule(n_nodes_1d)
  idx <- expand.grid(replicate(d, seq_along(rule1$nodes), simplify = FALSE))
  f_mat <- matrix(0, nrow = nrow(idx), ncol = d)
  w <- rep(1, nrow(idx))
  for (j in seq_len(d)) {
    ids <- idx[[j]]
    f_mat[, j] <- rule1$nodes[ids]
    w <- w * rule1$weights[ids]
  }
  list(z = f_mat, w = w)
}

# -----------------------------
# Posterior summaries given a node/sample set
# Nodes are standard normal points z_r; actual f_r = sqrt(tau2) * z_r.
# Weights w_r integrate w.r.t. N(0, I_d).
# -----------------------------
subject_posterior_summary <- function(y, pars, z_mat, w, method_label = "QEM") {
  R <- nrow(z_mat)
  d <- pars$d
  K <- pars$K
  T_len <- nrow(y)

  f_mat <- sqrt(pars$tau2) * z_mat
  gamma_nodes <- vector("list", R)
  xi_nodes <- vector("list", R)
  logw_post <- numeric(R)

  for (r in seq_len(R)) {
    fb <- fb_conditional(y, f_mat[r, ], pars)
    gamma_nodes[[r]] <- fb$gamma
    xi_nodes[[r]] <- fb$xi
    logw_post[r] <- log(pmax(w[r], 1e-300)) + fb$logLik
  }
  post_w <- softmax_log(logw_post)

  gamma_bar <- matrix(0, nrow = T_len, ncol = K)
  xi_bar <- array(0, dim = c(max(T_len - 1, 1), K, K))
  post_mean <- rep(0, d)
  post_second <- matrix(0, d, d)

  for (r in seq_len(R)) {
    gamma_bar <- gamma_bar + post_w[r] * gamma_nodes[[r]]
    if (T_len >= 2) xi_bar <- xi_bar + post_w[r] * xi_nodes[[r]]
    post_mean <- post_mean + post_w[r] * f_mat[r, ]
    post_second <- post_second + post_w[r] * tcrossprod(f_mat[r, ])
  }
  post_cov <- post_second - tcrossprod(post_mean)
  diag(post_cov) <- pmax(diag(post_cov), 0)

  logLik_i <- max(logw_post) + log(sum(exp(logw_post - max(logw_post))))

  list(
    gamma_bar = gamma_bar,
    xi_bar = xi_bar,
    post_w = post_w,
    f_mat = f_mat,
    gamma_nodes = gamma_nodes,
    xi_nodes = xi_nodes,
    post_mean = post_mean,
    post_cov = post_cov,
    logLik = logLik_i,
    total_nodes = R,
    method_label = method_label
  )
}

# -----------------------------
# M-step from posterior summaries for QEM / MCEM
# -----------------------------
update_mu_sigma2_from_summaries <- function(Y_list, summaries, pars_old) {
  K <- pars_old$K
  d <- pars_old$d
  mu_new <- vector("list", K)
  sigma2_new <- numeric(K)

  for (k in seq_len(K)) {
    num <- rep(0, d)
    den <- 0
    for (i in seq_along(Y_list)) {
      y <- Y_list[[i]]
      summ <- summaries[[i]]
      R <- length(summ$post_w)
      T_len <- nrow(y)
      for (r in seq_len(R)) {
        wr <- summ$post_w[r]
        gk <- summ$gamma_nodes[[r]][, k]
        num <- num + wr * colSums(gk * (y - matrix(summ$f_mat[r, ], nrow = T_len, ncol = d, byrow = TRUE)))
        den <- den + wr * sum(gk)
      }
    }
    mu_new[[k]] <- num / max(den, 1e-12)
  }

  for (k in seq_len(K)) {
    numer <- 0
    denom <- 0
    mu_k <- mu_new[[k]]
    for (i in seq_along(Y_list)) {
      y <- Y_list[[i]]
      summ <- summaries[[i]]
      R <- length(summ$post_w)
      T_len <- nrow(y)
      for (r in seq_len(R)) {
        wr <- summ$post_w[r]
        gk <- summ$gamma_nodes[[r]][, k]
        resid <- y - matrix(mu_k + summ$f_mat[r, ], nrow = T_len, ncol = d, byrow = TRUE)
        numer <- numer + wr * sum(gk * rowSums(resid^2))
        denom <- denom + wr * d * sum(gk)
      }
    }
    sigma2_new[k] <- max(numer / max(denom, 1e-12), 1e-6)
  }

  list(mu = mu_new, sigma2 = sigma2_new)
}

update_tau2_from_summaries <- function(summaries, d) {
  n <- length(summaries)
  val <- 0
  for (i in seq_len(n)) {
    summ <- summaries[[i]]
    val <- val + sum(diag(summ$post_cov)) + sum(summ$post_mean^2)
  }
  max(val / (n * d), 1e-6)
}

# -----------------------------
# Generic exact/MC outer-integration EM fit
# -----------------------------
fit_outer_integration_em <- function(dat, init_pars = NULL,
                                     integration = c("gh", "mc"),
                                     n_nodes_1d = NULL,
                                     mc_points = NULL,
                                     max_iter = 60,
                                     tol = 1e-4,
                                     verbose = FALSE,
                                     seed = NULL) {
  integration <- match.arg(integration)
  if (!is.null(seed)) set.seed(seed)

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

  if (integration == "gh") {
    if (is.null(n_nodes_1d)) stop("For Gaussian quadrature, provide n_nodes_1d.")
    base_rule <- tensor_rule_std_normal(d = d, n_nodes_1d = n_nodes_1d)
    method_name <- sprintf("QEM-%d", n_nodes_1d)
    total_nodes <- nrow(base_rule$z)
  } else {
    if (is.null(mc_points)) stop("For MCEM, provide mc_points.")
    method_name <- sprintf("MCEM-%d", mc_points)
    total_nodes <- mc_points
  }

  prev_metric <- Inf
  iter_used <- 0
  gamma_list <- xi_list <- m_list <- V_list <- state_hat <- NULL

  for (iter in seq_len(max_iter)) {
    if (integration == "mc") {
      z_draws <- matrix(rnorm(mc_points * d), nrow = mc_points, ncol = d)
      base_w <- rep(1 / mc_points, mc_points)
    } else {
      z_draws <- base_rule$z
      base_w <- base_rule$w
    }

    summaries <- vector("list", n)
    gamma_list <- vector("list", n)
    xi_list <- vector("list", n)
    m_list <- vector("list", n)
    V_list <- vector("list", n)
    obs_logLik <- 0

    for (i in seq_len(n)) {
      summ <- subject_posterior_summary(Y_list[[i]], pars, z_mat = z_draws, w = base_w, method_label = method_name)
      summaries[[i]] <- summ
      gamma_list[[i]] <- summ$gamma_bar
      xi_list[[i]] <- summ$xi_bar
      m_list[[i]] <- summ$post_mean
      V_list[[i]] <- summ$post_cov
      obs_logLik <- obs_logLik + summ$logLik
    }

    pars_old <- copy_pars(pars)
    pars$delta <- update_delta(gamma_list)
    pars$Gamma <- update_Gamma(gamma_list, xi_list)
    upd <- update_mu_sigma2_from_summaries(Y_list, summaries, pars_old)
    pars$mu <- upd$mu
    pars$sigma2 <- upd$sigma2
    pars$tau2 <- update_tau2_from_summaries(summaries, d)

    dist <- parameter_distance(pars, pars_old)
    iter_used <- iter
    if (verbose) cat(sprintf("%s iter %d, obs_ll=%.4f, dist=%.6f, tau2=%.4f\n", method_name, iter, obs_logLik, dist, pars$tau2))
    if (iter > 1 && dist < tol) break
    prev_metric <- dist
  }

  state_hat <- vector("list", n)
  for (i in seq_len(n)) state_hat[[i]] <- viterbi_conditional(Y_list[[i]], m_list[[i]], pars)

  list(
    method = method_name,
    pars = pars,
    m_list = m_list,
    V_list = V_list,
    gamma_list = gamma_list,
    xi_list = xi_list,
    state_hat = state_hat,
    iterations = iter_used,
    total_nodes = total_nodes,
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

  mu_rmse <- sqrt(mean(unlist(lapply(seq_len(K), function(k) (est$mu[[k]] - truth$mu[[k]])^2))))
  sigma2_rmse <- sqrt(mean((est$sigma2 - truth$sigma2)^2))
  tau2_abs <- abs(est$tau2 - truth$tau2)
  gamma_abs <- mean(abs(est$Gamma - truth$Gamma))

  f_mse <- mean(sapply(seq_along(dat$F), function(i) mean((fit$m_list[[i]] - dat$F[[i]])^2)))
  f_cor_mean <- mean(sapply(seq_along(dat$F), function(i) safe_cor(fit$m_list[[i]], dat$F[[i]])), na.rm = TRUE)
  state_acc <- mean(unlist(lapply(seq_along(dat$U), function(i) fit$state_hat[[i]] == dat$U[[i]])))

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
  cat(sprintf("  %-8s | K=%d | iter=%3d | time=%8.3f sec | nodes=%6d | tau2=%.4f\n",
              fit$method, K, fit$iterations, fit$time_sec, fit$total_nodes, fit$pars$tau2))
  cat("    mean(mu_k):", paste(round(mu_means, 3), collapse = ", "), "\n")
}

# -----------------------------
# One replication under a common initialization
# -----------------------------
run_one_replication_compare <- function(K = 3, n = 60, T_len = 60, d = 2, tau2 = 1,
                                        max_iter = 60, tol = 1e-4, seed = 1,
                                        gh_nodes = c(3, 5),
                                        mc_points = c(50, 100),
                                        verbose = FALSE) {
  truth <- make_true_pars(K = K, d = d, tau2 = tau2)
  dat <- simulate_gaussian_mhmm(n = n, T_len = T_len, pars = truth, seed = seed)

  init_pars <- make_true_pars(K = K, d = d, tau2 = max(0.4, tau2 * 0.7))
  init_centers <- seq(from = 0.8, to = -0.8, length.out = K)
  init_pars$mu <- lapply(init_centers, function(a) rep(a, d))
  init_pars$sigma2 <- rep(1.2, K)
  offdiag <- 0.15 / max(K - 1, 1)
  init_pars$Gamma <- matrix(offdiag, K, K)
  diag(init_pars$Gamma) <- 0.85
  init_pars$Gamma <- init_pars$Gamma / rowSums(init_pars$Gamma)
  init_pars$delta <- stationary_dist(init_pars$Gamma)

  fits <- list()
  results <- list()
  idx <- 1

  fit_avem <- fit_mean_anchor_avem(dat, init_pars = init_pars, max_iter = max_iter, tol = tol, verbose = verbose)
  fit_avem <- match_labels_general(fit_avem, truth)
  fits[[idx]] <- fit_avem
  results[[idx]] <- evaluate_fit(fit_avem, truth = truth, dat = dat)
  idx <- idx + 1

  for (q in gh_nodes) {
    fit_q <- fit_outer_integration_em(dat, init_pars = init_pars, integration = "gh",
                                      n_nodes_1d = q, max_iter = max_iter, tol = tol,
                                      verbose = verbose, seed = seed + q)
    fit_q <- match_labels_general(fit_q, truth)
    fits[[idx]] <- fit_q
    results[[idx]] <- evaluate_fit(fit_q, truth = truth, dat = dat)
    idx <- idx + 1
  }

  for (m in mc_points) {
    fit_mc <- fit_outer_integration_em(dat, init_pars = init_pars, integration = "mc",
                                       mc_points = m, max_iter = max_iter, tol = tol,
                                       verbose = verbose, seed = seed + m)
    fit_mc <- match_labels_general(fit_mc, truth)
    fits[[idx]] <- fit_mc
    results[[idx]] <- evaluate_fit(fit_mc, truth = truth, dat = dat)
    idx <- idx + 1
  }

  list(data = dat, truth = truth, fits = fits, results = do.call(rbind, results))
}

# -----------------------------
# Experiment runner for the requested comparison
# -----------------------------
run_requested_comparison <- function(n_rep = 10,
                                     d_values = c(2, 3),
                                     K = 3,
                                     n = 60,
                                     T_len = 60,
                                     tau2 = 1,
                                     gh_nodes = c(3, 5),
                                     mc_points = c(50, 100),
                                     max_iter = 60,
                                     tol = 1e-4,
                                     base_seed = 2026,
                                     verbose = TRUE) {
  out <- list()
  fail_log <- list()
  idx <- 1
  fail_idx <- 1
  
  for (d in d_values) {
    for (r in seq_len(n_rep)) {
      if (verbose) {
        cat(sprintf("d=%d | rep=%d/%d | K=%d n=%d T=%d\n", d, r, n_rep, K, n, T_len))
      }
      
      seed_now <- base_seed + 1000 * d + r
      
      one <- tryCatch(
        {
          run_one_replication_compare(
            K = K, n = n, T_len = T_len, d = d, tau2 = tau2,
            max_iter = max_iter, tol = tol,
            seed = seed_now,
            gh_nodes = gh_nodes,
            mc_points = mc_points,
            verbose = FALSE
          )
        },
        error = function(e) {
          fail_log[[fail_idx]] <<- data.frame(
            d = d,
            rep = r,
            seed = seed_now,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
          fail_idx <<- fail_idx + 1
          
          if (verbose) {
            cat(sprintf("  failed at d=%d, rep=%d, seed=%d\n", d, r, seed_now))
            cat("  reason:", conditionMessage(e), "\n")
          }
          NULL
        }
      )
      
      if (!is.null(one)) {
        out[[idx]] <- one$results
        idx <- idx + 1
      }
    }
  }
  
  res <- if (length(out) > 0) do.call(rbind, out) else NULL
  fail_df <- if (length(fail_log) > 0) do.call(rbind, fail_log) else NULL
  
  list(results = res, failures = fail_df)
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
  out[order(out$d, out$method), ]
}

# -----------------------------
# Example usage for the requested study
# -----------------------------
# results_req <- run_requested_comparison(
#   n_rep = 10,
#   d_values = c(2, 3),
#   K = 3,
#   n = 60,
#   T_len = 60,
#   tau2 = 1,
#   gh_nodes = c(3, 5),
#   mc_points = c(50, 100),
#   max_iter = 60,
#   tol = 1e-4,
#   base_seed = 2026,
#   verbose = TRUE
# )
#
# summary_req <- aggregate_results(results_req)
# print(summary_req)

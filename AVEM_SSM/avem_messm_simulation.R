# ============================================================
# AVEM simulation for mixed-effects state-space models (MESSM)
# ------------------------------------------------------------
# This script implements:
#   (i)   data generation
#   (ii)  AVEM algorithm following the LaTeX description
#   (iii) summary metrics and runtime reporting
#
# The code is written to be self-contained in base R.
# You can change the hyperparameters in the CONFIG section.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
config <- list(
  seed = 123,
  n = 50,
  T = 100,
  p = 4,
  q = 2,
  max_iter = 200,
  rel_tol = 1e-4,
  verbose = TRUE,
  save_prefix = "avem_messm_run",
  # DGP hyperparameters
  tau_g2 = 0.05,
  tau_h2 = 0.05,
  r_diag = rep(0.25, 4),   # diagonal entries of R
  # Initialization noise levels
  init_mu_g_sd = 0.10,
  init_mu_h_sd = 0.10,
  init_logR_sd = 0.10
)

# -----------------------------
# BASIC LINEAR-ALGEBRA UTILITIES
# -----------------------------
symmetrize <- function(A) {
  0.5 * (A + t(A))
}

safe_chol <- function(A, base_jitter = 1e-8, max_tries = 8) {
  A <- symmetrize(A)
  p <- nrow(A)
  for (k in 0:max_tries) {
    jitter <- base_jitter * (10^k)
    out <- try(chol(A + diag(jitter, p)), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  stop("Cholesky decomposition failed even after adding jitter.")
}

safe_solve <- function(A, b = NULL, base_jitter = 1e-8, max_tries = 8) {
  A <- symmetrize(A)
  p <- nrow(A)
  for (k in 0:max_tries) {
    jitter <- base_jitter * (10^k)
    Achol <- try(chol(A + diag(jitter, p)), silent = TRUE)
    if (!inherits(Achol, "try-error")) {
      if (is.null(b)) {
        return(chol2inv(Achol))
      } else {
        return(backsolve(Achol, forwardsolve(t(Achol), b)))
      }
    }
  }
  stop("Linear solve failed even after adding jitter.")
}

logdet_spd <- function(A) {
  C <- safe_chol(A)
  2 * sum(log(diag(C)))
}

rmvnorm1 <- function(mu, Sigma) {
  L <- safe_chol(Sigma)
  as.numeric(mu + t(L) %*% rnorm(length(mu)))
}

vec <- function(A) as.vector(A)

matrix_from_vec <- function(x, nrow_, ncol_) {
  matrix(x, nrow = nrow_, ncol = ncol_)
}

# -----------------------------
# LOWER-TRIANGULAR LOADING MAPS
# -----------------------------
make_lower_maps <- function(p, q) {
  # Free entries are H[r, c] for r >= c
  positions <- list()
  idx <- 1
  for (c in 1:q) {
    for (r in c:p) {
      positions[[idx]] <- c(r = r, c = c)
      idx <- idx + 1
    }
  }
  d_h <- length(positions)

  S_H <- matrix(0, nrow = p * q, ncol = d_h)
  row_select <- vector("list", p)
  for (j in 1:p) row_select[[j]] <- matrix(0, nrow = q, ncol = d_h)

  for (ell in 1:d_h) {
    r <- positions[[ell]]["r"]
    c <- positions[[ell]]["c"]
    E <- matrix(0, nrow = p, ncol = q)
    E[r, c] <- 1
    S_H[, ell] <- vec(E)
    row_select[[r]][c, ell] <- 1
  }


  list(
    positions = positions,
    d_h = d_h,
    S_H = S_H,
    row_select = row_select
  )
}

vecl_lower <- function(H, maps) {
  out <- numeric(maps$d_h)
  for (ell in 1:maps$d_h) {
    out[ell] <- H[maps$positions[[ell]]["r"], maps$positions[[ell]]["c"]]
  }
  out
}

unvecl_lower <- function(h, p, q, maps) {
  H <- matrix(0, nrow = p, ncol = q)
  for (ell in 1:maps$d_h) {
    H[maps$positions[[ell]]["r"], maps$positions[[ell]]["c"]] <- h[ell]
  }
  H
}

# -----------------------------
# MODEL-SPECIFIC HELPERS
# -----------------------------
stabilize_G <- function(G, max_radius = 0.98) {
  eig <- eigen(G, only.values = TRUE)$values
  rad <- max(Mod(eig))
  if (rad < max_radius) return(G)
  G * (max_radius / rad)
}

make_default_truth <- function(p, q, tau_g2, tau_h2, r_diag, maps) {
  if (!(p >= q)) stop("Need p >= q for lower-triangular H.")

  if (q == 2) {
    G <- matrix(c(0.70, -0.10,
                  0.10,  0.60), nrow = 2, byrow = TRUE)
  } else {
    G <- diag(seq(0.55, 0.75, length.out = q))
    if (q >= 2) {
      G[1, 2] <- 0.08
      G[2, 1] <- -0.05
    }
  }
  G <- stabilize_G(G)

  H <- matrix(0, nrow = p, ncol = q)
  for (c in 1:q) {
    for (r in c:p) {
      if (r == c) {
        H[r, c] <- 0.9 + 0.1 * (q - c)
      } else {
        H[r, c] <- 0.2 + 0.1 * ((r + c) %% 3)
      }
    }
  }

  mu_g <- vec(G)
  mu_h <- vecl_lower(H, maps)
  d_g <- q * q
  d_h <- maps$d_h

  list(
    m0 = rep(0, q),
    P0 = diag(q),
    Q = diag(q),  # fixed for identifiability
    R = diag(r_diag[1:p]),
    mu_g = mu_g,
    mu_h = mu_h,
    Sigma_g = tau_g2 * diag(d_g),
    Sigma_h = tau_h2 * diag(d_h),
    G = G,
    H = H
  )
}

sample_subject_parameters <- function(truth, p, q, maps) {
  g_i <- rmvnorm1(truth$mu_g, truth$Sigma_g)
  G_i <- matrix_from_vec(g_i, q, q)
  G_i <- stabilize_G(G_i)
  g_i <- vec(G_i)

  h_i <- rmvnorm1(truth$mu_h, truth$Sigma_h)
  H_i <- unvecl_lower(h_i, p, q, maps)

  list(g = g_i, h = h_i, G = G_i, H = H_i)
}

# -----------------------------
# DATA GENERATION
# -----------------------------
simulate_messm_data <- function(config, truth, maps) {
  set.seed(config$seed)
  n <- config$n
  T_ <- config$T
  p <- config$p
  q <- config$q

  data_list <- vector("list", n)
  latent_list <- vector("list", n)
  subj_pars <- vector("list", n)

  for (i in 1:n) {
    pars_i <- sample_subject_parameters(truth, p, q, maps)
    G_i <- pars_i$G
    H_i <- pars_i$H

    U <- matrix(0, nrow = q, ncol = T_)
    D <- matrix(0, nrow = p, ncol = T_)

    U[, 1] <- rmvnorm1(truth$m0, truth$P0)
    D[, 1] <- as.numeric(H_i %*% U[, 1] + rmvnorm1(rep(0, p), truth$R))

    if (T_ >= 2) {
      for (t in 2:T_) {
        U[, t] <- as.numeric(G_i %*% U[, t - 1] + rnorm(q))
        D[, t] <- as.numeric(H_i %*% U[, t] + rmvnorm1(rep(0, p), truth$R))
      }
    }

    data_list[[i]] <- D
    latent_list[[i]] <- U
    subj_pars[[i]] <- pars_i
  }

  list(
    D = data_list,
    U = latent_list,
    subject_pars = subj_pars,
    truth = truth
  )
}

# -----------------------------
# KALMAN FILTER / RTS SMOOTHER
# -----------------------------
kalman_smoother_subject <- function(Y, G, H, m0, P0, R, Q = NULL) {
  if (is.null(Q)) Q <- diag(nrow(G))
  p <- nrow(Y)
  T_ <- ncol(Y)
  q <- nrow(G)

  m_pred <- matrix(0, nrow = q, ncol = T_)
  m_filt <- matrix(0, nrow = q, ncol = T_)
  m_smooth <- matrix(0, nrow = q, ncol = T_)

  P_pred <- array(0, dim = c(q, q, T_))
  P_filt <- array(0, dim = c(q, q, T_))
  P_smooth <- array(0, dim = c(q, q, T_))
  J_arr <- array(0, dim = c(q, q, max(T_ - 1, 1)))

  loglik <- 0

  # t = 1
  m_pred[, 1] <- m0
  P_pred[, , 1] <- P0
  S1 <- H %*% P_pred[, , 1] %*% t(H) + R
  S1_inv <- safe_solve(S1)
  K1 <- P_pred[, , 1] %*% t(H) %*% S1_inv
  innov1 <- Y[, 1] - H %*% m_pred[, 1]
  m_filt[, 1] <- as.numeric(m_pred[, 1] + K1 %*% innov1)
  P_filt[, , 1] <- symmetrize(P_pred[, , 1] - K1 %*% H %*% P_pred[, , 1])
  loglik <- loglik - 0.5 * (p * log(2 * pi) + logdet_spd(S1) + as.numeric(t(innov1) %*% S1_inv %*% innov1))

  # t >= 2
  if (T_ >= 2) {
    for (t in 2:T_) {
      m_pred[, t] <- as.numeric(G %*% m_filt[, t - 1])
      P_pred[, , t] <- symmetrize(G %*% P_filt[, , t - 1] %*% t(G) + Q)

      S_t <- H %*% P_pred[, , t] %*% t(H) + R
      S_t_inv <- safe_solve(S_t)
      K_t <- P_pred[, , t] %*% t(H) %*% S_t_inv
      innov_t <- Y[, t] - H %*% m_pred[, t]

      m_filt[, t] <- as.numeric(m_pred[, t] + K_t %*% innov_t)
      P_filt[, , t] <- symmetrize(P_pred[, , t] - K_t %*% H %*% P_pred[, , t])

      loglik <- loglik - 0.5 * (p * log(2 * pi) + logdet_spd(S_t) + as.numeric(t(innov_t) %*% S_t_inv %*% innov_t))
    }
  }

  # RTS smoother
  m_smooth[, T_] <- m_filt[, T_]
  P_smooth[, , T_] <- P_filt[, , T_]

  if (T_ >= 2) {
    for (t in (T_ - 1):1) {
      J_t <- P_filt[, , t] %*% t(G) %*% safe_solve(P_pred[, , t + 1])
      J_arr[, , t] <- J_t
      m_smooth[, t] <- as.numeric(m_filt[, t] + J_t %*% (m_smooth[, t + 1] - m_pred[, t + 1]))
      P_smooth[, , t] <- symmetrize(P_filt[, , t] + J_t %*% (P_smooth[, , t + 1] - P_pred[, , t + 1]) %*% t(J_t))
    }
  }

  lag_cov <- array(0, dim = c(q, q, T_))
  if (T_ >= 2) {
    for (t in 2:T_) {
      lag_cov[, , t] <- P_smooth[, , t] %*% t(J_arr[, , t - 1])
    }
  }

  list(
    m_pred = m_pred,
    P_pred = P_pred,
    m_filt = m_filt,
    P_filt = P_filt,
    m_smooth = m_smooth,
    P_smooth = P_smooth,
    lag_cov = lag_cov,
    J = J_arr,
    loglik = loglik
  )
}

# -----------------------------
# INITIALIZATION
# -----------------------------
initialize_avem <- function(sim_data, config, truth, maps) {
  n <- config$n
  p <- config$p
  q <- config$q

  d_g <- q * q
  d_h <- maps$d_h

  set.seed(config$seed + 1)

  mu_g <- truth$mu_g + rnorm(d_g, sd = config$init_mu_g_sd)
  mu_h <- truth$mu_h + rnorm(d_h, sd = config$init_mu_h_sd)
  Sigma_g <- truth$Sigma_g
  Sigma_h <- truth$Sigma_h

  m0 <- truth$m0
  P0 <- truth$P0
  R_diag <- diag(truth$R) * exp(rnorm(p, sd = config$init_logR_sd))
  R <- diag(pmax(R_diag, 1e-3))

  nu_g <- matrix(rep(mu_g, n), nrow = d_g, ncol = n)
  nu_h <- matrix(rep(mu_h, n), nrow = d_h, ncol = n)
  Omega_g <- replicate(n, Sigma_g, simplify = FALSE)
  Omega_h <- replicate(n, Sigma_h, simplify = FALSE)

  list(
    theta = list(
      mu_g = mu_g,
      mu_h = mu_h,
      Sigma_g = Sigma_g,
      Sigma_h = Sigma_h,
      m0 = m0,
      P0 = P0,
      R = R
    ),
    variational = list(
      nu_g = nu_g,
      nu_h = nu_h,
      Omega_g = Omega_g,
      Omega_h = Omega_h
    )
  )
}



# -----------------------------
# EXPECTATION HELPERS
# -----------------------------
compute_Hrow_second_moment <- function(M_h, row_select_j) {
  row_select_j %*% M_h %*% t(row_select_j)
}

compute_R_update <- function(sim_data, state_stats, variational, config, maps) {
  n <- config$n
  T_ <- config$T
  p <- config$p

  r_new <- numeric(p)

  for (j in 1:p) {
    accum <- 0
    Rj_sel <- maps$row_select[[j]]
    for (i in 1:n) {
      nu_h_i <- variational$nu_h[, i]
      Omega_h_i <- variational$Omega_h[[i]]
      M_h_i <- Omega_h_i + tcrossprod(nu_h_i)
      Hbar_i <- unvecl_lower(nu_h_i, config$p, config$q, maps)

      for (t in 1:T_) {
        yjt <- sim_data$D[[i]][j, t]
        m_t <- state_stats[[i]]$m_smooth[, t]
        Q_t <- state_stats[[i]]$Q_t[, , t]

        EHjHj <- compute_Hrow_second_moment(M_h_i, Rj_sel)
        accum <- accum + yjt^2 - 2 * yjt * sum(Hbar_i[j, ] * m_t) + sum(Q_t * EHjHj)
      }
    }
    r_new[j] <- accum / (n * T_)
  }

  diag(pmax(r_new, 1e-6))
}

compute_transition_quad <- function(Q_t, Q_lag, Q_prev, nu_g, Omega_g, q) {
  M_g <- Omega_g + tcrossprod(nu_g)
  term1 <- sum(diag(Q_t))
  term2 <- -2 * sum(vec(Q_lag) * nu_g)
  term3 <- sum((kronecker(Q_prev, diag(q))) * M_g)
  term1 + term2 + term3
}

compute_obs_quad <- function(y_t, Q_t, m_t, nu_h, Omega_h, R, maps, p, q) {
  Hbar <- unvecl_lower(nu_h, p, q, maps)
  M_h <- Omega_h + tcrossprod(nu_h)

  out <- 0
  for (j in 1:p) {
    Rj_sel <- maps$row_select[[j]]
    EHjHj <- compute_Hrow_second_moment(M_h, Rj_sel)
    out <- out + (1 / R[j, j]) * (
      y_t[j]^2 - 2 * y_t[j] * sum(Hbar[j, ] * m_t) + sum(Q_t * EHjHj)
    )
  }
  out
}

compute_anchor_obs_quad <- function(y_t, Q_t, m_t, H0, R) {
  R_inv <- safe_solve(R)
  resid_mean <- y_t - H0 %*% m_t
  P_t <- Q_t - tcrossprod(m_t)
  term1 <- as.numeric(t(resid_mean) %*% R_inv %*% resid_mean)
  term2 <- sum(diag(R_inv %*% H0 %*% P_t %*% t(H0)))
  term1 + term2
}

compute_anchor_trans_quad <- function(Q_t, Q_lag, Q_prev, G0) {
  term1 <- sum(diag(Q_t))
  term2 <- -2 * sum(vec(Q_lag) * vec(G0))
  term3 <- sum((kronecker(Q_prev, diag(nrow(G0)))) * tcrossprod(vec(G0)))
  term1 + term2 + term3
}

gaussian_kl <- function(mu_q, Sigma_q, mu_p, Sigma_p) {
  d <- length(mu_q)
  Sigma_p_inv <- safe_solve(Sigma_p)
  diff <- mu_q - mu_p
  0.5 * (
    logdet_spd(Sigma_p) - logdet_spd(Sigma_q) - d +
      sum(diag(Sigma_p_inv %*% Sigma_q)) +
      as.numeric(t(diff) %*% Sigma_p_inv %*% diff)
  )
}

compute_subject_elbo <- function(i, sim_data, state_stats_i, variational, theta, config, maps) {
  T_ <- config$T
  p <- config$p
  q <- config$q

  nu_g_i <- variational$nu_g[, i]
  nu_h_i <- variational$nu_h[, i]
  Omega_g_i <- variational$Omega_g[[i]]
  Omega_h_i <- variational$Omega_h[[i]]

  Y <- sim_data$D[[i]]
  Q_arr <- state_stats_i$Q_t
  Qlag_arr <- state_stats_i$Q_lag
  m_smooth <- state_stats_i$m_smooth

  # E_q,p0 [log p(f_i)] - E_q[log q_i]
  prior_minus_entropy <- -gaussian_kl(nu_g_i, Omega_g_i, theta$mu_g, theta$Sigma_g) -
    gaussian_kl(nu_h_i, Omega_h_i, theta$mu_h, theta$Sigma_h)

  # Initial-state contribution
  m1 <- m_smooth[, 1]
  Q1 <- Q_arr[, , 1]
  P1 <- Q1 - tcrossprod(m1)
  P0_inv <- safe_solve(theta$P0)
  init_quad <- sum(diag(P0_inv %*% (P1 + tcrossprod(m1 - theta$m0))))
  init_term <- -0.5 * (logdet_spd(theta$P0) + init_quad)

  # Observation and transition contributions under q and p0
  obs_term <- 0
  trans_term <- 0
  for (t in 1:T_) {
    obs_term <- obs_term - 0.5 * (logdet_spd(theta$R) + compute_obs_quad(
      Y[, t], Q_arr[, , t], m_smooth[, t],
      nu_h_i, Omega_h_i, theta$R, maps, p, q
    ))
    if (t >= 2) {
      trans_term <- trans_term - 0.5 * compute_transition_quad(
        Q_arr[, , t], Qlag_arr[, , t], Q_arr[, , t - 1],
        nu_g_i, Omega_g_i, q
      )
    }
  }

  # Anchor correction:
  # E_0[log p(U,D | f0)] - log p(D | f0)
  G0 <- matrix_from_vec(state_stats_i$anchor_g, q, q)
  H0 <- unvecl_lower(state_stats_i$anchor_h, p, q, maps)

  anchor_init_quad <- sum(diag(P0_inv %*% (P1 + tcrossprod(m1 - theta$m0))))
  anchor_part <- -0.5 * (logdet_spd(theta$P0) + anchor_init_quad)

  for (t in 1:T_) {
    anchor_part <- anchor_part - 0.5 * (logdet_spd(theta$R) +
      compute_anchor_obs_quad(Y[, t], Q_arr[, , t], m_smooth[, t], H0, theta$R))
    if (t >= 2) {
      anchor_part <- anchor_part - 0.5 * compute_anchor_trans_quad(
        Q_arr[, , t], Qlag_arr[, , t], Q_arr[, , t - 1], G0
      )
    }
  }

  prior_minus_entropy + init_term + obs_term + trans_term - anchor_part + state_stats_i$loglik
}

# -----------------------------
# AVEM ALGORITHM
# -----------------------------
run_avem_messm <- function(sim_data, config, maps, init = NULL) {
  if (is.null(init)) {
    init <- initialize_avem(sim_data, config, sim_data$truth, maps)
  }

  theta <- init$theta
  variational <- init$variational
  n <- config$n
  T_ <- config$T
  p <- config$p
  q <- config$q

  history <- list(elbo = numeric(0), rel_change = numeric(0))
  start_time <- proc.time()[3]

  for (iter in 1:config$max_iter) {
    theta_old <- theta
    state_stats <- vector("list", n)

    # E-step
    for (i in 1:n) {
      g0_i <- variational$nu_g[, i]
      h0_i <- variational$nu_h[, i]

      G0_i <- matrix_from_vec(g0_i, q, q)
      H0_i <- unvecl_lower(h0_i, p, q, maps)

      ks <- kalman_smoother_subject(
        Y = sim_data$D[[i]],
        G = G0_i,
        H = H0_i,
        m0 = theta$m0,
        P0 = theta$P0,
        R = theta$R,
        Q = diag(q)
      )

      Q_t <- array(0, dim = c(q, q, T_))
      Q_lag <- array(0, dim = c(q, q, T_))
      for (t in 1:T_) {
        Q_t[, , t] <- symmetrize(ks$P_smooth[, , t] + tcrossprod(ks$m_smooth[, t]))
      }
      if (T_ >= 2) {
        for (t in 2:T_) {
          Q_lag[, , t] <- ks$lag_cov[, , t] + tcrossprod(ks$m_smooth[, t], ks$m_smooth[, t - 1])
        }
      }

      state_stats[[i]] <- list(
        m_smooth = ks$m_smooth,
        P_smooth = ks$P_smooth,
        lag_cov = ks$lag_cov,
        Q_t = Q_t,
        Q_lag = Q_lag,
        loglik = ks$loglik,
        anchor_g = g0_i,
        anchor_h = h0_i
      )

      # Update q_i(g_i)
      Lambda_g_i <- safe_solve(theta$Sigma_g)
      eta_g_i <- safe_solve(theta$Sigma_g, theta$mu_g)
      if (T_ >= 2) {
        for (t in 2:T_) {
          Lambda_g_i <- Lambda_g_i + kronecker(Q_t[, , t - 1], diag(q))
          eta_g_i <- eta_g_i + vec(Q_lag[, , t])
        }
      }
      Omega_g_i <- safe_solve(Lambda_g_i)
      nu_g_i <- as.numeric(Omega_g_i %*% eta_g_i)

      # Update q_i(h_i)
      R_inv <- safe_solve(theta$R)
      Lambda_h_i <- safe_solve(theta$Sigma_h)
      eta_h_i <- safe_solve(theta$Sigma_h, theta$mu_h)

      for (t in 1:T_) {
        Lambda_h_i <- Lambda_h_i + t(maps$S_H) %*% kronecker(Q_t[, , t], R_inv) %*% maps$S_H
        eta_h_i <- eta_h_i + t(maps$S_H) %*% kronecker(ks$m_smooth[, t], as.numeric(R_inv %*% sim_data$D[[i]][, t]))
      }

      Omega_h_i <- safe_solve(Lambda_h_i)
      nu_h_i <- as.numeric(Omega_h_i %*% eta_h_i)

      variational$nu_g[, i] <- nu_g_i
      variational$Omega_g[[i]] <- symmetrize(Omega_g_i)
      variational$nu_h[, i] <- nu_h_i
      variational$Omega_h[[i]] <- symmetrize(Omega_h_i)
    }

    # M-step
    theta$mu_g <- rowMeans(variational$nu_g)
    theta$mu_h <- rowMeans(variational$nu_h)

    d_g <- length(theta$mu_g)
    d_h <- length(theta$mu_h)
    Sigma_g_new <- matrix(0, nrow = d_g, ncol = d_g)
    Sigma_h_new <- matrix(0, nrow = d_h, ncol = d_h)
    for (i in 1:n) {
      dg_i <- variational$nu_g[, i] - theta$mu_g
      dh_i <- variational$nu_h[, i] - theta$mu_h
      Sigma_g_new <- Sigma_g_new + variational$Omega_g[[i]] + tcrossprod(dg_i)
      Sigma_h_new <- Sigma_h_new + variational$Omega_h[[i]] + tcrossprod(dh_i)
    }
    theta$Sigma_g <- symmetrize(Sigma_g_new / n)
    theta$Sigma_h <- symmetrize(Sigma_h_new / n)

    # Update m0 and P0
    m0_new <- Reduce("+", lapply(state_stats, function(z) z$m_smooth[, 1])) / n
    P0_new <- matrix(0, nrow = q, ncol = q)
    for (i in 1:n) {
      diff_i <- state_stats[[i]]$m_smooth[, 1] - m0_new
      P0_new <- P0_new + state_stats[[i]]$P_smooth[, , 1] + tcrossprod(diff_i)
    }
    theta$m0 <- as.numeric(m0_new)
    theta$P0 <- symmetrize(P0_new / n)

    # Update R
    theta$R <- compute_R_update(sim_data, state_stats, variational, config, maps)

    # Monitored anchored ELBO
    elbo_iter <- 0
    for (i in 1:n) {
      elbo_iter <- elbo_iter + compute_subject_elbo(
        i = i,
        sim_data = sim_data,
        state_stats_i = state_stats[[i]],
        variational = variational,
        theta = theta,
        config = config,
        maps = maps
      )
    }

    # Relative parameter change
    num <- sqrt(sum((theta$mu_g - theta_old$mu_g)^2) +
                  sum((theta$mu_h - theta_old$mu_h)^2) +
                  sum((theta$m0 - theta_old$m0)^2) +
                  sum((diag(theta$R) - diag(theta_old$R))^2))
    den <- sqrt(sum(theta_old$mu_g^2) + sum(theta_old$mu_h^2) +
                  sum(theta_old$m0^2) + sum(diag(theta_old$R)^2)) + 1e-8
    rel_change <- num / den

    history$elbo <- c(history$elbo, elbo_iter)
    history$rel_change <- c(history$rel_change, rel_change)

    if (config$verbose) {
      cat(sprintf("Iter %3d | ELBO = % .6f | rel.change = %.6e\n", iter, elbo_iter, rel_change))
    }

    if (iter >= 2 && rel_change < config$rel_tol) break
  }

  total_time <- proc.time()[3] - start_time

  list(
    theta = theta,
    variational = variational,
    state_stats = state_stats,
    history = history,
    total_time = total_time,
    n_iter = length(history$elbo)
  )
}

# -----------------------------
# EVALUATION METRICS
# -----------------------------
compute_metrics <- function(fit, sim_data, config, maps) {
  n <- config$n
  p <- config$p
  q <- config$q
  truth <- sim_data$truth

  # Population-level metrics
  rmse_mu_g <- sqrt(mean((fit$theta$mu_g - truth$mu_g)^2))
  rmse_mu_h <- sqrt(mean((fit$theta$mu_h - truth$mu_h)^2))
  mse_Sigma_g <- mean((fit$theta$Sigma_g - truth$Sigma_g)^2)
  mse_Sigma_h <- mean((fit$theta$Sigma_h - truth$Sigma_h)^2)
  rmse_R <- sqrt(mean((diag(fit$theta$R) - diag(truth$R))^2))
  rmse_m0 <- sqrt(mean((fit$theta$m0 - truth$m0)^2))
  mse_P0 <- mean((fit$theta$P0 - truth$P0)^2)

  # Subject-specific random effects
  mse_g_subject <- 0
  mse_h_subject <- 0
  for (i in 1:n) {
    g_true <- sim_data$subject_pars[[i]]$g
    h_true <- sim_data$subject_pars[[i]]$h
    mse_g_subject <- mse_g_subject + mean((fit$variational$nu_g[, i] - g_true)^2)
    mse_h_subject <- mse_h_subject + mean((fit$variational$nu_h[, i] - h_true)^2)
  }
  mse_g_subject <- mse_g_subject / n
  mse_h_subject <- mse_h_subject / n

  # Latent-state recovery
  latent_mse <- 0
  for (i in 1:n) {
    U_true <- sim_data$U[[i]]
    U_hat <- fit$state_stats[[i]]$m_smooth
    latent_mse <- latent_mse + mean((U_hat - U_true)^2)
  }
  latent_mse <- latent_mse / n

  # Matrix-level summaries for the population mean G and H
  G_hat <- matrix_from_vec(fit$theta$mu_g, q, q)
  H_hat <- unvecl_lower(fit$theta$mu_h, p, q, maps)
  mean_G_mse <- mean((G_hat - truth$G)^2)
  mean_H_mse <- mean((H_hat - truth$H)^2)

  data.frame(
    rmse_mu_g = rmse_mu_g,
    rmse_mu_h = rmse_mu_h,
    mse_Sigma_g = mse_Sigma_g,
    mse_Sigma_h = mse_Sigma_h,
    rmse_R = rmse_R,
    rmse_m0 = rmse_m0,
    mse_P0 = mse_P0,
    mse_g_subject = mse_g_subject,
    mse_h_subject = mse_h_subject,
    latent_mse = latent_mse,
    mean_G_mse = mean_G_mse,
    mean_H_mse = mean_H_mse,
    total_time_sec = fit$total_time,
    n_iter = fit$n_iter
  )
}

# -----------------------------
# MAIN DRIVER
# -----------------------------
run_one_simulation <- function(config) {
  maps <- make_lower_maps(config$p, config$q)
  truth <- make_default_truth(
    p = config$p,
    q = config$q,
    tau_g2 = config$tau_g2,
    tau_h2 = config$tau_h2,
    r_diag = config$r_diag,
    maps = maps
  )

  sim_data <- simulate_messm_data(config, truth, maps)
  fit <- run_avem_messm(sim_data, config, maps)
  metrics <- compute_metrics(fit, sim_data, config, maps)

  out <- list(
    config = config,
    truth = truth,
    sim_data = sim_data,
    fit = fit,
    metrics = metrics
  )

  prefix <- config$save_prefix
  saveRDS(out, paste0(prefix, "_full_output.rds"))
  write.csv(metrics, paste0(prefix, "_metrics.csv"), row.names = FALSE)

  cat("\n================ Summary ================\n")
  print(metrics)
  cat(sprintf("\nSaved full output to: %s\n", paste0(prefix, "_full_output.rds")))
  cat(sprintf("Saved metrics to:     %s\n", paste0(prefix, "_metrics.csv")))
  cat("=========================================\n")

  invisible(out)
}

# -----------------------------
# RUN
# -----------------------------
#result <- run_one_simulation(config)

# ================================================================
# Bernoulli mixed HMM: Laplace-AVEM vs QEM vs MCEM
#
# Model:
#   U_it in {1, 2} follows a common 2-state Markov chain.
#   f_i ~ N(0, tau2).
#   Y_it | (U_it = k, f_i) ~ Bernoulli(p_ikt),
#   p_ikt = expit(beta_k + f_i).
#
# Here K = 2 and d = 1.
#
# Methods:
#   1) Laplace-AVEM:
#        - Anchor f_{0i} at current posterior mean nu_i.
#        - Run forward-backward conditional on f_{0i}.
#        - Update q_i(f_i) by Laplace approximation.
#
#   2) QEM:
#        - Approximate the outer integral over f_i by Gaussian-Hermite
#          quadrature under f_i ~ N(0, tau2).
#
#   3) MCEM:
#        - Approximate the outer integral over f_i by Monte Carlo samples
#          from f_i ~ N(0, tau2).
#
# Required package:
#   install.packages("statmod")
# ================================================================


# ================================================================
# Utilities
# ================================================================

expit <- function(x) {
  1 / (1 + exp(-x))
}

log1pexp <- function(x) {
  # Stable computation of log(1 + exp(x)).
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

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

copy_pars <- function(pars) {
  list(
    beta = as.numeric(pars$beta),
    tau2 = as.numeric(pars$tau2),
    Gamma = matrix(pars$Gamma, nrow = 2, ncol = 2),
    delta = as.numeric(pars$delta),
    K = 2
  )
}

parameter_distance <- function(pars_new, pars_old) {
  max(
    max(abs(pars_new$beta - pars_old$beta)),
    abs(pars_new$tau2 - pars_old$tau2),
    max(abs(pars_new$Gamma - pars_old$Gamma)),
    max(abs(pars_new$delta - pars_old$delta))
  )
}


# ================================================================
# True parameters and data generation
# ================================================================

make_true_pars_bernoulli <- function(
    beta = c(-1.0, 1.0),
    tau2 = 0.7,
    Gamma = NULL
) {
  K <- 2
  
  if (is.null(Gamma)) {
    Gamma <- matrix(
      c(
        0.92, 0.08,
        0.08, 0.92
      ),
      nrow = K,
      byrow = TRUE
    )
  } else {
    stopifnot(all(dim(Gamma) == c(K, K)))
    Gamma <- Gamma / rowSums(Gamma)
  }
  
  list(
    beta = as.numeric(beta),
    tau2 = tau2,
    Gamma = Gamma,
    delta = stationary_dist(Gamma),
    K = K
  )
}

simulate_bernoulli_mhmm <- function(n, T_len, pars, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  K <- pars$K
  beta <- pars$beta
  tau2 <- pars$tau2
  Gamma <- pars$Gamma
  delta <- pars$delta
  
  Y_list <- vector("list", n)
  U_list <- vector("list", n)
  F_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    f_i <- rnorm(1, mean = 0, sd = sqrt(tau2))
    
    U <- integer(T_len)
    Y <- integer(T_len)
    
    U[1] <- sample.int(K, size = 1, prob = delta)
    p1 <- expit(beta[U[1]] + f_i)
    Y[1] <- rbinom(1, size = 1, prob = p1)
    
    if (T_len >= 2) {
      for (t in 2:T_len) {
        U[t] <- sample.int(K, size = 1, prob = Gamma[U[t - 1], ])
        p_t <- expit(beta[U[t]] + f_i)
        Y[t] <- rbinom(1, size = 1, prob = p_t)
      }
    }
    
    Y_list[[i]] <- Y
    U_list[[i]] <- U
    F_list[[i]] <- f_i
  }
  
  list(
    Y = Y_list,
    U = U_list,
    F = F_list,
    pars = pars,
    n = n,
    T_len = T_len,
    K = K
  )
}


# ================================================================
# Bernoulli emission log-density
# ================================================================

log_emission_matrix_bernoulli <- function(y, f, pars) {
  T_len <- length(y)
  K <- pars$K
  beta <- pars$beta
  
  out <- matrix(0, nrow = T_len, ncol = K)
  
  for (k in seq_len(K)) {
    eta <- beta[k] + f
    out[, k] <- y * eta - log1pexp(eta)
  }
  
  out
}


# ================================================================
# Forward-backward conditional on fixed f
# ================================================================

fb_conditional_bernoulli <- function(y, f, pars) {
  Gamma <- pars$Gamma
  delta <- pars$delta
  
  T_len <- length(y)
  K <- length(delta)
  
  log_emit <- log_emission_matrix_bernoulli(y, f, pars)
  emit <- exp(log_emit)
  emit[emit < 1e-300] <- 1e-300
  
  alpha <- matrix(0, nrow = T_len, ncol = K)
  beta_fb <- matrix(0, nrow = T_len, ncol = K)
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
  
  beta_fb[T_len, ] <- 1
  
  if (T_len >= 2) {
    for (t in (T_len - 1):1) {
      beta_fb[t, ] <-
        as.numeric(Gamma %*% (emit[t + 1, ] * beta_fb[t + 1, ])) / cscale[t + 1]
    }
  }
  
  gamma <- alpha * beta_fb
  gamma <- gamma / rowSums(gamma)
  
  xi <- array(0, dim = c(max(T_len - 1, 1), K, K))
  
  if (T_len >= 2) {
    for (t in 1:(T_len - 1)) {
      numer <- outer(alpha[t, ], emit[t + 1, ] * beta_fb[t + 1, ]) * Gamma
      xi[t, , ] <- numer / sum(numer)
    }
  }
  
  list(
    logLik = sum(log(cscale)),
    gamma = gamma,
    xi = xi
  )
}


# ================================================================
# Viterbi conditional on fixed f
# ================================================================

viterbi_conditional_bernoulli <- function(y, f, pars) {
  T_len <- length(y)
  K <- pars$K
  
  log_emit <- log_emission_matrix_bernoulli(y, f, pars)
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
    for (t in (T_len - 1):1) {
      path[t] <- psi[t + 1, path[t + 1]]
    }
  }
  
  path
}


# ================================================================
# Updates for delta and Gamma
# ================================================================

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


# ================================================================
# Gaussian quadrature rule
# ================================================================

std_normal_gh_rule <- function(n_nodes) {
  if (length(n_nodes) != 1 ||
      !is.numeric(n_nodes) ||
      n_nodes < 1 ||
      n_nodes != as.integer(n_nodes)) {
    stop("n_nodes must be a positive integer.")
  }
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required. Please run install.packages('statmod').")
  }
  
  gh <- statmod::gauss.quad.prob(n_nodes, dist = "normal")
  
  list(
    nodes = gh$nodes,
    weights = gh$weights
  )
}


# ================================================================
# Laplace update for q_i(f_i)
# ================================================================
# Anchored log-density:
#
# ell_i^A(f)
# =
# - f^2 / (2 tau2)
# + sum_t sum_k zeta_ikt
#     { y_it (beta_k + f) - log(1 + exp(beta_k + f)) }
# + const.
#
# The mode nu_i solves ell'(f) = 0.
# The variance is Omega_i = {-ell''(nu_i)}^{-1}.
#
# In this logistic random-intercept case, ell_i^A(f) is concave
# because the prior term and Bernoulli log-likelihood terms are concave.
# ================================================================

laplace_update_f <- function(y, gamma_i, pars, lower = -5, upper = 5) {
  beta <- pars$beta
  tau2 <- pars$tau2
  K <- pars$K
  
  ell <- function(f) {
    val <- -0.5 * f^2 / tau2
    
    for (k in seq_len(K)) {
      eta <- beta[k] + f
      val <- val + sum(gamma_i[, k] * (y * eta - log1pexp(eta)))
    }
    
    val
  }
  
  grad <- function(f) {
    val <- -f / tau2
    
    for (k in seq_len(K)) {
      eta <- beta[k] + f
      p <- expit(eta)
      val <- val + sum(gamma_i[, k] * (y - p))
    }
    
    val
  }
  
  hess <- function(f) {
    val <- -1 / tau2
    
    for (k in seq_len(K)) {
      eta <- beta[k] + f
      p <- expit(eta)
      val <- val - sum(gamma_i[, k] * p * (1 - p))
    }
    
    val
  }
  
  opt <- optimize(
    f = function(x) -ell(x),
    interval = c(lower, upper)
  )
  
  nu <- opt$minimum
  H <- hess(nu)
  
  # Negative Hessian should be positive.
  Omega <- 1 / max(-H, 1e-6)
  
  list(
    mean = nu,
    var = Omega,
    mode = nu,
    hess = H,
    obj = ell(nu),
    grad = grad(nu)
  )
}


# ================================================================
# M-step for beta in Laplace-AVEM
# ================================================================
# We maximize, for each k,
#
# sum_i sum_t zeta_ikt E_{q_i}
# [ y_it (beta_k + f_i) - log{1 + exp(beta_k + f_i)} ].
#
# Since q_i(f_i) = N(nu_i, Omega_i), the expectation is approximated
# by one-dimensional Gaussian quadrature.
# ================================================================

update_beta_avem <- function(Y_list,
                             gamma_list,
                             m_list,
                             V_list,
                             beta_old,
                             n_quad = 15,
                             lower = -5,
                             upper = 5) {
  K <- length(beta_old)
  rule <- std_normal_gh_rule(n_quad)
  
  beta_new <- numeric(K)
  
  for (k in seq_len(K)) {
    obj_k <- function(beta_k) {
      val <- 0
      
      for (i in seq_along(Y_list)) {
        y <- Y_list[[i]]
        gk <- gamma_list[[i]][, k]
        
        nu_i <- m_list[[i]]
        sd_i <- sqrt(max(V_list[[i]], 1e-10))
        
        f_nodes <- nu_i + sd_i * rule$nodes
        
        for (r in seq_along(f_nodes)) {
          eta <- beta_k + f_nodes[r]
          val <- val + rule$weights[r] * sum(gk * (y * eta - log1pexp(eta)))
        }
      }
      
      val
    }
    
    opt <- optimize(
      f = function(b) -obj_k(b),
      interval = c(lower, upper)
    )
    
    beta_new[k] <- opt$minimum
  }
  
  beta_new
}

update_tau2_avem <- function(m_list, V_list) {
  n <- length(m_list)
  
  val <- 0
  for (i in seq_len(n)) {
    val <- val + m_list[[i]]^2 + V_list[[i]]
  }
  
  max(val / n, 1e-6)
}


# ================================================================
# Laplace-AVEM fitting function
# ================================================================

fit_laplace_avem_bernoulli <- function(dat,
                                       init_pars = NULL,
                                       max_iter = 100,
                                       tol = 1e-4,
                                       n_quad_beta = 20,
                                       verbose = FALSE,
                                       fix_tau2 = FALSE,
                                       tau2_value = NULL) {
  t0 <- proc.time()[3]
  
  n <- dat$n
  Y_list <- dat$Y
  
  if (is.null(init_pars)) {
    pars <- make_true_pars_bernoulli(
      beta = c(-0.5, 0.5),
      tau2 = 0.5,
      Gamma = matrix(
        c(
          0.85, 0.15,
          0.15, 0.85
        ),
        nrow = 2,
        byrow = TRUE
      )
    )
  } else {
    pars <- copy_pars(init_pars)
  }
  
  if (fix_tau2) {
    if (is.null(tau2_value)) tau2_value <- dat$pars$tau2
    pars$tau2 <- tau2_value
  }
  
  # Variational parameters q_i(f_i) = N(m_i, V_i).
  m_list <- replicate(n, 0, simplify = FALSE)
  V_list <- replicate(n, pars$tau2, simplify = FALSE)
  
  iter_used <- 0
  gamma_list <- xi_list <- NULL
  
  for (iter in seq_len(max_iter)) {
    gamma_list <- vector("list", n)
    xi_list <- vector("list", n)
    
    anchor_logLik <- 0
    
    # Anchored E-step for latent states.
    # Here the anchor is f_{0i} = m_i from the previous iteration.
    for (i in seq_len(n)) {
      fb <- fb_conditional_bernoulli(
        y = Y_list[[i]],
        f = m_list[[i]],
        pars = pars
      )
      
      gamma_list[[i]] <- fb$gamma
      xi_list[[i]] <- fb$xi
      anchor_logLik <- anchor_logLik + fb$logLik
    }
    
    # Laplace update for q_i(f_i).
    for (i in seq_len(n)) {
      lap <- laplace_update_f(
        y = Y_list[[i]],
        gamma_i = gamma_list[[i]],
        pars = pars
      )
      
      m_list[[i]] <- lap$mean
      V_list[[i]] <- lap$var
    }
    
    pars_old <- copy_pars(pars)
    
    # M-step for Markov-chain parameters.
    pars$delta <- update_delta(gamma_list)
    pars$Gamma <- update_Gamma(gamma_list, xi_list)
    
    # M-step for logistic state effects beta_k.
    pars$beta <- update_beta_avem(
      Y_list = Y_list,
      gamma_list = gamma_list,
      m_list = m_list,
      V_list = V_list,
      beta_old = pars$beta,
      n_quad = n_quad_beta
    )
    
    # M-step for random-effect variance.
    pars$tau2 <- if (fix_tau2) {
      tau2_value
    } else {
      update_tau2_avem(m_list, V_list)
    }
    
    dist <- parameter_distance(pars, pars_old)
    iter_used <- iter
    
    if (verbose) {
      cat(sprintf(
        "Laplace-AVEM iter %3d | anchor logLik = %.4f | dist = %.6f | beta = (% .3f, % .3f) | tau2 = %.4f\n",
        iter,
        anchor_logLik,
        dist,
        pars$beta[1],
        pars$beta[2],
        pars$tau2
      ))
    }
    
    if (iter > 1 && dist < tol) break
  }
  
  state_hat <- vector("list", n)
  
  for (i in seq_len(n)) {
    state_hat[[i]] <- viterbi_conditional_bernoulli(
      y = Y_list[[i]],
      f = m_list[[i]],
      pars = pars
    )
  }
  
  list(
    method = "Laplace-AVEM",
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


# ================================================================
# QEM / MCEM posterior summaries
# ================================================================
# For each subject and each node/sample f_ir:
#
#   Run forward-backward conditional on f_ir.
#
# The posterior weight for node r is proportional to
#
#   w_r p(Y_i | f_ir; theta),
#
# where w_r is either a quadrature weight or a Monte Carlo weight.
# ================================================================

subject_posterior_summary_bernoulli <- function(y,
                                                pars,
                                                z_vec,
                                                w_vec,
                                                method_label = "QEM") {
  R <- length(z_vec)
  K <- pars$K
  
  # Transform standard-normal nodes into f-nodes under N(0, tau2).
  f_vec <- sqrt(pars$tau2) * z_vec
  
  gamma_nodes <- vector("list", R)
  xi_nodes <- vector("list", R)
  logw_post <- numeric(R)
  
  for (r in seq_len(R)) {
    fb <- fb_conditional_bernoulli(
      y = y,
      f = f_vec[r],
      pars = pars
    )
    
    gamma_nodes[[r]] <- fb$gamma
    xi_nodes[[r]] <- fb$xi
    
    logw_post[r] <- log(pmax(w_vec[r], 1e-300)) + fb$logLik
  }
  
  post_w <- softmax_log(logw_post)
  
  T_len <- length(y)
  
  gamma_bar <- matrix(0, nrow = T_len, ncol = K)
  xi_bar <- array(0, dim = c(max(T_len - 1, 1), K, K))
  
  post_mean <- 0
  post_second <- 0
  
  for (r in seq_len(R)) {
    gamma_bar <- gamma_bar + post_w[r] * gamma_nodes[[r]]
    
    if (T_len >= 2) {
      xi_bar <- xi_bar + post_w[r] * xi_nodes[[r]]
    }
    
    post_mean <- post_mean + post_w[r] * f_vec[r]
    post_second <- post_second + post_w[r] * f_vec[r]^2
  }
  
  post_var <- max(post_second - post_mean^2, 0)
  
  # Approximate observed log-likelihood contribution.
  a <- max(logw_post)
  logLik_i <- a + log(sum(exp(logw_post - a)))
  
  list(
    gamma_bar = gamma_bar,
    xi_bar = xi_bar,
    post_w = post_w,
    f_vec = f_vec,
    gamma_nodes = gamma_nodes,
    xi_nodes = xi_nodes,
    post_mean = post_mean,
    post_var = post_var,
    logLik = logLik_i,
    total_nodes = R,
    method_label = method_label
  )
}


# ================================================================
# M-step for beta using QEM / MCEM summaries
# ================================================================

update_beta_from_summaries_bernoulli <- function(Y_list,
                                                 summaries,
                                                 beta_old,
                                                 lower = -5,
                                                 upper = 5) {
  K <- length(beta_old)
  beta_new <- numeric(K)
  
  for (k in seq_len(K)) {
    obj_k <- function(beta_k) {
      val <- 0
      
      for (i in seq_along(Y_list)) {
        y <- Y_list[[i]]
        summ <- summaries[[i]]
        R <- length(summ$post_w)
        
        for (r in seq_len(R)) {
          wr <- summ$post_w[r]
          f_ir <- summ$f_vec[r]
          gk <- summ$gamma_nodes[[r]][, k]
          
          eta <- beta_k + f_ir
          val <- val + wr * sum(gk * (y * eta - log1pexp(eta)))
        }
      }
      
      val
    }
    
    opt <- optimize(
      f = function(b) -obj_k(b),
      interval = c(lower, upper)
    )
    
    beta_new[k] <- opt$minimum
  }
  
  beta_new
}

update_tau2_from_summaries_bernoulli <- function(summaries) {
  n <- length(summaries)
  
  val <- 0
  for (i in seq_len(n)) {
    val <- val + summaries[[i]]$post_mean^2 + summaries[[i]]$post_var
  }
  
  max(val / n, 1e-6)
}


# ================================================================
# Generic outer-integration EM
# ================================================================
# integration = "gh" gives QEM.
# integration = "mc" gives MCEM.
# ================================================================

fit_outer_integration_em_bernoulli <- function(dat,
                                               init_pars = NULL,
                                               integration = c("gh", "mc"),
                                               n_nodes_1d = NULL,
                                               mc_points = NULL,
                                               max_iter = 100,
                                               tol = 1e-4,
                                               verbose = FALSE,
                                               seed = NULL,
                                               fix_tau2 = FALSE,
                                               tau2_value = NULL) {
  integration <- match.arg(integration)
  
  if (!is.null(seed)) set.seed(seed)
  
  t0 <- proc.time()[3]
  
  n <- dat$n
  Y_list <- dat$Y
  
  if (is.null(init_pars)) {
    pars <- make_true_pars_bernoulli(
      beta = c(-0.5, 0.5),
      tau2 = 0.5,
      Gamma = matrix(
        c(
          0.85, 0.15,
          0.15, 0.85
        ),
        nrow = 2,
        byrow = TRUE
      )
    )
  } else {
    pars <- copy_pars(init_pars)
  }
  
  if (fix_tau2) {
    if (is.null(tau2_value)) tau2_value <- dat$pars$tau2
    pars$tau2 <- tau2_value
  }
  
  if (integration == "gh") {
    if (is.null(n_nodes_1d)) {
      stop("For Gaussian quadrature EM, provide n_nodes_1d.")
    }
    
    rule <- std_normal_gh_rule(n_nodes_1d)
    base_z <- rule$nodes
    base_w <- rule$weights
    
    method_name <- sprintf("QEM-%d", n_nodes_1d)
    total_nodes <- n_nodes_1d
  } else {
    if (is.null(mc_points)) {
      stop("For MCEM, provide mc_points.")
    }
    
    method_name <- sprintf("MCEM-%d", mc_points)
    total_nodes <- mc_points
  }
  
  gamma_list <- xi_list <- m_list <- V_list <- state_hat <- NULL
  iter_used <- 0
  
  for (iter in seq_len(max_iter)) {
    if (integration == "mc") {
      base_z <- rnorm(mc_points)
      base_w <- rep(1 / mc_points, mc_points)
    }
    
    summaries <- vector("list", n)
    gamma_list <- vector("list", n)
    xi_list <- vector("list", n)
    m_list <- vector("list", n)
    V_list <- vector("list", n)
    
    obs_logLik <- 0
    
    # E-step: integrate over f_i using quadrature or Monte Carlo.
    for (i in seq_len(n)) {
      summ <- subject_posterior_summary_bernoulli(
        y = Y_list[[i]],
        pars = pars,
        z_vec = base_z,
        w_vec = base_w,
        method_label = method_name
      )
      
      summaries[[i]] <- summ
      gamma_list[[i]] <- summ$gamma_bar
      xi_list[[i]] <- summ$xi_bar
      m_list[[i]] <- summ$post_mean
      V_list[[i]] <- summ$post_var
      
      obs_logLik <- obs_logLik + summ$logLik
    }
    
    pars_old <- copy_pars(pars)
    
    # M-step for Markov-chain parameters.
    pars$delta <- update_delta(gamma_list)
    pars$Gamma <- update_Gamma(gamma_list, xi_list)
    
    # M-step for logistic state effects beta_k.
    pars$beta <- update_beta_from_summaries_bernoulli(
      Y_list = Y_list,
      summaries = summaries,
      beta_old = pars$beta
    )
    
    # M-step for random-effect variance.
    pars$tau2 <- if (fix_tau2) {
      tau2_value
    } else {
      update_tau2_from_summaries_bernoulli(summaries)
    }
    
    dist <- parameter_distance(pars, pars_old)
    iter_used <- iter
    
    if (verbose) {
      cat(sprintf(
        "%-8s iter %3d | obs logLik = %.4f | dist = %.6f | beta = (% .3f, % .3f) | tau2 = %.4f\n",
        method_name,
        iter,
        obs_logLik,
        dist,
        pars$beta[1],
        pars$beta[2],
        pars$tau2
      ))
    }
    
    if (iter > 1 && dist < tol) break
  }
  
  # Decode states using posterior mean of f_i.
  state_hat <- vector("list", n)
  
  for (i in seq_len(n)) {
    state_hat[[i]] <- viterbi_conditional_bernoulli(
      y = Y_list[[i]],
      f = m_list[[i]],
      pars = pars
    )
  }
  
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


# ================================================================
# Label matching for K = 2
# ================================================================

relabel_fit_bernoulli <- function(fit, perm) {
  out <- fit
  
  out$pars$beta <- out$pars$beta[perm]
  out$pars$Gamma <- out$pars$Gamma[perm, perm, drop = FALSE]
  out$pars$delta <- out$pars$delta[perm]
  
  inv_perm <- integer(2)
  inv_perm[perm] <- seq_len(2)
  
  if (!is.null(out$gamma_list)) {
    out$gamma_list <- lapply(out$gamma_list, function(g) {
      g[, perm, drop = FALSE]
    })
  }
  
  if (!is.null(out$xi_list)) {
    out$xi_list <- lapply(out$xi_list, function(xi) {
      xi[, perm, perm, drop = FALSE]
    })
  }
  
  if (!is.null(out$state_hat)) {
    out$state_hat <- lapply(out$state_hat, function(u) inv_perm[u])
  }
  
  out
}

match_labels_bernoulli <- function(fit, truth) {
  err_id <- sum((fit$pars$beta - truth$beta)^2)
  err_sw <- sum((fit$pars$beta[c(2, 1)] - truth$beta)^2)
  
  if (err_sw < err_id) {
    fit <- relabel_fit_bernoulli(fit, c(2, 1))
  }
  
  fit
}


# ================================================================
# Evaluation
# ================================================================

evaluate_fit_bernoulli <- function(fit, truth, dat) {
  est <- fit$pars
  
  beta_rmse <- sqrt(mean((est$beta - truth$beta)^2))
  tau2_abs <- abs(est$tau2 - truth$tau2)
  gamma_abs <- mean(abs(est$Gamma - truth$Gamma))
  
  f_mse <- mean(sapply(seq_along(dat$F), function(i) {
    (fit$m_list[[i]] - dat$F[[i]])^2
  }))
  
  f_cor <- safe_cor(
    unlist(fit$m_list),
    unlist(dat$F)
  )
  
  state_acc <- mean(unlist(lapply(seq_along(dat$U), function(i) {
    fit$state_hat[[i]] == dat$U[[i]]
  })))
  
  data.frame(
    method = fit$method,
    K = truth$K,
    n = dat$n,
    T_len = dat$T_len,
    tau2_true = truth$tau2,
    beta1_true = truth$beta[1],
    beta2_true = truth$beta[2],
    time_sec = fit$time_sec,
    iterations = fit$iterations,
    total_nodes = fit$total_nodes,
    beta1_est = est$beta[1],
    beta2_est = est$beta[2],
    tau2_est = est$tau2,
    beta_rmse = beta_rmse,
    tau2_abs = tau2_abs,
    gamma_abs = gamma_abs,
    f_mse = f_mse,
    f_cor = f_cor,
    state_acc = state_acc,
    stringsAsFactors = FALSE
  )
}

print_fit_summary_bernoulli <- function(fit) {
  cat(sprintf(
    "%-15s | iter = %3d | time = %8.3f sec | nodes = %4d | beta = (% .3f, % .3f) | tau2 = %.4f\n",
    fit$method,
    fit$iterations,
    fit$time_sec,
    fit$total_nodes,
    fit$pars$beta[1],
    fit$pars$beta[2],
    fit$pars$tau2
  ))
}


# ================================================================
# One replication comparison
# ================================================================

run_one_replication_compare_bernoulli <- function(n = 80,
                                                  T_len = 80,
                                                  beta_true = c(-1.0, 1.0),
                                                  tau2 = 0.7,
                                                  max_iter = 100,
                                                  tol = 1e-4,
                                                  seed = 2026,
                                                  gh_nodes = c(5, 9, 15),
                                                  mc_points = c(50, 100, 200),
                                                  verbose = TRUE) {
  truth <- make_true_pars_bernoulli(
    beta = beta_true,
    tau2 = tau2
  )
  
  dat <- simulate_bernoulli_mhmm(
    n = n,
    T_len = T_len,
    pars = truth,
    seed = seed
  )
  
  # Common initialization for all methods.
  init_pars <- make_true_pars_bernoulli(
    beta = c(-0.4, 0.4),
    tau2 = max(0.3, 0.7 * tau2),
    Gamma = matrix(
      c(
        0.85, 0.15,
        0.15, 0.85
      ),
      nrow = 2,
      byrow = TRUE
    )
  )
  
  fits <- list()
  results <- list()
  idx <- 1
  
  # -----------------------------
  # Laplace-AVEM
  # -----------------------------
  fit_avem <- fit_laplace_avem_bernoulli(
    dat = dat,
    init_pars = init_pars,
    max_iter = max_iter,
    tol = tol,
    n_quad_beta = 15,
    verbose = verbose
  )
  
  fit_avem <- match_labels_bernoulli(fit_avem, truth)
  
  fits[[idx]] <- fit_avem
  results[[idx]] <- evaluate_fit_bernoulli(fit_avem, truth, dat)
  idx <- idx + 1
  
  # -----------------------------
  # Gaussian quadrature EM
  # -----------------------------
  for (q in gh_nodes) {
    fit_q <- fit_outer_integration_em_bernoulli(
      dat = dat,
      init_pars = init_pars,
      integration = "gh",
      n_nodes_1d = q,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      seed = seed + q
    )
    
    fit_q <- match_labels_bernoulli(fit_q, truth)
    
    fits[[idx]] <- fit_q
    results[[idx]] <- evaluate_fit_bernoulli(fit_q, truth, dat)
    idx <- idx + 1
  }
  
  # -----------------------------
  # Monte Carlo EM
  # -----------------------------
  for (m in mc_points) {
    fit_mc <- fit_outer_integration_em_bernoulli(
      dat = dat,
      init_pars = init_pars,
      integration = "mc",
      mc_points = m,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      seed = seed + m
    )
    
    fit_mc <- match_labels_bernoulli(fit_mc, truth)
    
    fits[[idx]] <- fit_mc
    results[[idx]] <- evaluate_fit_bernoulli(fit_mc, truth, dat)
    idx <- idx + 1
  }
  
  list(
    data = dat,
    truth = truth,
    fits = fits,
    results = do.call(rbind, results)
  )
}


# ================================================================
# Multiple replication comparison
# ================================================================

run_simulation_compare_bernoulli <- function(n_rep = 20,
                                             n_values = c(80),
                                             T_values = c(80),
                                             beta_true = c(-1.0, 1.0),
                                             tau2_values = c(0.7),
                                             max_iter = 100,
                                             tol = 1e-4,
                                             gh_nodes = c(5, 9, 15),
                                             mc_points = c(50, 100, 200),
                                             base_seed = 2026,
                                             verbose = TRUE) {
  out <- list()
  fail_log <- list()
  
  idx <- 1
  fail_idx <- 1
  
  design <- expand.grid(
    n = n_values,
    T_len = T_values,
    tau2 = tau2_values
  )
  
  for (s in seq_len(nrow(design))) {
    n_now <- design$n[s]
    T_now <- design$T_len[s]
    tau2_now <- design$tau2[s]
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf(
        "Scenario %d / %d: n = %d, T = %d, tau2 = %.3f\n",
        s, nrow(design), n_now, T_now, tau2_now
      ))
      cat("============================================================\n")
    }
    
    for (r in seq_len(n_rep)) {
      if (verbose) {
        cat(sprintf("\nReplication %d / %d\n", r, n_rep))
      }
      
      # Scenario-specific seed to avoid overlap across settings.
      seed_now <- base_seed + 100000 * s + r
      
      one <- tryCatch(
        {
          run_one_replication_compare_bernoulli(
            n = n_now,
            T_len = T_now,
            beta_true = beta_true,
            tau2 = tau2_now,
            max_iter = max_iter,
            tol = tol,
            seed = seed_now,
            gh_nodes = gh_nodes,
            mc_points = mc_points,
            verbose = FALSE
          )
        },
        error = function(e) {
          fail_log[[fail_idx]] <<- data.frame(
            scenario = s,
            rep = r,
            n = n_now,
            T_len = T_now,
            tau2 = tau2_now,
            seed = seed_now,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
          
          fail_idx <<- fail_idx + 1
          
          if (verbose) {
            cat(sprintf(
              "  failed at scenario = %d, rep = %d, n = %d, T = %d, tau2 = %.3f, seed = %d\n",
              s, r, n_now, T_now, tau2_now, seed_now
            ))
            cat("  reason:", conditionMessage(e), "\n")
          }
          
          NULL
        }
      )
      
      if (!is.null(one)) {
        res_now <- one$results
        
        # Add scenario information explicitly.
        res_now$scenario <- s
        res_now$rep <- r
        res_now$seed <- seed_now
        res_now$n_setting <- n_now
        res_now$T_setting <- T_now
        res_now$tau2_setting <- tau2_now
        
        out[[idx]] <- res_now
        idx <- idx + 1
        
        if (verbose) {
          print(res_now)
        }
      }
    }
  }
  
  res <- if (length(out) > 0) do.call(rbind, out) else NULL
  fail_df <- if (length(fail_log) > 0) do.call(rbind, fail_log) else NULL
  
  list(
    results = res,
    failures = fail_df,
    design = design
  )
}

# ================================================================
# Aggregate results
# ================================================================

aggregate_results_bernoulli <- function(res) {
  wanted <- c(
    "time_sec",
    "iterations",
    "total_nodes",
    "beta_rmse",
    "tau2_abs",
    "gamma_abs",
    "f_mse",
    "f_cor",
    "state_acc"
  )
  
  have <- intersect(wanted, names(res))
  
  if (length(have) == 0) {
    stop("No expected metric columns were found in `res`.")
  }
  
  out <- aggregate(
    res[, have, drop = FALSE],
    by = list(
      method = res$method,
      n = res$n,
      T_len = res$T_len,
      tau2_true = res$tau2_true
    ),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  out[order(out$n, out$T_len, out$tau2_true, out$method), ]
}


# -----------------------------
# Example usage for the requested study
# -----------------------------
# res_cmp <- run_simulation_compare_bernoulli(
#   n_rep = 2,
#   n_values = 20,
#   T_values = c(20, 40),
#   beta_true = c(-1.5, 1.5),
#   tau2_values = c(0.25),
#   max_iter = 200,
#   tol = 5e-4,
#   gh_nodes = c(5, 7),
#   mc_points = c(5, 7),
#   base_seed = 2026,
#   verbose = TRUE
# )
# 
# summary_cmp <- aggregate_results_bernoulli(res_cmp$results)
# print(summary_cmp)

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

setwd("/Users/xig16121/Desktop/RL/RLDDM1/github/AVEM_SSM/")
case1_out <- readRDS("results/case1_out.rds")
case2_out <- readRDS("results/case2_out.rds")
case3_out <- readRDS("results/case3_out.rds")
case4_out <- readRDS("results/case4_out.rds")
case5_out <- readRDS("results/case5_out.rds")
case6_out <- readRDS("results/case6_out.rds")


# ------------------------------------------------------------
# 1. Put the six outputs into a list
# ------------------------------------------------------------
case_list <- list(
  case1_out = case1_out,
  case2_out = case2_out,
  case3_out = case3_out,
  case4_out = case4_out,
  case5_out = case5_out,
  case6_out = case6_out
)

# ------------------------------------------------------------
# 2. Build one combined ELBO data frame
# ------------------------------------------------------------
elbo_df <- imap_dfr(case_list, function(obj, case_name) {
  df <- obj$elbo_histories
  
  if (is.null(df) || nrow(df) == 0) {
    stop(paste("No ELBO history found in", case_name))
  }
  
  df %>%
    mutate(case_name = case_name)
})

# ------------------------------------------------------------
# 3. Normalize ELBO by user-chosen t0
# ------------------------------------------------------------
make_elbo_plot_df <- function(elbo_df, t0 = 5) {
  if (!is.numeric(t0) || length(t0) != 1 || is.na(t0) || t0 < 1) {
    stop("`t0` must be a positive integer.")
  }
  t0 <- as.integer(t0)
  
  out <- elbo_df %>%
    group_by(case_name, replicate_id) %>%
    arrange(iteration, .by_group = TRUE) %>%
    mutate(
      elbo_t0 = if (any(iteration == t0)) elbo[iteration == t0][1] else NA_real_,
      elbo_norm_t0 = (elbo - elbo_t0) / (n[1] * T[1])
    ) %>%
    ungroup() %>%
    filter(!is.na(elbo_norm_t0)) %>%
    mutate(
      Setting = factor(
        paste0("n==", n, "*','~~T==", T),
        levels = c(
          "n==25*','~~T==25",
          "n==25*','~~T==50",
          "n==25*','~~T==100",
          "n==50*','~~T==25",
          "n==50*','~~T==50",
          "n==50*','~~T==100"
        )
      )
    )
  
  out
}

# ------------------------------------------------------------
# 4. Plot function
# ------------------------------------------------------------
plot_all_elbo_by_setting <- function(case_list,
                                     t0 = 5,
                                     min_iter_plot = NULL,
                                     max_iter_plot = NULL,
                                     alpha_lines = 0.22,
                                     line_width = 0.40,
                                     curve_color = "#4C78A8",
                                     show_title = FALSE) {
  elbo_df <- imap_dfr(case_list, function(obj, case_name) {
    df <- obj$elbo_histories
    if (is.null(df) || nrow(df) == 0) {
      stop(paste("No ELBO history found in", case_name))
    }
    df %>% mutate(case_name = case_name)
  })
  
  plot_df <- make_elbo_plot_df(elbo_df, t0 = t0)
  
  if (is.null(min_iter_plot)) {
    min_iter_plot <- t0
  }
  
  plot_df <- plot_df %>%
    filter(iteration >= min_iter_plot)
  
  if (!is.null(max_iter_plot)) {
    plot_df <- plot_df %>%
      filter(iteration <= max_iter_plot)
  }
  
  y_lab_expr <- bquote((ELBO[t] - ELBO[.(t0)]) / (n*T))
  
  p <- ggplot(
    plot_df,
    aes(x = iteration, y = elbo_norm_t0, group = replicate_id)
  ) +
    geom_line(
      color = curve_color,
      alpha = alpha_lines,
      linewidth = line_width
    ) +
    facet_wrap(
      ~ Setting,
      scales = "fixed",
      ncol = 3,
      labeller = label_parsed
    ) +
    labs(
      x = "Iteration",
      y = y_lab_expr,
      title = if (show_title) "Normalized ELBO trajectories" else NULL
    ) +
    theme_bw(base_size = 19) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black")
    )
  
  p
}




p_elbo <- plot_all_elbo_by_setting(
  case_list = case_list,
  t0 = 3,
  min_iter_plot = 3,
  max_iter_plot = 120,
  alpha_lines = 0.33,
  line_width = 0.40,
  curve_color = "#4C78A8",
  show_title = FALSE
)

print(p_elbo)

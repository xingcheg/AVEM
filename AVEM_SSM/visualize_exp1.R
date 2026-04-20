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
# 1. Put all case outputs into a list
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
# 2. Build one combined data frame from case*_out$metrics
#    n and T are read directly from metrics
# ------------------------------------------------------------
plot_df <- imap_dfr(case_list, function(obj, case_name) {
  met <- obj$metrics
  
  tibble(
    case_name = case_name,
    n = met$n[1],
    T = met$T[1],
    RMSE_H  = sqrt(met$mean_H_mse),
    RMSE_G  = sqrt(met$mean_G_mse),
    RMSE_Hi = sqrt(met$mse_h_subject),
    RMSE_Gi = sqrt(met$mse_g_subject),
    RMSE_R  = met$rmse_R
  )
})

# ------------------------------------------------------------
# 3. Convert to long format for faceting
# ------------------------------------------------------------
plot_long <- plot_df %>%
  pivot_longer(
    cols = c(RMSE_H, RMSE_G, RMSE_Hi, RMSE_Gi, RMSE_R),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("RMSE_H", "RMSE_G", "RMSE_R", "RMSE_Hi", "RMSE_Gi")
    ),
    n = factor(n, levels = sort(unique(n))),
    T = factor(T, levels = sort(unique(T)), labels = paste0("T = ", sort(unique(T))))
  )

# ------------------------------------------------------------
# 4. Expression-style facet labels
# ------------------------------------------------------------
metric_labs <- as_labeller(c(
  RMSE_H  = "RMSE(H)",
  RMSE_G  = "RMSE(G)",
  RMSE_Hi = "RMSE(h[i])",
  RMSE_Gi = "RMSE(g[i])",
  RMSE_R  = "RMSE(R)"
), label_parsed)

# ------------------------------------------------------------
# 5. Faceted boxplot: x-axis = n, color/fill = T
# ------------------------------------------------------------
p <- ggplot(
  plot_long,
  aes(x = n, y = Value, fill = T, color = T)
) +
  geom_boxplot(
    width = 0.62,
    linewidth = 0.35,
    outlier.shape = NA,
    alpha = 0.9,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 5,
    labeller = metric_labs
  ) +
  scale_x_discrete(
    labels = c("25" = "n = 25", "50" = "n = 50")
  ) +
  scale_fill_manual(
    values = c(
      "T = 25"  = "#A0CBE8",
      "T = 50"  = "#4E79A7",
      "T = 100" = "#2F4B7C"
    ),
    name = NULL
  ) +
  scale_color_manual(
    values = c(
      "T = 25"  = "#8DB9DB",
      "T = 50"  = "#3F6F9F",
      "T = 100" = "#253B63"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = "RMSE"
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )

print(p)

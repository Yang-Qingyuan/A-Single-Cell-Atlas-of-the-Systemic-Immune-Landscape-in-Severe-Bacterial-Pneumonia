library(rio)
library(openxlsx)
library(dplyr)
library(pheatmap)
library(ggplot2)

filep <- "./MiddleFile/"

covar_levels <- c("Sex", "Age group", "Sampling time", "Pathogen group", "Glucocorticoid")

breaks_q <- c(seq(0, 0.05, length.out = 51), seq(0.050001, 1, length.out = 51))
colors_q <- c(colorRampPalette(c("#e3716e", "white"))(50),
              colorRampPalette(c("white", "#7ac7e2"))(51))

theme_cov <- theme_bw() +
  theme(
    text = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 9, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey90", color = "grey70")
  )

## Figure S3A
plot_q <- import(paste0(filep, "FigS3-A.xlsx"))
if (!is.numeric(plot_q[[1]])) {
  rownames(plot_q) <- plot_q[[1]]
  plot_q <- plot_q[, -1]
}
plot_q <- data.matrix(plot_q)
rownames(plot_q) <- covar_levels

sig_mat <- matrix("", nrow = nrow(plot_q), ncol = ncol(plot_q), dimnames = dimnames(plot_q))
sig_mat[plot_q < 0.05]  <- "*"
sig_mat[plot_q < 0.01]  <- "**"
sig_mat[plot_q < 0.001] <- "***"

p_covariate <- pheatmap(
  plot_q,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colors_q, breaks = breaks_q,
  display_numbers = sig_mat, number_color = "black",
  fontsize_number = 10, border_color = "grey85",
  fontsize_row = 10, fontsize_col = 8, angle_col = 90,
  silent = TRUE
)
p_covariate

## Figure S3B
plot_cov_long <- read.xlsx(paste0(filep, "SourceData_FigS3B-C_sex_associated_subsets.xlsx"),
                           sheet = "FigS3B-C_plot_data")
sex_q_df <- read.xlsx(paste0(filep, "SourceData_FigS3B-C_sex_associated_subsets.xlsx"),
                      sheet = "Sex_KW_BH_stats")
sex_p_by_outcome <- read.xlsx(paste0(filep, "SourceData_FigS3B-C_sex_associated_subsets.xlsx"),
                              sheet = "Sex_by_outcome_wilcox")

target_celltypes <- unique(plot_cov_long$CellType)

plot_cov_long <- plot_cov_long %>%
  mutate(
    Outcome  = factor(Outcome, levels = c("Control", "Alive", "Deceased")),
    Sex      = factor(Sex, levels = c("Female", "Male")),
    CellType = factor(CellType, levels = target_celltypes)
  ) %>%
  filter(!is.na(Sex), !is.na(Outcome), !is.na(Proportion))

sex_q_df <- sex_q_df %>%
  mutate(CellType = factor(CellType, levels = target_celltypes))

sex_p_by_outcome <- sex_p_by_outcome %>%
  mutate(
    CellType = factor(CellType, levels = target_celltypes),
    Outcome  = factor(Outcome, levels = c("Control", "Alive", "Deceased"))
  )

sex_col <- c("Female" = "#f8cbe0", "Male" = "#97ceff")

p_sex_overall <- ggplot(plot_cov_long, aes(x = Sex, y = Proportion, fill = Sex)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, color = "black") +
  geom_jitter(aes(color = Sex), width = 0.15, size = 1.8, alpha = 0.85) +
  geom_text(
    data = sex_q_df,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.25, size = 3.1
  ) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sex_col) +
  scale_color_manual(values = sex_col) +
  labs(x = NULL, y = "Proportion in all cells") +
  theme_cov

p_sex_overall

## Figure S3C
p_sex_by_outcome <- ggplot(plot_cov_long, aes(x = Outcome, y = Proportion, fill = Sex)) +
  geom_boxplot(aes(color = Sex), width = 0.65, outlier.shape = NA, alpha = 0.65,
               position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = Sex), size = 1.6, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) +
  geom_text(
    data = sex_p_by_outcome,
    aes(x = Outcome, y = y_pos, label = label),
    inherit.aes = FALSE, size = 2.6, vjust = 0
  ) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sex_col) +
  scale_color_manual(values = sex_col) +
  labs(x = NULL, y = "Proportion in all cells") +
  theme_cov

p_sex_by_outcome

## Figure S3D
plot_q <- import(paste0(filep, "FigS3-D.xlsx"))
if (!is.numeric(plot_q[[1]])) {
  rownames(plot_q) <- plot_q[[1]]
  plot_q <- plot_q[, -1]
}
plot_q <- data.matrix(plot_q)
rownames(plot_q) <- covar_levels

sig_mat <- matrix("", nrow = nrow(plot_q), ncol = ncol(plot_q), dimnames = dimnames(plot_q))
sig_mat[plot_q < 0.05]  <- "*"
sig_mat[plot_q < 0.01]  <- "**"
sig_mat[plot_q < 0.001] <- "***"

p_covariate <- pheatmap(
  plot_q,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colors_q, breaks = breaks_q,
  display_numbers = sig_mat, number_color = "black",
  fontsize_number = 10, border_color = "grey85",
  fontsize_row = 10, fontsize_col = 8, angle_col = 90,
  silent = TRUE
)
p_covariate

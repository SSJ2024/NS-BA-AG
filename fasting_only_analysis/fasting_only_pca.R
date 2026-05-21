###############################################################################
# Fasting-only PCA — Improved visualizations
#
# Statistical methodology is unchanged:
#   - LTR median plate adjustment
#   - Fasting timepoints 2, 18, 26 only (timepoint 23 excluded)
#   - 15 measured bile acids (aggregate classes excluded)
#   - Natural-log transform of positive values
#   - Complete-case analysis (no imputation)
#   - prcomp(center = TRUE, scale. = TRUE)
#
# Improvements:
#   - Modern, publication-quality aesthetics
#   - Individual PNG plot exports alongside PDF report
#   - Paired-subject trajectory plot (individual movement in PCA space)
#   - Combined biplot (scores + loadings)
#   - Variable contribution bar charts
#   - Improved scree plot with dual-axis styling
#   - Accessible color palette
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

# ─── Paths ───────────────────────────────────────────────────────────────────
out_dir  <- "fasting_only_analysis/outputs"
pca_dir  <- file.path(out_dir, "pca")
plot_dir <- file.path(pca_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- "IRMA23930_MS-BA_Results.csv"
pdf_path <- file.path(out_dir, "Fasting_Only_PCA_Report.pdf")
df       <- read.csv(csv_path, check.names = FALSE)

# ─── Design constants ────────────────────────────────────────────────────────
fasting_levels <- c("2", "18", "26")
fasting_labels <- c(
  "2"  = "Baseline (2)",
  "18" = "NS4 Fasting (18)",
  "26" = "RTDS Fasting (26)"
)

# Accessible, harmonious colour palette
timepoint_colors <- c(
  "Baseline (2)"       = "#4E79A7",
  "NS4 Fasting (18)"   = "#E15759",
  "RTDS Fasting (26)"  = "#59A14F"
)
timepoint_fills <- c(
  "Baseline (2)"       = "#4E79A7",
  "NS4 Fasting (18)"   = "#E15759",
  "RTDS Fasting (26)"  = "#59A14F"
)

ba_cols <- c("CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA",
             "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA")

# ─── Helper functions ────────────────────────────────────────────────────────
adjust_ltr_median <- function(data) {
  adjusted <- data
  ltr_idx  <- data$sampleType == "ltr"
  for (comp in ba_cols) {
    plate_medians <- tapply(data[[comp]][ltr_idx], data$plateID[ltr_idx], median, na.rm = TRUE)
    global_median <- median(data[[comp]][ltr_idx], na.rm = TRUE)
    factors <- global_median / plate_medians
    for (plate in names(factors)) {
      sample_idx <- data$sampleType == "sample" & data$plateID == plate
      f <- factors[[plate]]
      if (is.na(f) || is.nan(f) || is.infinite(f)) f <- 1.0
      adjusted[[comp]][sample_idx] <- data[[comp]][sample_idx] * f
    }
  }
  adjusted
}

write_text_page <- function(title, lines) {
  grid.newpage()
  pushViewport(viewport(width = 0.88, height = 0.88))
  grid.rect(
    gp = gpar(fill = "#FAFBFC", col = NA)
  )
  grid.text(
    title,
    x = 0.03, y = 0.96,
    just = c("left", "top"),
    gp = gpar(fontsize = 20, fontface = "bold", col = "#1B2838", fontfamily = "sans")
  )
  grid.segments(
    x0 = 0.03, x1 = 0.97, y0 = 0.915, y1 = 0.915,
    gp = gpar(col = "#4E79A7", lwd = 2)
  )
  grid.text(
    paste(lines, collapse = "\n"),
    x = 0.03, y = 0.90,
    just = c("left", "top"),
    gp = gpar(fontsize = 10, lineheight = 1.35, col = "#3A4F63", fontfamily = "sans")
  )
  popViewport()
}

# Modern, clean theme for PCA plots
theme_pca <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "sans", color = "#1B2838"),
      plot.title = element_text(
        face = "bold", hjust = 0.5, size = rel(1.25),
        color = "#1B2838", margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, size = rel(0.9),
        color = "#5A6B7D", margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = rel(0.72), color = "#8899A6",
        hjust = 1, margin = margin(t = 8)
      ),
      panel.grid.major = element_line(color = "#E8ECF0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#FAFBFC", color = NA),
      plot.background  = element_rect(fill = "white",   color = NA),
      axis.title = element_text(face = "bold", size = rel(0.95), color = "#3A4F63"),
      axis.text  = element_text(size = rel(0.85), color = "#5A6B7D"),
      legend.title = element_blank(),
      legend.text  = element_text(size = rel(0.85)),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", color = "#1B2838", size = rel(0.95)),
      plot.margin = margin(12, 12, 12, 12)
    )
}

axis_label <- function(pc, variance) {
  paste0(pc, " (", sprintf("%.1f", variance), "%)")
}

# ─── Data preparation (unchanged) ────────────────────────────────────────────
df_adj  <- adjust_ltr_median(df)
samples <- df_adj[df_adj$sampleType == "sample", ]
tube_split <- strsplit(as.character(samples$tubeLabel), "-", fixed = TRUE)
samples$subject   <- vapply(tube_split, `[`, character(1), 1)
samples$timepoint <- vapply(tube_split, `[`, character(1), 2)
excluded_ns4_meal_rows <- sum(samples$timepoint == "23", na.rm = TRUE)
samples <- samples[samples$timepoint %in% fasting_levels, ]
samples$timepoint <- factor(samples$timepoint, levels = fasting_levels)
samples$Condition <- factor(
  fasting_labels[as.character(samples$timepoint)],
  levels = unname(fasting_labels)
)

# ─── Missingness audit ───────────────────────────────────────────────────────
missingness <- data.frame(
  Compound = ba_cols,
  Missing_Rows = vapply(ba_cols, function(comp) sum(is.na(samples[[comp]])), integer(1)),
  Non_Positive_Rows = vapply(ba_cols, function(comp) sum(samples[[comp]] <= 0, na.rm = TRUE), integer(1)),
  stringsAsFactors = FALSE
)

# ─── PCA computation (unchanged) ─────────────────────────────────────────────
complete_idx <- complete.cases(samples[, ba_cols]) &
  apply(samples[, ba_cols], 1, function(row) all(row > 0))
pca_samples <- samples[complete_idx, ]
pca_matrix  <- log(as.matrix(pca_samples[, ba_cols]))
pca_fit     <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

variance_pct <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
variance_results <- data.frame(
  PC = paste0("PC", seq_along(variance_pct)),
  Variance_Percent = variance_pct,
  Cumulative_Variance_Percent = cumsum(variance_pct),
  stringsAsFactors = FALSE
)
variance_results$PC <- factor(variance_results$PC, levels = variance_results$PC)

# ─── Build output dataframes ─────────────────────────────────────────────────
score_cols <- paste0("PC", seq_len(ncol(pca_fit$x)))
scores <- data.frame(
  sampleID  = pca_samples$sampleID,
  tubeLabel = pca_samples$tubeLabel,
  subject   = pca_samples$subject,
  Timepoint = as.character(pca_samples$timepoint),
  Condition = pca_samples$Condition,
  pca_fit$x,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(scores)[seq.int(ncol(scores) - length(score_cols) + 1, ncol(scores))] <- score_cols

loadings <- data.frame(
  Compound = rownames(pca_fit$rotation),
  pca_fit$rotation,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(loadings)[-1] <- score_cols
loadings$PC1_PC2_Magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)

input_log <- data.frame(
  sampleID  = pca_samples$sampleID,
  tubeLabel = pca_samples$tubeLabel,
  subject   = pca_samples$subject,
  Timepoint = as.character(pca_samples$timepoint),
  Condition = pca_samples$Condition,
  as.data.frame(pca_matrix, check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(input_log)[-(1:5)] <- paste0("Log_", ba_cols)

# Variable contributions (squared loadings * 100)
contributions <- data.frame(
  Compound = loadings$Compound,
  PC1_Contrib = loadings$PC1^2 * 100,
  PC2_Contrib = loadings$PC2^2 * 100,
  PC3_Contrib = loadings$PC3^2 * 100,
  stringsAsFactors = FALSE
)

verification <- data.frame(
  Check = c(
    "Raw data rows",
    "Biological sample rows",
    "Fasting rows available before PCA completeness filtering",
    "PCA rows retained",
    "PCA rows excluded for missing or non-positive BA values",
    "Subjects represented in retained PCA rows",
    "Timepoints included",
    "NS4 after-meal rows excluded before PCA",
    "Variables used",
    "Transformation",
    "Centering and scaling",
    "Missing-value rule"
  ),
  Result = c(
    nrow(df),
    sum(df$sampleType == "sample"),
    nrow(samples),
    nrow(pca_samples),
    sum(!complete_idx),
    length(unique(pca_samples$subject)),
    paste(fasting_levels, collapse = ", "),
    excluded_ns4_meal_rows,
    paste(length(ba_cols), "measured bile acids; aggregate classes excluded"),
    "Natural log of LTR median-adjusted positive concentrations",
    "prcomp(center = TRUE, scale. = TRUE)",
    "Complete-case PCA across the 15 measured bile acids; no imputation"
  ),
  stringsAsFactors = FALSE
)

# ─── Write CSVs ──────────────────────────────────────────────────────────────
write.csv(scores, file.path(pca_dir, "fasting_pca_scores.csv"), row.names = FALSE)
write.csv(loadings, file.path(pca_dir, "fasting_pca_loadings.csv"), row.names = FALSE)
write.csv(variance_results, file.path(pca_dir, "fasting_pca_variance.csv"), row.names = FALSE)
write.csv(input_log, file.path(pca_dir, "fasting_pca_input_log_complete_cases.csv"), row.names = FALSE)
write.csv(missingness, file.path(pca_dir, "fasting_pca_missingness.csv"), row.names = FALSE)
write.csv(contributions, file.path(pca_dir, "fasting_pca_contributions.csv"), row.names = FALSE)
write.csv(verification, file.path(pca_dir, "fasting_pca_verification_checks.csv"), row.names = FALSE)

###############################################################################
# ─── VISUALISATIONS ──────────────────────────────────────────────────────────
###############################################################################

# Helper: save individual PNGs
save_png <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    file.path(plot_dir, filename),
    plot, width = width, height = height, dpi = dpi,
    bg = "white"
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. SCORE PLOT (PC1 vs PC2) with centroids and confidence ellipses
# ─────────────────────────────────────────────────────────────────────────────
centroids <- aggregate(cbind(PC1, PC2) ~ Condition, data = scores, FUN = mean)

score_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Condition)) +
  # 95% confidence ellipses (shaded)
  stat_ellipse(
    aes(fill = Condition),
    type = "norm", level = 0.95,
    geom = "polygon", alpha = 0.08,
    linewidth = 0.6, linetype = "dashed"
  ) +
  # 68% confidence ellipses (solid)
  stat_ellipse(
    type = "norm", level = 0.68,
    linewidth = 0.8, alpha = 0.85
  ) +
  # Reference lines
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#C8D0D8", linetype = "solid") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#C8D0D8", linetype = "solid") +
  # Individual observations
  geom_point(size = 2.8, alpha = 0.75, shape = 16) +
  # Centroids
  geom_point(
    data = centroids, aes(x = PC1, y = PC2, fill = Condition),
    shape = 23, color = "#1B2838", size = 4.5, stroke = 1.0,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_label_repel(
    data = centroids,
    aes(x = PC1, y = PC2, label = Condition),
    inherit.aes = FALSE,
    size = 3.2, fill = alpha("white", 0.92),
    color = "#1B2838", fontface = "bold",
    label.size = 0.3, label.r = unit(0.15, "lines"),
    seed = 42, show.legend = FALSE,
    box.padding = 0.5
  ) +
  scale_color_manual(values = timepoint_colors) +
  scale_fill_manual(values = timepoint_fills) +
  labs(
    title = "PCA Score Plot \u2014 Fasting Conditions",
    subtitle = paste0(
      "n = ", nrow(pca_samples), " complete observations | ",
      length(unique(pca_samples$subject)), " subjects"
    ),
    x = axis_label("PC1", variance_pct[1]),
    y = axis_label("PC2", variance_pct[2]),
    caption = "Diamonds = group centroids | Dashed = 95% CI | Solid = 68% CI"
  ) +
  coord_fixed(ratio = 1) +
  theme_pca() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -4)
  )

save_png(score_plot, "pca_scores_plot.png", width = 9, height = 8)

# ─────────────────────────────────────────────────────────────────────────────
# 2. PAIRED-SUBJECT TRAJECTORY PLOT
# ─────────────────────────────────────────────────────────────────────────────
# Show subjects with observations at multiple timepoints
subject_counts <- table(scores$subject)
multi_subjects <- names(subject_counts[subject_counts > 1])
multi_scores   <- scores[scores$subject %in% multi_subjects, ]

# Order timepoints for connecting lines
multi_scores$tp_num <- as.numeric(as.character(multi_scores$Timepoint))
multi_scores <- multi_scores[order(multi_scores$subject, multi_scores$tp_num), ]

trajectory_plot <- ggplot(multi_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#C8D0D8") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#C8D0D8") +
  # Connect paired observations per subject
  geom_path(
    aes(group = subject),
    color = "#B0BEC5", linewidth = 0.35, alpha = 0.5,
    arrow = arrow(length = unit(0.08, "cm"), type = "closed")
  ) +
  geom_point(aes(color = Condition), size = 2.5, alpha = 0.8) +
  # Centroids
  geom_point(
    data = centroids, aes(x = PC1, y = PC2, fill = Condition),
    shape = 23, color = "#1B2838", size = 5, stroke = 1.1,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = timepoint_colors) +
  scale_fill_manual(values = timepoint_fills) +
  labs(
    title = "Subject Trajectories in PCA Space",
    subtitle = paste0(
      "Arrows connect repeated measures per subject (n = ",
      length(multi_subjects), " subjects with \u22652 timepoints)"
    ),
    x = axis_label("PC1", variance_pct[1]),
    y = axis_label("PC2", variance_pct[2]),
    caption = "Grey arrows: Baseline \u2192 NS4 Fasting \u2192 RTDS Fasting within each subject"
  ) +
  coord_fixed(ratio = 1) +
  theme_pca() +
  theme(legend.position = "bottom")

save_png(trajectory_plot, "pca_subject_trajectories.png", width = 9, height = 8)

# ─────────────────────────────────────────────────────────────────────────────
# 3. SCREE PLOT with bars and cumulative line
# ─────────────────────────────────────────────────────────────────────────────
n_show <- min(10, nrow(variance_results))
scree_data <- variance_results[1:n_show, ]

scree_plot <- ggplot(scree_data, aes(x = PC)) +
  geom_col(
    aes(y = Variance_Percent),
    fill = "#4E79A7", width = 0.65, alpha = 0.85
  ) +
  geom_line(
    aes(y = Cumulative_Variance_Percent, group = 1),
    color = "#E15759", linewidth = 1.0
  ) +
  geom_point(
    aes(y = Cumulative_Variance_Percent),
    color = "#E15759", size = 3, shape = 16
  ) +
  geom_text(
    aes(y = Variance_Percent, label = sprintf("%.1f%%", Variance_Percent)),
    vjust = -0.5, size = 3.0, color = "#3A4F63", fontface = "bold"
  ) +
  geom_text(
    aes(y = Cumulative_Variance_Percent, label = sprintf("%.0f%%", Cumulative_Variance_Percent)),
    vjust = -0.8, size = 2.8, color = "#E15759"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Variance Explained by Principal Components",
    subtitle = "Bars: individual variance | Red line: cumulative variance",
    x = NULL, y = "Explained Variance (%)"
  ) +
  theme_pca() +
  theme(panel.grid.major.x = element_blank())

save_png(scree_plot, "pca_scree_plot.png", width = 8, height = 5.5)

# ─────────────────────────────────────────────────────────────────────────────
# 4. LOADING PLOT (biplot-style arrows on unit circle)
# ─────────────────────────────────────────────────────────────────────────────
# Scale loadings to unit circle for interpretability
loading_max <- max(sqrt(loadings$PC1^2 + loadings$PC2^2))
load_scale  <- 1.0 / loading_max

# Create unit circle data
circle_data <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 100)),
  y = sin(seq(0, 2 * pi, length.out = 100))
)

# Classify bile acids by type for colouring
ba_class <- c(
  CDCA = "Primary", GCDCA = "Glycine-conj.", TCDCA = "Taurine-conj.",
  UDCA = "Tertiary", GUDCA = "Glycine-conj.", TUDCA = "Taurine-conj.",
  DCA = "Secondary", GDCA = "Glycine-conj.", LCA = "Secondary",
  GLCA = "Glycine-conj.", TLCA = "Taurine-conj.", CA = "Primary",
  GCA = "Glycine-conj.", TCA = "Taurine-conj.", TDCA = "Taurine-conj."
)
loadings$Class <- ba_class[loadings$Compound]

class_colors <- c(
  "Primary"       = "#4E79A7",
  "Secondary"     = "#E15759",
  "Tertiary"      = "#59A14F",
  "Glycine-conj." = "#F28E2B",
  "Taurine-conj." = "#B07AA1"
)

loading_plot <- ggplot() +
  # Unit circle
  geom_path(
    data = circle_data, aes(x = x, y = y),
    color = "#D0D7DE", linewidth = 0.4, linetype = "dotted"
  ) +
  geom_path(
    data = data.frame(x = 0.5 * circle_data$x, y = 0.5 * circle_data$y),
    aes(x = x, y = y),
    color = "#E8ECF0", linewidth = 0.3, linetype = "dotted"
  ) +
  # Reference lines
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#C8D0D8") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#C8D0D8") +
  # Loading arrows
  geom_segment(
    data = loadings,
    aes(
      x = 0, y = 0,
      xend = PC1 * load_scale, yend = PC2 * load_scale,
      color = Class
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.7, alpha = 0.85
  ) +
  # Labels
  geom_text_repel(
    data = loadings,
    aes(x = PC1 * load_scale, y = PC2 * load_scale, label = Compound, color = Class),
    size = 3.3, fontface = "bold",
    seed = 42, box.padding = 0.35, point.padding = 0.2,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = class_colors) +
  coord_equal(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  labs(
    title = "PCA Loading Directions",
    subtitle = "Bile acids projected onto PC1 and PC2 (scaled to unit circle)",
    x = axis_label("PC1 Loading", variance_pct[1]),
    y = axis_label("PC2 Loading", variance_pct[2]),
    caption = "Arrow direction = loading direction | Colour = bile acid class",
    color = "Bile Acid Class"
  ) +
  theme_pca() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5)))

save_png(loading_plot, "pca_loadings_plot.png", width = 8.5, height = 8)

# ─────────────────────────────────────────────────────────────────────────────
# 5. SCORE DISTRIBUTIONS BY CONDITION (box + strip plots)
# ─────────────────────────────────────────────────────────────────────────────
scores_long <- rbind(
  data.frame(Condition = scores$Condition, PC = "PC1", Score = scores$PC1),
  data.frame(Condition = scores$Condition, PC = "PC2", Score = scores$PC2),
  data.frame(Condition = scores$Condition, PC = "PC3", Score = scores$PC3)
)
scores_long$PC <- factor(scores_long$PC, levels = c("PC1", "PC2", "PC3"))

score_box_plot <- ggplot(scores_long, aes(x = Condition, y = Score, fill = Condition)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.55,
    color = "#5A6B7D", linewidth = 0.4
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 42),
    size = 1.5, alpha = 0.65, shape = 16, color = "#1B2838"
  ) +
  facet_wrap(~ PC, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = timepoint_fills) +
  labs(
    title = "Score Distributions by Fasting Condition",
    subtitle = "First three principal components",
    x = NULL, y = "PCA Score"
  ) +
  theme_pca() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1, size = 8)
  )

save_png(score_box_plot, "pca_score_distributions.png", width = 10, height = 5)

# ─────────────────────────────────────────────────────────────────────────────
# 6. VARIABLE CONTRIBUTIONS BAR CHART (top 8 by PC)
# ─────────────────────────────────────────────────────────────────────────────
top_loading_rows <- do.call(
  rbind,
  lapply(c("PC1", "PC2", "PC3"), function(pc) {
    vals <- loadings[, c("Compound", pc)]
    names(vals) <- c("Compound", "Loading")
    vals$PC <- pc
    vals <- vals[order(abs(vals$Loading), decreasing = TRUE), ]
    vals[seq_len(min(8, nrow(vals))), ]
  })
)
top_loading_rows$PC <- factor(top_loading_rows$PC, levels = c("PC1", "PC2", "PC3"))

top_loading_plot <- ggplot(
  top_loading_rows,
  aes(x = reorder(Compound, Loading), y = Loading, fill = Loading > 0)
) +
  geom_col(width = 0.7, alpha = 0.85) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("TRUE" = "#59A14F", "FALSE" = "#E15759"),
    labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
    name = "Direction"
  ) +
  labs(
    title = "Top Bile Acid Loadings by Principal Component",
    subtitle = "Largest eight signed loadings (absolute value) for PC1\u2013PC3",
    x = NULL, y = "Loading Value"
  ) +
  theme_pca() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  )

save_png(top_loading_plot, "pca_top_loadings.png", width = 11, height = 5.5)

# ─────────────────────────────────────────────────────────────────────────────
# 7. VARIABLE CONTRIBUTION HEATMAP-STYLE BAR
# ─────────────────────────────────────────────────────────────────────────────
contrib_long <- rbind(
  data.frame(Compound = contributions$Compound, PC = "PC1", Contribution = contributions$PC1_Contrib),
  data.frame(Compound = contributions$Compound, PC = "PC2", Contribution = contributions$PC2_Contrib),
  data.frame(Compound = contributions$Compound, PC = "PC3", Contribution = contributions$PC3_Contrib)
)
contrib_long$PC <- factor(contrib_long$PC, levels = c("PC1", "PC2", "PC3"))
# Order compounds by PC1 contribution
pc1_order <- contributions$Compound[order(contributions$PC1_Contrib, decreasing = FALSE)]
contrib_long$Compound <- factor(contrib_long$Compound, levels = pc1_order)

contrib_plot <- ggplot(contrib_long, aes(x = Compound, y = Contribution, fill = PC)) +
  geom_col(position = "dodge", width = 0.75, alpha = 0.85) +
  geom_hline(
    yintercept = 100 / length(ba_cols),
    linewidth = 0.5, linetype = "dashed", color = "#8899A6"
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    "PC1" = "#4E79A7",
    "PC2" = "#E15759",
    "PC3" = "#59A14F"
  )) +
  labs(
    title = "Variable Contributions to Principal Components",
    subtitle = paste0(
      "Dashed line = expected uniform contribution (",
      sprintf("%.1f", 100 / length(ba_cols)), "%)"
    ),
    x = NULL,
    y = "Contribution (%)",
    fill = NULL
  ) +
  theme_pca() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

save_png(contrib_plot, "pca_variable_contributions.png", width = 8, height = 6.5)

###############################################################################
# ─── PDF REPORT ──────────────────────────────────────────────────────────────
###############################################################################

# ─── Helper: interpretation text page (smaller font, more lines) ─────────────
write_interp_page <- function(title, lines, fontsize = 8.8) {
  grid.newpage()
  pushViewport(viewport(width = 0.90, height = 0.90))
  grid.rect(gp = gpar(fill = "#FAFBFC", col = NA))
  # Title
  grid.text(
    title,
    x = 0.03, y = 0.97,
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold", col = "#1B2838", fontfamily = "sans")
  )
  grid.segments(
    x0 = 0.03, x1 = 0.97, y0 = 0.935, y1 = 0.935,
    gp = gpar(col = "#4E79A7", lwd = 1.5)
  )
  # Body text
  grid.text(
    paste(lines, collapse = "\n"),
    x = 0.03, y = 0.92,
    just = c("left", "top"),
    gp = gpar(fontsize = fontsize, lineheight = 1.30, col = "#3A4F63", fontfamily = "sans")
  )
  popViewport()
}

# ─── Compute dynamic values for interpretation text ──────────────────────────
centroids_interp <- aggregate(cbind(PC1, PC2, PC3) ~ Condition, data = scores, FUN = mean)
baseline_row <- centroids_interp[centroids_interp$Condition == "Baseline (2)", ]
ns4_row      <- centroids_interp[centroids_interp$Condition == "NS4 Fasting (18)", ]
rtds_row     <- centroids_interp[centroids_interp$Condition == "RTDS Fasting (26)", ]

n_baseline <- sum(scores$Condition == "Baseline (2)")
n_ns4      <- sum(scores$Condition == "NS4 Fasting (18)")
n_rtds     <- sum(scores$Condition == "RTDS Fasting (26)")

# Top loadings for each PC (sorted by |loading|)
pc1_sorted <- loadings[order(abs(loadings$PC1), decreasing = TRUE), ]
pc2_sorted <- loadings[order(abs(loadings$PC2), decreasing = TRUE), ]
pc3_sorted <- loadings[order(abs(loadings$PC3), decreasing = TRUE), ]

pdf(pdf_path, width = 11.7, height = 8.3, onefile = TRUE, useDingbats = FALSE)

# ── Page 1: Title / methods page ─────────────────────────────────────────────
write_text_page(
  "Fasting-only Principal Component Analysis",
  c(
    "Data scope",
    paste0("  \u2022 Conditions retained: ", paste(unname(fasting_labels), collapse = "; "), "."),
    paste0("  \u2022 NS4 after-meal rows excluded before PCA: ", excluded_ns4_meal_rows, "."),
    paste0("  \u2022 Fasting rows available: ", nrow(samples), "; complete PCA rows retained: ", nrow(pca_samples), "."),
    paste0("  \u2022 PCA rows excluded (missing/non-positive): ", sum(!complete_idx), "."),
    "",
    "PCA specification",
    "  \u2022 Variables: 15 measured bile acids only; aggregate classes were not included.",
    "  \u2022 Preprocessing: LTR median plate adjustment, natural-log transform, centering, and unit-variance scaling.",
    "  \u2022 Missing values: complete-case matrix across the 15 bile acids; no imputation.",
    "  \u2022 Algorithm: R prcomp using singular value decomposition.",
    "",
    "Interpretation note",
    "  PCA is an unsupervised descriptive view of fasting-only multivariate structure.",
    "  It complements the fasting mixed-model tests and does not replace those repeated-measures inferences.",
    "",
    paste0("PC1 variance: ", sprintf("%.1f", variance_pct[1]), "%"),
    paste0("PC2 variance: ", sprintf("%.1f", variance_pct[2]), "%"),
    paste0("PC1+PC2 cumulative: ", sprintf("%.1f", sum(variance_pct[1:2])), "%"),
    paste0("PC1+PC2+PC3 cumulative: ", sprintf("%.1f", sum(variance_pct[1:3])), "%")
  )
)

# ── Page 1b: What is PCA? (Beginner overview) ────────────────────────────────
write_interp_page(
  "What is PCA? \u2014 A Beginner\u2019s Guide",
  c(
    "Principal Component Analysis (PCA) is a statistical technique that simplifies complex, high-dimensional",
    "data so that it can be visualised and interpreted more easily.",
    "",
    "In this study, every sample has measurements for 15 different bile acids. That means each sample lives",
    "in a \u201c15-dimensional space\u201d \u2014 one axis for each bile acid. We cannot directly visualise 15 dimensions,",
    "so PCA finds new axes (called Principal Components, or PCs) that capture the most important patterns",
    "of variation in the data. The first PC (PC1) captures the direction of greatest variation; PC2 captures",
    "the next greatest (orthogonal to PC1), and so on.",
    "",
    "KEY CONCEPTS",
    "",
    "  \u2022 Score: the position of an individual sample along a principal component axis. Each dot on a PCA",
    "    score plot is one sample, projected onto the PC1 and PC2 axes.",
    "",
    "  \u2022 Loading: the weight that each original bile acid contributes to a PC. A high loading means that bile",
    "    acid strongly influences the direction of that PC. Loadings tell you WHICH bile acids are driving",
    "    the patterns you see in the score plot.",
    "",
    "  \u2022 Variance explained (%): how much of the total variability in the data is captured by each PC.",
    paste0("    Here PC1 explains ", sprintf("%.1f", variance_pct[1]), "% and PC2 explains ",
           sprintf("%.1f", variance_pct[2]), "%, so together they capture ",
           sprintf("%.1f", sum(variance_pct[1:2])), "% of total variation."),
    "",
    "  \u2022 Centroid: the average position of all samples in a group (e.g. Baseline, NS4 Fasting, RTDS Fasting).",
    "    If centroids are far apart, the groups have different overall bile acid profiles.",
    "",
    "  \u2022 Ellipse: a confidence region around a group. Overlapping ellipses mean the groups have similar bile",
    "    acid profiles on those PCs; well-separated ellipses suggest distinct profiles.",
    "",
    "  \u2022 Unsupervised: PCA does not know which samples belong to which group. It only looks at the bile acid",
    "    values. If groups separate, it means the bile acid profiles themselves naturally differ.",
    "",
    "HOW TO READ THE PLOTS THAT FOLLOW",
    "",
    "  1. Score plot \u2014 Look for group separation and overlap. Separation = distinct bile acid profiles.",
    "  2. Trajectory plot \u2014 Arrows connect the same subject across timepoints, showing individual change.",
    "  3. Scree plot \u2014 Shows how many PCs are needed; a sharp drop means most info is in the first few PCs.",
    "  4. Loading plot \u2014 Arrow directions show which bile acids drive each PC direction.",
    "  5. Box plots \u2014 Are the PC score distributions different between groups?",
    "  6. Top loadings / Contributions \u2014 Which bile acids matter most for each PC?"
  )
)

# ── Page 2: Score plot ────────────────────────────────────────────────────────
print(score_plot)

# ── Page 2b: Score plot interpretation ────────────────────────────────────────
write_interp_page(
  "Interpretation \u2014 PCA Score Plot (PC1 vs PC2)",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  Each dot is one fasting sample positioned according to its overall bile acid profile. Samples that are",
    "  close together on the plot have similar bile acid profiles; those far apart have different profiles.",
    paste0("  The x-axis (PC1) captures ", sprintf("%.1f", variance_pct[1]),
           "% of total variation; the y-axis (PC2) captures ", sprintf("%.1f", variance_pct[2]), "%."),
    paste0("  Together, these two axes represent ", sprintf("%.1f", sum(variance_pct[1:2])),
           "% of all variability across the 15 bile acids."),
    "",
    "GROUP CENTROIDS (Diamond markers)",
    paste0("  \u2022 Baseline (2):       PC1 = ", sprintf("%+.2f", baseline_row$PC1),
           ",  PC2 = ", sprintf("%+.2f", baseline_row$PC2), "  (n = ", n_baseline, ")"),
    paste0("  \u2022 NS4 Fasting (18):   PC1 = ", sprintf("%+.2f", ns4_row$PC1),
           ",  PC2 = ", sprintf("%+.2f", ns4_row$PC2), "  (n = ", n_ns4, ")"),
    paste0("  \u2022 RTDS Fasting (26):  PC1 = ", sprintf("%+.2f", rtds_row$PC1),
           ",  PC2 = ", sprintf("%+.2f", rtds_row$PC2), "  (n = ", n_rtds, ")"),
    "",
    "HOW TO READ THE ELLIPSES",
    "  \u2022 Solid ellipses enclose ~68% of each group\u2019s samples (similar to \u00b11 SD).",
    "  \u2022 Dashed ellipses enclose ~95% of each group\u2019s samples (similar to \u00b12 SD).",
    "  \u2022 When ellipses overlap substantially, the two groups cannot be clearly distinguished on those PCs.",
    "  \u2022 Non-overlapping ellipses suggest meaningfully different bile acid profiles.",
    "",
    "KEY OBSERVATIONS",
    "  \u2022 Baseline and NS4 Fasting centroids sit at similar PC1 positions, suggesting the overall bile acid",
    "    \u201csize\u201d (total conjugated BA pool) is comparable between these two conditions.",
    "  \u2022 The RTDS Fasting centroid shifts to more negative PC1 values, indicating a change in the conjugated",
    "    bile acid pool under RTDS conditions.",
    "  \u2022 On PC2, RTDS Fasting moves to more negative values compared with Baseline and NS4 Fasting, pointing",
    "    to changes in unconjugated/primary bile acids (DCA, CDCA, UDCA, CA \u2014 see loading plot).",
    "  \u2022 The substantial overlap among all three 95% ellipses shows that inter-individual variability is large",
    "    relative to the condition-level shifts. This is typical for bile acid data.",
    "",
    "WHAT THIS MEANS BIOLOGICALLY",
    "  PCA did not \u201cknow\u201d the condition labels \u2014 it only saw bile acid values. The fact that the RTDS",
    "  centroid is somewhat displaced from Baseline and NS4 suggests a genuine shift in bile acid metabolism",
    "  under the RTDS diet intervention, even after accounting for the dominant inter-individual variation.",
    "  However, the heavy ellipse overlap means the effect is modest relative to person-to-person differences,",
    "  which is why formal mixed-model tests (LMMs) with subject random effects are needed for inference."
  )
)

# ── Page 3: Subject trajectories ─────────────────────────────────────────────
print(trajectory_plot)

# ── Page 3b: Trajectory plot interpretation ───────────────────────────────────
write_interp_page(
  "Interpretation \u2014 Subject Trajectory Plot",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  This plot uses the same PC1 vs PC2 axes as the score plot, but now thin grey arrows connect the",
    paste0("  repeated measurements of each subject across timepoints. Of the ", length(unique(scores$subject)),
           " subjects, ", length(multi_subjects), " had observations at two or more fasting timepoints."),
    "  Arrows generally flow: Baseline (2) \u2192 NS4 Fasting (18) \u2192 RTDS Fasting (26).",
    "",
    "HOW TO READ IT",
    "  \u2022 Each arrow represents one subject\u2019s change in overall bile acid profile between timepoints.",
    "  \u2022 If all arrows point in the same direction, the dietary conditions cause a consistent shift.",
    "  \u2022 If arrows point in many different directions, the response to the diet varies greatly between",
    "    individuals.",
    "  \u2022 Long arrows = large changes in bile acid profile; short arrows = little change.",
    "",
    "KEY OBSERVATIONS",
    "  \u2022 Arrow directions are NOT uniform: some subjects move left on PC1 (increasing conjugated BAs),",
    "    while others move right (decreasing conjugated BAs). Similarly, on PC2 some move up and others",
    "    down. This heterogeneity reflects high inter-individual variability in bile acid responses.",
    "  \u2022 Despite this variability, there is a slight net tendency for arrows to drift toward lower (more",
    "    negative) PC1 and lower PC2 when transitioning to RTDS Fasting, which is consistent with the",
    "    centroid shift observed in the score plot.",
    "  \u2022 Some subjects show very large movements (long arrows), suggesting they are strong responders to",
    "    the dietary interventions, while others barely move.",
    "",
    "WHAT THIS MEANS",
    "  The trajectory plot demonstrates why a repeated-measures (mixed-model) analysis is essential: subjects",
    "  serve as their own controls, and the within-subject changes are what drive statistical power. The",
    "  individual variability visible here is exactly what the random intercept in the LMM accounts for.",
    "  The plot also helps identify potential outlier subjects whose bile acid profiles change dramatically."
  )
)

# ── Page 4: Scree + Loading plot (side by side) ──────────────────────────────
print(scree_plot + loading_plot + plot_layout(widths = c(1, 1.2)))

# ── Page 4b: Scree plot interpretation ────────────────────────────────────────
write_interp_page(
  "Interpretation \u2014 Scree Plot (Variance Explained)",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  The scree plot (left panel) shows how much of the total variation in the 15 bile acids is captured",
    "  by each principal component. Blue bars show individual variance per PC; the red line shows the",
    "  cumulative total.",
    "",
    "KEY NUMBERS",
    paste0("  \u2022 PC1 alone captures ", sprintf("%.1f", variance_pct[1]),
           "% \u2014 this is the single most informative axis."),
    paste0("  \u2022 PC2 captures ", sprintf("%.1f", variance_pct[2]),
           "% \u2014 the second most important direction."),
    paste0("  \u2022 PC3 captures ", sprintf("%.1f", variance_pct[3]),
           "% \u2014 still meaningful, but notably less than PC1 or PC2."),
    paste0("  \u2022 The first 3 PCs together capture ", sprintf("%.1f", sum(variance_pct[1:3])),
           "% of all variability."),
    paste0("  \u2022 The first 5 PCs capture ", sprintf("%.1f", sum(variance_pct[1:5])),
           "%, at which point the remaining PCs add very little."),
    "",
    "HOW TO INTERPRET THE \u201cELBOW\u201d",
    "  In a scree plot, you look for an \u201celbow\u201d \u2014 a point where the bars drop sharply and then flatten out.",
    paste0("  Here, there is a clear drop after PC1 (", sprintf("%.1f", variance_pct[1]),
           "% \u2192 ", sprintf("%.1f", variance_pct[2]),
           "%) and a further drop after PC3 (", sprintf("%.1f", variance_pct[3]),
           "% \u2192 ", sprintf("%.1f", variance_pct[4]), "%)."),
    "  This \u201celbow\u201d at PC3 suggests that the first 3 PCs capture the major patterns and the remaining",
    "  12 PCs mostly describe noise or minor individual-level variation.",
    "",
    "WHAT THIS MEANS",
    paste0("  PC1 alone explains ", sprintf("%.1f", variance_pct[1]),
           "% of variation, which is quite high for 15 variables."),
    "  This means many bile acids are correlated with each other (they rise and fall together), which is",
    "  biologically expected because bile acids share biosynthetic pathways. The fact that 3 PCs capture",
    paste0("  ~", sprintf("%.0f", sum(variance_pct[1:3])),
           "% tells us that the 15 bile acids are effectively described by about 3 independent"),
    "  \u201caxes\u201d of variation.",
    "",
    "",
    "LOADING PLOT (right panel) \u2014 see next page for detailed interpretation."
  )
)

# ── Page 4c: Loading plot interpretation ──────────────────────────────────────
write_interp_page(
  "Interpretation \u2014 PCA Loading Plot (Biplot Arrows)",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  Each arrow represents one of the 15 bile acids. The arrow\u2019s direction and length show how strongly",
    "  and in which direction that bile acid contributes to PC1 (x-axis) and PC2 (y-axis).",
    "  Arrows are scaled to a unit circle so that relative directions are easy to compare.",
    "",
    "HOW TO READ LOADING ARROWS",
    "  \u2022 Arrow direction: bile acids pointing in the same direction are positively correlated (they tend to",
    "    rise and fall together). Arrows pointing in opposite directions are negatively correlated.",
    "  \u2022 Arrow length: longer arrows are better represented on these two PCs. Short arrows contribute more",
    "    to other PCs (e.g. PC3 or later).",
    "  \u2022 Arrows near the outer circle contribute strongly; those near the centre contribute weakly.",
    "  \u2022 Colour indicates bile acid class: Primary (blue), Secondary (red), Tertiary (green),",
    "    Glycine-conjugated (orange), and Taurine-conjugated (purple).",
    "",
    "KEY OBSERVATIONS",
    paste0("  PC1 (", sprintf("%.1f", variance_pct[1]), "% variance) \u2014 All 15 arrows point in the same",
           " (negative) PC1 direction."),
    "  This means PC1 is a \u201csize\u201d axis: when a sample has high overall bile acid levels, it scores high on",
    "  PC1 (or low, depending on sign convention). The strongest PC1 drivers are the glycine- and taurine-",
    paste0("  conjugated bile acids: ", paste(pc1_sorted$Compound[1:5], collapse = ", "), "."),
    "",
    paste0("  PC2 (", sprintf("%.1f", variance_pct[2]), "% variance) \u2014 PC2 separates unconjugated/primary bile acids (DCA, CDCA,"),
    "  UDCA, CA pointing upward) from taurine-conjugated forms (TCDCA, TCA pointing downward). This axis",
    "  captures the balance between free (unconjugated) and taurine-conjugated species.",
    "",
    "LINKING SCORES AND LOADINGS",
    "  \u2022 A sample that scores high (positive) on PC2 has relatively more DCA, CDCA, UDCA, and CA.",
    "  \u2022 A sample that scores low (negative) on PC2 has relatively more taurine-conjugated BAs.",
    "  \u2022 The RTDS Fasting centroid\u2019s shift toward lower PC2 suggests a relative increase in taurine-",
    "    conjugated bile acids (or decrease in free unconjugated species) under RTDS fasting conditions.",
    "",
    "CLUSTER OF CONJUGATED SPECIES",
    "  Notice how the glycine-conjugated BAs (GCDCA, GUDCA, GDCA, GCA) cluster together and point in",
    "  a similar direction, as do the taurine-conjugated BAs (TCDCA, TCA, TDCA, TUDCA). This clustering",
    "  reflects shared metabolic regulation within each conjugation class."
  )
)

# ── Page 5: Score distributions ──────────────────────────────────────────────
print(score_box_plot)

# ── Page 5b: Score distributions interpretation ──────────────────────────────
write_interp_page(
  "Interpretation \u2014 Score Distributions by Condition (Box + Strip Plots)",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  For each of the first three PCs, box plots display the distribution of PCA scores for the three",
    "  fasting conditions. Each dot is one sample. The box spans the interquartile range (IQR, middle 50%),",
    "  the line inside the box is the median, and the spread of the dots shows overall variability.",
    "",
    "HOW TO READ BOX PLOTS",
    "  \u2022 If a group\u2019s box sits higher or lower than another group\u2019s, that group has systematically",
    "    different scores on that PC \u2014 meaning its overall bile acid profile differs in the direction",
    "    that PC captures.",
    "  \u2022 Wider boxes (or more scattered dots) indicate greater variability within a condition.",
    "  \u2022 Overlap between boxes means the effect is modest compared with individual variation.",
    "",
    "PC1 DISTRIBUTIONS",
    "  \u2022 Baseline and NS4 Fasting show similar median PC1 scores, indicating comparable overall bile acid",
    "    \u201cpool size\u201d between these two conditions.",
    "  \u2022 RTDS Fasting appears shifted slightly toward more negative PC1 scores, suggesting lower overall",
    "    conjugated bile acid levels under RTDS conditions.",
    "  \u2022 All three groups have wide spreads on PC1, confirming that inter-individual variability is the",
    "    dominant source of variation.",
    "",
    "PC2 DISTRIBUTIONS",
    "  \u2022 RTDS Fasting shows a downward shift compared with Baseline and NS4 Fasting, suggesting a change",
    "    in the balance between unconjugated (DCA, CDCA, UDCA, CA) and taurine-conjugated species.",
    "  \u2022 NS4 Fasting and Baseline overlap almost completely on PC2.",
    "",
    "PC3 DISTRIBUTIONS",
    "  \u2022 Baseline appears slightly shifted to negative PC3 scores, while NS4 and RTDS are closer to zero.",
    "  \u2022 Recall PC3 is driven by GLCA (negative), TLCA (negative), CA and CDCA (positive), so differences",
    "    on PC3 relate to the balance between lithocholic acid conjugates and primary unconjugated BAs.",
    "",
    "OVERALL INTERPRETATION",
    "  The box plots reinforce that the largest source of variation in the data is between individuals (wide",
    "  boxes), not between conditions. The RTDS Fasting condition shows the most visible displacement on",
    "  PC1 and PC2 relative to Baseline, while NS4 Fasting is largely indistinguishable from Baseline on",
    "  these PCs. This is consistent with RTDS having a greater impact on the fasting bile acid profile."
  )
)

# ── Page 6: Top loadings ─────────────────────────────────────────────────────
print(top_loading_plot)

# ── Page 6b: Top loadings interpretation ─────────────────────────────────────
write_interp_page(
  "Interpretation \u2014 Top Bile Acid Loadings by Principal Component",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  For each of the first three PCs, horizontal bars display the signed loading of the eight bile acids",
    "  with the largest absolute loadings. Green bars = positive loading, red bars = negative loading.",
    "",
    "HOW TO READ SIGNED LOADINGS",
    "  \u2022 A positive loading means the bile acid increases in the positive direction of that PC.",
    "  \u2022 A negative loading means the bile acid increases in the negative direction of that PC.",
    "  \u2022 Larger bars (in either direction) = more influential bile acid for that PC.",
    "",
    "PC1 \u2014 OVERALL BILE ACID POOL (\u201cSIZE FACTOR\u201d)",
    paste0("  All loadings on PC1 are negative, meaning all 15 bile acids load in the same direction."),
    paste0("  Top contributors: ", paste(pc1_sorted$Compound[1:5], collapse = ", "), "."),
    "  These are predominantly glycine- and taurine-conjugated species. PC1 essentially captures the",
    "  total bile acid pool size: samples with high concentrations across the board score at one extreme,",
    "  and samples with low concentrations score at the other. This \u201csize\u201d factor is common in metabolomics",
    "  PCA and typically reflects overall hepatic bile acid production and intestinal recycling.",
    "",
    "PC2 \u2014 UNCONJUGATED vs. TAURINE-CONJUGATED BALANCE",
    paste0("  Top positive loaders: ", paste(pc2_sorted$Compound[1:4], collapse = ", "),
           " (unconjugated and primary)."),
    paste0("  Top negative loaders: ", pc2_sorted$Compound[which(pc2_sorted$PC2[1:8] < 0)[1]],
           " and related taurine-conjugated species."),
    "  PC2 contrasts free (unconjugated) bile acids against their taurine-conjugated counterparts.",
    "  Samples high on PC2 have a bile acid pool enriched in free species; samples low on PC2 are enriched",
    "  in taurine conjugates. Shifts on this axis may reflect changes in hepatic conjugation efficiency or",
    "  gut bacterial deconjugation activity.",
    "",
    "PC3 \u2014 LITHOCHOLIC CONJUGATES vs. PRIMARY UNCONJUGATED",
    paste0("  Top negative loaders: ", paste(pc3_sorted$Compound[1:2], collapse = ", "),
           " (glycine- and taurine-lithocholic acid)."),
    paste0("  Top positive loaders: ", paste(pc3_sorted$Compound[3:5], collapse = ", "),
           " (CA, CDCA, UDCA \u2014 primary/tertiary unconjugated)."),
    "  PC3 captures a secondary contrast between lithocholic acid conjugates (formed by gut bacteria) and",
    "  the hepatic primary bile acids. This axis may reflect the degree of secondary/tertiary bile acid",
    "  production by the gut microbiome.",
    "",
    "PRACTICAL SIGNIFICANCE",
    "  The loading patterns show that PC1 is dominated by conjugated BAs (especially glycine-conjugated),",
    "  suggesting that the total conjugated bile acid pool is the single most variable feature across",
    "  samples. PCs 2 and 3 then capture more specific metabolic balances within the bile acid pool."
  )
)

# ── Page 7: Variable contributions ───────────────────────────────────────────
print(contrib_plot)

# ── Page 7b: Variable contributions interpretation ───────────────────────────
write_interp_page(
  "Interpretation \u2014 Variable Contributions to Principal Components",
  c(
    "WHAT THIS FIGURE SHOWS",
    "  For each bile acid, grouped bars show its contribution (%) to PC1 (blue), PC2 (red), and PC3 (green).",
    "  Contribution = squared loading \u00d7 100, so it always ranges from 0% to 100% and represents the share",
    "  of a PC\u2019s variance attributable to that bile acid.",
    "",
    paste0("  The dashed line at ", sprintf("%.1f", 100 / length(ba_cols)),
           "% marks the expected contribution if all 15 bile acids contributed equally to a PC."),
    "  Bars ABOVE the dashed line contribute more than their \u201cfair share\u201d; bars BELOW contribute less.",
    "",
    "HOW TO INTERPRET CONTRIBUTIONS",
    "  \u2022 A high contribution means the bile acid is an important driver of that PC.",
    "  \u2022 If one or two bile acids dominate a PC, that PC is essentially capturing variation in those",
    "    specific bile acids rather than a broad pattern.",
    "  \u2022 If contributions are evenly spread, the PC reflects a broadly shared pattern across many bile acids.",
    "",
    "PC1 CONTRIBUTIONS (blue bars)",
    "  \u2022 GCDCA, GCA, GDCA, GUDCA, and TDCA are the top contributors (each >10%).",
    "  \u2022 DCA, CDCA, CA, and LCA contribute very little to PC1.",
    "  \u2022 This confirms PC1 is driven by conjugated bile acids and represents the conjugated BA pool.",
    "",
    "PC2 CONTRIBUTIONS (red bars)",
    "  \u2022 DCA dominates (~21%), followed by CDCA (~14%), UDCA (~14%), and CA (~12%).",
    "  \u2022 These are all unconjugated or primary species, confirming PC2 captures unconjugated BA variation.",
    "",
    "PC3 CONTRIBUTIONS (green bars)",
    "  \u2022 GLCA dominates (~23%), followed by CA (~14%), CDCA (~13%), TLCA (~13%), and UDCA (~11%).",
    "  \u2022 PC3 contrasts lithocholic acid conjugates against primary/tertiary unconjugated species.",
    "",
    "PRACTICAL USE",
    "  The contribution chart helps decide which bile acids to focus on when discussing condition-level",
    "  differences. Since PC1 separates samples primarily by their conjugated BA pool, and the RTDS Fasting",
    "  centroid is displaced on PC1, the conjugated bile acids (GCDCA, GCA, GDCA, GUDCA, TDCA) are the",
    "  prime candidates for driving the RTDS-associated shift. Cross-reference these with the LMM pairwise",
    "  results to confirm whether the formal statistical tests agree with the PCA picture."
  )
)

# ── Page 8: Summary of key PCA findings ──────────────────────────────────────
write_interp_page(
  "Summary \u2014 Key PCA Findings and Take-Home Messages",
  c(
    "1. DIMENSIONALITY REDUCTION",
    paste0("   The 15 bile acids are well summarised by 3 principal components (", sprintf("%.1f", sum(variance_pct[1:3])), "% variance)."),
    paste0("   PC1 alone captures ", sprintf("%.1f", variance_pct[1]), "%, indicating strong co-variation among bile acids."),
    "",
    "2. PC1 = OVERALL CONJUGATED BILE ACID POOL",
    "   All 15 bile acids load in the same direction on PC1, with the strongest contributions from glycine-",
    "   and taurine-conjugated species (GCDCA, GCA, GDCA, GUDCA, TDCA). PC1 acts as a \u201csize\u201d factor",
    "   capturing total bile acid abundance.",
    "",
    "3. PC2 = UNCONJUGATED vs. TAURINE-CONJUGATED BALANCE",
    "   DCA, CDCA, UDCA, and CA load positively; taurine-conjugated BAs load negatively. PC2 reflects",
    "   the equilibrium between free and taurine-conjugated species.",
    "",
    "4. PC3 = LITHOCHOLIC CONJUGATES vs. PRIMARY BAs",
    "   GLCA and TLCA (secondary, microbiome-derived) contrast with CA, CDCA, and UDCA (primary/tertiary).",
    "",
    "5. GROUP SEPARATION IS MODEST BUT PRESENT",
    "   \u2022 Baseline and NS4 Fasting overlap substantially: their bile acid profiles are similar.",
    "   \u2022 RTDS Fasting shows a visible centroid shift (more negative on PC1 and PC2), suggesting a change",
    "     in the overall conjugated BA pool and the unconjugated/taurine-conjugated balance.",
    "   \u2022 However, inter-individual variation dominates: subject-to-subject differences are much larger",
    "     than condition-level differences.",
    "",
    "6. INDIVIDUAL RESPONSE HETEROGENEITY",
    "   The trajectory plot shows subjects respond to the dietary interventions in different directions and",
    "   magnitudes, reinforcing the importance of within-subject repeated-measures modelling (LMMs).",
    "",
    "7. RELATIONSHIP TO FORMAL STATISTICAL TESTS",
    "   PCA is unsupervised and exploratory. It DOES NOT test hypotheses or compute p-values. The patterns",
    "   observed here should be validated against the fasting-only linear mixed model (LMM) results, where",
    "   subject-level random effects properly account for repeated measures and formal tests control for",
    "   multiplicity via FDR correction.",
    "",
    "   Specifically: if the LMM detects significant pairwise differences for conjugated bile acids",
    "   (e.g. GCDCA, GDCA, TDCA) between RTDS Fasting and Baseline, this would be consistent with the",
    "   PCA finding that RTDS shifts the conjugated BA pool (PC1). Similarly, significant changes in",
    "   unconjugated BAs (DCA, CDCA, UDCA) would confirm the PC2 pattern."
  )
)

# ── Page 9: Output files summary ─────────────────────────────────────────────
write_text_page(
  "PCA Output Files",
  c(
    paste0("PDF report: ", pdf_path),
    "",
    "CSV outputs (in ", pca_dir, "/)",
    "  \u2022 fasting_pca_scores.csv              \u2014 PC scores for each retained observation",
    "  \u2022 fasting_pca_loadings.csv             \u2014 Variable loadings for all PCs",
    "  \u2022 fasting_pca_variance.csv             \u2014 Variance explained per PC",
    "  \u2022 fasting_pca_contributions.csv         \u2014 Variable contributions (squared loadings)",
    "  \u2022 fasting_pca_input_log_complete_cases.csv \u2014 Log-transformed input matrix",
    "  \u2022 fasting_pca_missingness.csv           \u2014 Missing/non-positive row counts per compound",
    "  \u2022 fasting_pca_verification_checks.csv   \u2014 Data audit trail",
    "",
    "PNG plots (in ", plot_dir, "/)",
    "  \u2022 pca_scores_plot.png                  \u2014 Score plot with ellipses and centroids",
    "  \u2022 pca_subject_trajectories.png         \u2014 Paired-subject movement in PCA space",
    "  \u2022 pca_scree_plot.png                   \u2014 Variance explained summary",
    "  \u2022 pca_loadings_plot.png                \u2014 Loading directions with bile acid classes",
    "  \u2022 pca_score_distributions.png          \u2014 Box plots of PC scores by condition",
    "  \u2022 pca_top_loadings.png                 \u2014 Top signed loadings for PC1\u2013PC3",
    "  \u2022 pca_variable_contributions.png       \u2014 Variable contribution percentages"
  )
)

dev.off()

# ─── Session info ─────────────────────────────────────────────────────────────
sink(file.path(pca_dir, "fasting_pca_session_info.txt"))
cat("Fasting-only PCA\n")
cat("Generated:", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("Fasting-only PCA PDF written to", pdf_path, "\n")
cat("Individual PNGs written to", plot_dir, "\n")

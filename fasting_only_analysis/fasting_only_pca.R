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

pdf(pdf_path, width = 11.7, height = 8.3, onefile = TRUE, useDingbats = FALSE)

# Page 1: Title / methods page
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

# Page 2: Score plot
print(score_plot)

# Page 3: Subject trajectories
print(trajectory_plot)

# Page 4: Scree + Loading plot (side by side)
print(scree_plot + loading_plot + plot_layout(widths = c(1, 1.2)))

# Page 5: Score distributions
print(score_box_plot)

# Page 6: Top loadings
print(top_loading_plot)

# Page 7: Variable contributions
print(contrib_plot)

# Page 8: Output files summary
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

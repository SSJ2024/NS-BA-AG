suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

out_dir <- "comparison_work/fixed_csv_report_outputs"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- "IRMA23930_MS-BA_Results.csv"
df <- read.csv(csv_path, check.names = FALSE)

ba_cols <- c("CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA",
             "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA")
glycine_cols <- c("GCA", "GCDCA", "GDCA", "GLCA", "GUDCA")
taurine_cols <- c("TCA", "TCDCA", "TDCA", "TLCA", "TUDCA")
unconjugated_cols <- c("CA", "CDCA", "DCA", "LCA", "UDCA")
total_cols <- c(glycine_cols, taurine_cols, unconjugated_cols)
aggregate_cols <- c("Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated")
all_result_cols <- c(ba_cols, aggregate_cols)

report_order <- c(
  "GCA", "TCA", "GCDCA", "TCDCA", "CA", "CDCA",
  "GDCA", "TDCA", "GLCA", "TLCA", "DCA", "LCA",
  "UDCA", "GUDCA", "TUDCA",
  "Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated"
)

plot_grid_order <- c(
  "GCA", "TCA", "GCDCA", "TCDCA", "CA", "CDCA",
  "GDCA", "TDCA", "GLCA", "TLCA", "DCA", "LCA",
  "UDCA", "GUDCA", "TUDCA"
)

compound_name <- function(comp) {
  names <- c(
    CDCA = "Chenodeoxycholic Acid (CDCA)",
    GCDCA = "Glycochenodeoxycholic Acid (GCDCA)",
    TCDCA = "Taurochenodeoxycholic Acid (TCDCA)",
    UDCA = "Ursodeoxycholic Acid (UDCA)",
    GUDCA = "Glycoursodeoxycholic Acid (GUDCA)",
    TUDCA = "Tauroursodeoxycholic Acid (TUDCA)",
    DCA = "Deoxycholic Acid (DCA)",
    GDCA = "Glycodeoxycholic Acid (GDCA)",
    LCA = "Lithocholic Acid (LCA)",
    GLCA = "Glycolithocholic Acid (GLCA)",
    TLCA = "Taurolithocholic Acid (TLCA)",
    CA = "Cholic Acid (CA)",
    GCA = "Glycocholic Acid (GCA)",
    TCA = "Taurocholic Acid (TCA)",
    TDCA = "Taurodeoxycholic Acid (TDCA)",
    Total_BA = "Total Bile Acids",
    Glycine_Conjugated = "Glycine-Conjugated Bile Acids",
    Taurine_Conjugated = "Taurine-Conjugated Bile Acids",
    Unconjugated = "Unconjugated Bile Acids"
  )
  unname(names[[comp]])
}

compound_class <- function(comp) {
  if (comp == "Total_BA") return("Aggregate Total Pool")
  if (comp == "Glycine_Conjugated") return("Aggregate Glycine Conjugates")
  if (comp == "Taurine_Conjugated") return("Aggregate Taurine Conjugates")
  if (comp == "Unconjugated") return("Aggregate Unconjugated Pool")
  if (comp %in% c("CA", "CDCA")) return("Primary Unconjugated")
  if (comp %in% c("GCA", "TCA", "GCDCA", "TCDCA")) return("Primary Conjugated")
  if (comp %in% c("DCA", "LCA")) return("Secondary Unconjugated")
  if (comp %in% c("GDCA", "TDCA", "GLCA", "TLCA")) return("Secondary Conjugated")
  if (comp %in% c("UDCA", "GUDCA", "TUDCA")) return("Tertiary / Microbial")
  "Bile Acid"
}

sig_star <- function(q) {
  if (is.na(q)) return("NA")
  if (q < 0.001) return("***")
  if (q < 0.01) return("**")
  if (q < 0.05) return("*")
  "ns"
}

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}

safe_rowsum <- function(data, cols) {
  vals <- data[, cols]
  sums <- rowSums(vals, na.rm = TRUE)
  sums[rowSums(!is.na(vals)) == 0] <- NA_real_
  sums
}

adjust_ltr_median <- function(data) {
  adjusted <- data
  for (comp in ba_cols) {
    ltr_idx <- data$sampleType == "ltr"
    plate_medians <- tapply(data[[comp]][ltr_idx], data$plateID[ltr_idx], median, na.rm = TRUE)
    global_median <- median(data[[comp]][ltr_idx], na.rm = TRUE)
    factors <- global_median / plate_medians
    for (plate in names(factors)) {
      sample_idx <- data$sampleType == "sample" & data$plateID == plate
      factor <- factors[[plate]]
      if (is.na(factor) || is.nan(factor) || is.infinite(factor)) factor <- 1.0
      adjusted[[comp]][sample_idx] <- data[[comp]][sample_idx] * factor
    }
  }
  adjusted
}

df_adj <- adjust_ltr_median(df)
df_adj$Total_BA <- safe_rowsum(df_adj, total_cols)
df_adj$Glycine_Conjugated <- safe_rowsum(df_adj, glycine_cols)
df_adj$Taurine_Conjugated <- safe_rowsum(df_adj, taurine_cols)
df_adj$Unconjugated <- safe_rowsum(df_adj, unconjugated_cols)

samples <- df_adj[df_adj$sampleType == "sample", ]
tube_split <- strsplit(as.character(samples$tubeLabel), "-", fixed = TRUE)
samples$subject <- vapply(tube_split, `[`, character(1), 1)
samples$timepoint <- vapply(tube_split, `[`, character(1), 2)
samples <- samples[samples$timepoint %in% c("2", "18", "23", "26"), ]
samples$timepoint <- factor(samples$timepoint, levels = c("2", "18", "23", "26"))

model_rows <- list()
emm_rows <- list()
plot_list <- list()

for (comp in all_result_cols) {
  dat <- samples[!is.na(samples[[comp]]) & samples[[comp]] > 0, ]
  dat$y_log <- log(dat[[comp]])

  model <- lmer(y_log ~ timepoint + (1 | subject), data = dat, REML = TRUE)
  emm <- emmeans(model, ~ timepoint, lmer.df = "satterthwaite")
  pair_df <- as.data.frame(pairs(emm, adjust = "none"))
  emm_ci <- as.data.frame(confint(emm))

  c_2 <- pair_df[pair_df$contrast == "timepoint2 - timepoint23", ]
  c_18 <- pair_df[pair_df$contrast == "timepoint18 - timepoint23", ]
  c_26 <- pair_df[pair_df$contrast == "timepoint23 - timepoint26", ]

  model_rows[[comp]] <- data.frame(
    Compound = comp,
    Compound_Name = compound_name(comp),
    Class = compound_class(comp),
    Result_Family = ifelse(comp %in% ba_cols, "Individual BA", "Aggregate class"),
    N_Obs = nrow(dat),
    N_Subjects = length(unique(dat$subject)),
    LFC_23_vs_2 = -c_2$estimate,
    SE_23_vs_2 = c_2$SE,
    DF_23_vs_2 = c_2$df,
    P_23_vs_2 = c_2$p.value,
    LFC_23_vs_18 = -c_18$estimate,
    SE_23_vs_18 = c_18$SE,
    DF_23_vs_18 = c_18$df,
    P_23_vs_18 = c_18$p.value,
    LFC_23_vs_26 = c_26$estimate,
    SE_23_vs_26 = c_26$SE,
    DF_23_vs_26 = c_26$df,
    P_23_vs_26 = c_26$p.value,
    stringsAsFactors = FALSE
  )

  emm_rows[[comp]] <- data.frame(
    Compound = comp,
    Timepoint = as.character(emm_ci$timepoint),
    N_Obs = as.integer(table(dat$timepoint)[as.character(emm_ci$timepoint)]),
    EMM_Log = emm_ci$emmean,
    SE_Log = emm_ci$SE,
    DF = emm_ci$df,
    CI_Lower_Log = emm_ci$lower.CL,
    CI_Upper_Log = emm_ci$upper.CL,
    EMM = exp(emm_ci$emmean),
    CI_Lower = exp(emm_ci$lower.CL),
    CI_Upper = exp(emm_ci$upper.CL),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, model_rows)
results$Order <- match(results$Compound, report_order)
results <- results[order(results$Order), ]

q_cols <- c("Q_23_vs_2", "Q_23_vs_18", "Q_23_vs_26")
p_cols <- c("P_23_vs_2", "P_23_vs_18", "P_23_vs_26")

for (family in c("Individual BA", "Aggregate class")) {
  idx <- results$Result_Family == family
  p_values <- as.numeric(t(results[idx, p_cols]))
  q_values <- p.adjust(p_values, method = "BH")
  results[idx, q_cols] <- matrix(q_values, ncol = 3, byrow = TRUE)
}

results$Sig_23_vs_2 <- vapply(results$Q_23_vs_2, sig_star, character(1))
results$Sig_23_vs_18 <- vapply(results$Q_23_vs_18, sig_star, character(1))
results$Sig_23_vs_26 <- vapply(results$Q_23_vs_26, sig_star, character(1))
results$Fold_23_vs_2 <- exp(results$LFC_23_vs_2)
results$Fold_23_vs_18 <- exp(results$LFC_23_vs_18)
results$Fold_23_vs_26 <- exp(results$LFC_23_vs_26)

emm_results <- do.call(rbind, emm_rows)
emm_results$Order <- match(emm_results$Compound, report_order)
emm_results$Timepoint_Order <- match(emm_results$Timepoint, c("2", "18", "23", "26"))
emm_results <- emm_results[order(emm_results$Order, emm_results$Timepoint_Order), ]

write.csv(results, file.path(out_dir, "corrected_lmm_results.csv"), row.names = FALSE)
write.csv(emm_results, file.path(out_dir, "corrected_emm_results.csv"), row.names = FALSE)

axis_labels <- c(
  "2" = "Baseline\n(2)",
  "18" = "NS4 Fasting\n(18)",
  "23" = "NS4 After-Meal\n(23)",
  "26" = "RTDS Fasting\n(26)"
)

format_y <- function(x) {
  format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
}

make_plot <- function(comp, compact = FALSE) {
  d <- emm_results[emm_results$Compound == comp, ]
  d$Timepoint <- factor(d$Timepoint, levels = c("2", "18", "23", "26"))
  d <- d[order(d$Timepoint), ]
  fasting_path <- d[d$Timepoint %in% c("2", "18", "26"), ]
  fasting_path$Path <- "Fasting branch"
  meal_path <- d[d$Timepoint %in% c("18", "23", "26"), ]
  meal_path$Path <- "Meal branch"
  path_data <- rbind(fasting_path, meal_path)
  r <- results[results$Compound == comp, ]

  q_text <- paste0(
    "FDR q: 23 vs 2 ", fmt_p(r$Q_23_vs_2), " ", r$Sig_23_vs_2,
    "; 23 vs 18 ", fmt_p(r$Q_23_vs_18), " ", r$Sig_23_vs_18,
    "; 23 vs 26 ", fmt_p(r$Q_23_vs_26), " ", r$Sig_23_vs_26
  )

  title <- if (compact) comp else compound_name(comp)

  ggplot() +
    geom_line(data = path_data, aes(x = Timepoint, y = EMM, group = Path, color = Path), linewidth = 0.95) +
    geom_point(
      data = d,
      aes(x = Timepoint, y = EMM),
      color = "#1F4E79",
      size = ifelse(compact, 1.8, 2.8),
      show.legend = FALSE
    ) +
    geom_point(
      data = d[d$Timepoint == "23", ],
      aes(x = Timepoint, y = EMM),
      color = "#C55A11",
      size = ifelse(compact, 1.8, 2.8),
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = d,
      aes(x = Timepoint, ymin = CI_Lower, ymax = CI_Upper),
      width = 0.12,
      color = "#384B5A",
      linewidth = 0.4
    ) +
    scale_x_discrete(labels = if (compact) c("2" = "2", "18" = "18", "23" = "23", "26" = "26") else axis_labels) +
    scale_y_log10(labels = format_y) +
    scale_color_manual(values = c("Fasting branch" = "#1F4E79", "Meal branch" = "#C55A11")) +
    labs(
      title = title,
      subtitle = if (compact) NULL else q_text,
      x = NULL,
      y = if (compact) NULL else "Adjusted concentration (log10 scale)",
      color = NULL,
      caption = if (compact) NULL else "Divergence: the meal branch starts from 18 and merges at 26; only point 23 is colored orange. Points are LMM EMMs; error bars are 95% CIs."
    ) +
    theme_minimal(base_size = ifelse(compact, 8, 10)) +
    theme(
      plot.title = element_text(face = "bold", size = ifelse(compact, 9, 11), hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 7, color = "#546A7B"),
      legend.position = ifelse(compact, "none", "bottom")
    )
}

for (comp in all_result_cols) {
  p <- make_plot(comp, compact = FALSE)
  plot_path <- file.path(plot_dir, paste0(comp, "_corrected_trajectory.png"))
  ggsave(plot_path, p, width = 6.8, height = 4.6, dpi = 220)
  plot_list[[comp]] <- make_plot(comp, compact = TRUE)
}

grid_plot <- wrap_plots(plot_list[plot_grid_order], ncol = 3, nrow = 5)
ggsave(file.path(plot_dir, "all_15_corrected_trajectories_grid.png"), grid_plot, width = 9.5, height = 12.0, dpi = 220)

agg_grid <- wrap_plots(plot_list[aggregate_cols], ncol = 2, nrow = 2)
ggsave(file.path(plot_dir, "aggregate_corrected_trajectories_grid.png"), agg_grid, width = 8.5, height = 6.2, dpi = 220)

all_missing_rows <- samples[rowSums(!is.na(samples[, ba_cols])) == 0, c("sampleID", "tubeLabel", "plateID", "subject", "timepoint")]
verification <- data.frame(
  Check = c(
    "Raw data rows",
    "Biological sample rows",
    "Longitudinal rows used before compound-specific missing filtering",
    "Subjects",
    "Timepoint 10 excluded",
    "All-missing BA rows preserved as missing",
    "Finite-sample df method",
    "Multiple testing correction",
    "Non-positive modeled values"
  ),
  Result = c(
    nrow(df),
    sum(df$sampleType == "sample"),
    nrow(samples),
    length(unique(samples$subject)),
    "Yes",
    paste0(nrow(all_missing_rows), " row: ", paste(all_missing_rows$tubeLabel, collapse = ", ")),
    "Satterthwaite via emmeans(lmer.df = 'satterthwaite')",
    "BH FDR, separately for 45 individual BA contrasts and 12 aggregate contrasts",
    sum(sapply(all_result_cols, function(comp) sum(samples[[comp]] <= 0, na.rm = TRUE)))
  ),
  stringsAsFactors = FALSE
)
write.csv(verification, file.path(out_dir, "verification_checks.csv"), row.names = FALSE)

sink(file.path(out_dir, "session_info.txt"))
cat("Corrected CSV/R bile acid analysis\n")
cat("Generated:", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("Corrected analysis outputs written to", out_dir, "\n")

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(patchwork)
})

out_dir <- "fasting_only_analysis/outputs"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- "IRMA23930_MS-BA_Results.csv"
df <- read.csv(csv_path, check.names = FALSE)

fasting_levels <- c("2", "18", "26")
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
  ltr_idx <- data$sampleType == "ltr"
  for (comp in ba_cols) {
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

get_pair_row <- function(pair_df, label) {
  pair_df[pair_df$contrast == label, , drop = FALSE]
}

add_pairwise_fdr <- function(data, idx, p_cols, q_cols) {
  p_values <- as.numeric(t(data[idx, p_cols]))
  q_values <- p.adjust(p_values, method = "BH")
  data[idx, q_cols] <- matrix(q_values, ncol = length(q_cols), byrow = TRUE)
  data
}

add_omnibus_fdr <- function(data, idx) {
  data$Q_Timepoint[idx] <- p.adjust(data$P_Timepoint[idx], method = "BH")
  data
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
excluded_ns4_meal_rows <- sum(samples$timepoint == "23", na.rm = TRUE)
samples <- samples[samples$timepoint %in% fasting_levels, ]
samples$timepoint <- factor(samples$timepoint, levels = fasting_levels)

pair_rows <- list()
omnibus_rows <- list()
emm_rows <- list()
plot_list <- list()

for (comp in all_result_cols) {
  dat <- samples[!is.na(samples[[comp]]) & samples[[comp]] > 0, ]
  dat$y_log <- log(dat[[comp]])

  model <- lmer(y_log ~ timepoint + (1 | subject), data = dat, REML = TRUE)
  emm <- emmeans(model, ~ timepoint, lmer.df = "satterthwaite")
  pair_df <- as.data.frame(pairs(emm, adjust = "none"))
  emm_ci <- as.data.frame(confint(emm))
  omnibus_df <- as.data.frame(anova(model, ddf = "Satterthwaite"))
  omnibus_timepoint <- omnibus_df["timepoint", , drop = FALSE]

  c_2_18 <- get_pair_row(pair_df, "timepoint2 - timepoint18")
  c_2_26 <- get_pair_row(pair_df, "timepoint2 - timepoint26")
  c_18_26 <- get_pair_row(pair_df, "timepoint18 - timepoint26")

  pair_rows[[comp]] <- data.frame(
    Compound = comp,
    Compound_Name = compound_name(comp),
    Class = compound_class(comp),
    Result_Family = ifelse(comp %in% ba_cols, "Individual BA", "Aggregate class"),
    N_Obs = nrow(dat),
    N_Subjects = length(unique(dat$subject)),
    LFC_18_vs_2 = -c_2_18$estimate,
    SE_18_vs_2 = c_2_18$SE,
    DF_18_vs_2 = c_2_18$df,
    P_18_vs_2 = c_2_18$p.value,
    LFC_26_vs_2 = -c_2_26$estimate,
    SE_26_vs_2 = c_2_26$SE,
    DF_26_vs_2 = c_2_26$df,
    P_26_vs_2 = c_2_26$p.value,
    LFC_26_vs_18 = -c_18_26$estimate,
    SE_26_vs_18 = c_18_26$SE,
    DF_26_vs_18 = c_18_26$df,
    P_26_vs_18 = c_18_26$p.value,
    stringsAsFactors = FALSE
  )

  omnibus_rows[[comp]] <- data.frame(
    Compound = comp,
    Compound_Name = compound_name(comp),
    Class = compound_class(comp),
    Result_Family = ifelse(comp %in% ba_cols, "Individual BA", "Aggregate class"),
    N_Obs = nrow(dat),
    N_Subjects = length(unique(dat$subject)),
    NumDF_Timepoint = omnibus_timepoint$NumDF,
    DenDF_Timepoint = omnibus_timepoint$DenDF,
    F_Timepoint = omnibus_timepoint[["F value"]],
    P_Timepoint = omnibus_timepoint[["Pr(>F)"]],
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

pair_results <- do.call(rbind, pair_rows)
pair_results$Order <- match(pair_results$Compound, report_order)
pair_results <- pair_results[order(pair_results$Order), ]

pair_p_cols <- c("P_18_vs_2", "P_26_vs_2", "P_26_vs_18")
pair_q_cols <- c("Q_18_vs_2", "Q_26_vs_2", "Q_26_vs_18")
pair_results <- add_pairwise_fdr(
  pair_results,
  pair_results$Result_Family == "Individual BA",
  pair_p_cols,
  pair_q_cols
)
pair_results <- add_pairwise_fdr(
  pair_results,
  pair_results$Result_Family == "Aggregate class",
  pair_p_cols,
  pair_q_cols
)
pair_results$Sig_18_vs_2 <- vapply(pair_results$Q_18_vs_2, sig_star, character(1))
pair_results$Sig_26_vs_2 <- vapply(pair_results$Q_26_vs_2, sig_star, character(1))
pair_results$Sig_26_vs_18 <- vapply(pair_results$Q_26_vs_18, sig_star, character(1))
pair_results$Fold_18_vs_2 <- exp(pair_results$LFC_18_vs_2)
pair_results$Fold_26_vs_2 <- exp(pair_results$LFC_26_vs_2)
pair_results$Fold_26_vs_18 <- exp(pair_results$LFC_26_vs_18)

omnibus_results <- do.call(rbind, omnibus_rows)
omnibus_results$Order <- match(omnibus_results$Compound, report_order)
omnibus_results <- omnibus_results[order(omnibus_results$Order), ]
omnibus_results$Q_Timepoint <- NA_real_
omnibus_results <- add_omnibus_fdr(
  omnibus_results,
  omnibus_results$Result_Family == "Individual BA"
)
omnibus_results <- add_omnibus_fdr(
  omnibus_results,
  omnibus_results$Result_Family == "Aggregate class"
)
omnibus_results$Sig_Timepoint <- vapply(omnibus_results$Q_Timepoint, sig_star, character(1))

emm_results <- do.call(rbind, emm_rows)
emm_results$Order <- match(emm_results$Compound, report_order)
emm_results$Timepoint_Order <- match(emm_results$Timepoint, fasting_levels)
emm_results <- emm_results[order(emm_results$Order, emm_results$Timepoint_Order), ]

write.csv(pair_results, file.path(out_dir, "fasting_pairwise_results.csv"), row.names = FALSE)
write.csv(omnibus_results, file.path(out_dir, "fasting_omnibus_results.csv"), row.names = FALSE)
write.csv(emm_results, file.path(out_dir, "fasting_emm_results.csv"), row.names = FALSE)

axis_labels <- c(
  "2" = "Baseline\n(2)",
  "18" = "NS4 Fasting\n(18)",
  "26" = "RTDS Fasting\n(26)"
)

format_y <- function(x) {
  format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
}

make_plot <- function(comp, compact = FALSE) {
  d <- emm_results[emm_results$Compound == comp, ]
  d$Timepoint <- factor(d$Timepoint, levels = fasting_levels)
  d <- d[order(d$Timepoint), ]
  r <- pair_results[pair_results$Compound == comp, ]
  title <- if (compact) comp else compound_name(comp)
  q_text <- paste0(
    "FDR q: 18 vs 2 ", fmt_p(r$Q_18_vs_2), " ", r$Sig_18_vs_2,
    "; 26 vs 2 ", fmt_p(r$Q_26_vs_2), " ", r$Sig_26_vs_2,
    "; 26 vs 18 ", fmt_p(r$Q_26_vs_18), " ", r$Sig_26_vs_18
  )

  ggplot(d, aes(x = Timepoint, y = EMM, group = 1)) +
    geom_line(color = "#1F4E79", linewidth = 0.95) +
    geom_point(color = "#1F4E79", size = ifelse(compact, 1.8, 2.8)) +
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper),
      width = 0.12,
      color = "#384B5A",
      linewidth = 0.4
    ) +
    scale_x_discrete(labels = if (compact) c("2" = "2", "18" = "18", "26" = "26") else axis_labels) +
    scale_y_log10(labels = format_y) +
    labs(
      title = title,
      subtitle = if (compact) NULL else q_text,
      x = NULL,
      y = if (compact) NULL else "Adjusted concentration (log10 scale)",
      caption = if (compact) NULL else "Fasting-only LMM EMMs with 95% CIs. Timepoint 23 is excluded before fitting."
    ) +
    theme_minimal(base_size = ifelse(compact, 8, 10)) +
    theme(
      plot.title = element_text(face = "bold", size = ifelse(compact, 9, 11), hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 7, color = "#546A7B")
    )
}

for (comp in all_result_cols) {
  p <- make_plot(comp, compact = FALSE)
  ggsave(
    file.path(plot_dir, paste0(comp, "_fasting_trajectory.png")),
    p,
    width = 6.8,
    height = 4.6,
    dpi = 220
  )
  plot_list[[comp]] <- make_plot(comp, compact = TRUE)
}

ba_grid <- wrap_plots(plot_list[plot_grid_order], ncol = 3, nrow = 5)
ggsave(
  file.path(plot_dir, "all_15_fasting_trajectories_grid.png"),
  ba_grid,
  width = 9.5,
  height = 12.0,
  dpi = 220
)

aggregate_grid <- wrap_plots(plot_list[aggregate_cols], ncol = 2, nrow = 2)
ggsave(
  file.path(plot_dir, "aggregate_fasting_trajectories_grid.png"),
  aggregate_grid,
  width = 8.5,
  height = 6.2,
  dpi = 220
)

all_missing_rows <- samples[rowSums(!is.na(samples[, ba_cols])) == 0, c("sampleID", "tubeLabel", "plateID", "subject", "timepoint")]
verification <- data.frame(
  Check = c(
    "Raw data rows",
    "Biological sample rows",
    "Fasting rows used before compound-specific missing filtering",
    "Subjects",
    "Timepoints included",
    "NS4 after-meal rows excluded before model fitting",
    "Timepoint 10 excluded",
    "All-missing BA rows preserved as missing",
    "Finite-sample df method",
    "Pairwise multiple testing correction",
    "Omnibus multiple testing correction",
    "Non-positive modeled values"
  ),
  Result = c(
    nrow(df),
    sum(df$sampleType == "sample"),
    nrow(samples),
    length(unique(samples$subject)),
    paste(fasting_levels, collapse = ", "),
    excluded_ns4_meal_rows,
    "Yes",
    paste0(nrow(all_missing_rows), " row: ", paste(all_missing_rows$tubeLabel, collapse = ", ")),
    "Satterthwaite via lmerTest and emmeans(lmer.df = 'satterthwaite')",
    "BH FDR, separately for 45 individual BA fasting contrasts and 12 aggregate fasting contrasts",
    "BH FDR, separately for 15 individual BA timepoint tests and 4 aggregate timepoint tests",
    sum(sapply(all_result_cols, function(comp) sum(samples[[comp]] <= 0, na.rm = TRUE)))
  ),
  stringsAsFactors = FALSE
)
write.csv(verification, file.path(out_dir, "verification_checks.csv"), row.names = FALSE)

sink(file.path(out_dir, "session_info.txt"))
cat("Fasting-only CSV/R bile acid analysis\n")
cat("Generated:", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("Fasting-only analysis outputs written to", out_dir, "\n")

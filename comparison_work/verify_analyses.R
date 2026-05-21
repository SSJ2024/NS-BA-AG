library(lme4)
library(lmerTest)
library(emmeans)
library(readxl)

ba_cols <- c("CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA",
             "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA")

glycine_cols <- c("GCA", "GCDCA", "GDCA", "GLCA", "GUDCA")
taurine_cols <- c("TCA", "TCDCA", "TDCA", "TLCA", "TUDCA")
unconjugated_cols <- c("CA", "CDCA", "DCA", "LCA", "UDCA")

star <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

adjust_ltr <- function(df, stat = c("mean", "median")) {
  stat <- match.arg(stat)
  out <- df
  ltr <- df[df$sampleType == "ltr", ]
  for (comp in ba_cols) {
    plate_stat <- tapply(ltr[[comp]], ltr$plateID, if (stat == "mean") mean else median, na.rm = TRUE)
    global_stat <- (if (stat == "mean") mean else median)(ltr[[comp]], na.rm = TRUE)
    factors <- global_stat / plate_stat
    for (plate in names(factors)) {
      idx <- df$sampleType == "sample" & df$plateID == plate
      out[[comp]][idx] <- df[[comp]][idx] * factors[[plate]]
    }
  }
  out
}

prep_samples <- function(df, levels) {
  samples <- df[df$sampleType == "sample", ]
  split <- strsplit(as.character(samples$tubeLabel), "-", fixed = TRUE)
  samples$subject <- vapply(split, `[`, character(1), 1)
  samples$timepoint <- vapply(split, `[`, character(1), 2)
  samples <- samples[samples$timepoint %in% levels, ]
  samples$timepoint <- factor(samples$timepoint, levels = levels)
  samples
}

csv <- read.csv("IRMA23930_MS-BA_Results.csv", check.names = FALSE)
xlsx_results <- read_excel("Bile_Acid_LMM_Results.xlsx", sheet = "LMM_Results")

mean_adj <- adjust_ltr(csv, "mean")
mean_samples <- prep_samples(mean_adj, c("2", "18", "26", "23"))

lmm_rows <- list()
for (comp in ba_cols) {
  dat <- mean_samples[!is.na(mean_samples[[comp]]) & mean_samples[[comp]] > 0, ]
  dat$y <- log10(dat[[comp]])
  model <- lmer(y ~ timepoint + (1 | subject), data = dat, REML = TRUE)

  emm_satt <- emmeans(model, ~ timepoint, lmer.df = "satterthwaite")
  post_satt <- as.data.frame(contrast(
    emm_satt,
    list(post_vs_fasting_avg = c(-1/3, -1/3, -1/3, 1)),
    adjust = "none"
  ))

  emm_kr <- emmeans(model, ~ timepoint, lmer.df = "kenward-roger")
  post_kr <- as.data.frame(contrast(
    emm_kr,
    list(post_vs_fasting_avg = c(-1/3, -1/3, -1/3, 1)),
    adjust = "none"
  ))

  emm_asym <- emmeans(model, ~ timepoint, lmer.df = "asymptotic")
  post_asym <- as.data.frame(contrast(
    emm_asym,
    list(post_vs_fasting_avg = c(-1/3, -1/3, -1/3, 1)),
    adjust = "none"
  ))

  x <- xlsx_results[xlsx_results$Compound == comp, ]
  lmm_rows[[comp]] <- data.frame(
    Compound = comp,
    N_Obs = nrow(dat),
    Xlsx_Post_p = as.numeric(x$`Post p`),
    R_Asymptotic_Post_p = post_asym$p.value,
    R_Satterthwaite_Post_p = post_satt$p.value,
    R_KenwardRoger_Post_p = post_kr$p.value,
    Xlsx_Fold_23_vs_Fasting = as.numeric(x$`Fold Change 23 / Fasting Avg`),
    R_Satterthwaite_log10_diff = post_satt$estimate,
    R_Satterthwaite_Fold = 10^post_satt$estimate,
    stringsAsFactors = FALSE
  )
}
lmm_check <- do.call(rbind, lmm_rows)
lmm_check$Xlsx_Post_FDR_q <- p.adjust(lmm_check$Xlsx_Post_p, method = "BH")
lmm_check$R_Asymptotic_FDR_q <- p.adjust(lmm_check$R_Asymptotic_Post_p, method = "BH")
lmm_check$R_Satterthwaite_FDR_q <- p.adjust(lmm_check$R_Satterthwaite_Post_p, method = "BH")
lmm_check$R_KenwardRoger_FDR_q <- p.adjust(lmm_check$R_KenwardRoger_Post_p, method = "BH")
lmm_check$Xlsx_Sig <- vapply(lmm_check$Xlsx_Post_FDR_q, star, character(1))
lmm_check$R_Satterthwaite_Sig <- vapply(lmm_check$R_Satterthwaite_FDR_q, star, character(1))
write.csv(lmm_check, "comparison_work/lmm_composite_check.csv", row.names = FALSE)

median_adj <- adjust_ltr(csv, "median")

safe_rowsum <- function(df, cols) {
  vals <- df[, cols]
  sums <- rowSums(vals, na.rm = TRUE)
  sums[rowSums(!is.na(vals)) == 0] <- NA_real_
  sums
}

median_adj_report_style <- median_adj
median_adj_report_style$Total_BA <- rowSums(median_adj_report_style[, c(glycine_cols, taurine_cols, unconjugated_cols)], na.rm = TRUE)
median_adj_report_style$Glycine_Conjugated <- rowSums(median_adj_report_style[, glycine_cols], na.rm = TRUE)
median_adj_report_style$Taurine_Conjugated <- rowSums(median_adj_report_style[, taurine_cols], na.rm = TRUE)
median_adj_report_style$Unconjugated <- rowSums(median_adj_report_style[, unconjugated_cols], na.rm = TRUE)

median_adj_corrected <- median_adj
median_adj_corrected$Total_BA <- safe_rowsum(median_adj_corrected, c(glycine_cols, taurine_cols, unconjugated_cols))
median_adj_corrected$Glycine_Conjugated <- safe_rowsum(median_adj_corrected, glycine_cols)
median_adj_corrected$Taurine_Conjugated <- safe_rowsum(median_adj_corrected, taurine_cols)
median_adj_corrected$Unconjugated <- safe_rowsum(median_adj_corrected, unconjugated_cols)

all_pair_cols <- c(ba_cols, "Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated")

run_pairwise <- function(source_df, output_path, report_style = FALSE, lmer_df = "kenward-roger") {
  pair_samples <- prep_samples(source_df, c("2", "18", "23", "26"))
  pair_rows <- list()
  for (comp in all_pair_cols) {
    if (report_style) {
      dat <- pair_samples[!is.na(pair_samples[[comp]]), ]
      offset <- if (min(dat[[comp]], na.rm = TRUE) <= 0) 1 else 0
      dat$y <- log(dat[[comp]] + offset)
    } else {
      dat <- pair_samples[!is.na(pair_samples[[comp]]) & pair_samples[[comp]] > 0, ]
      dat$y <- log(dat[[comp]])
    }
    model <- lmer(y ~ timepoint + (1 | subject), data = dat, REML = TRUE)
    pairs_all <- as.data.frame(pairs(emmeans(model, ~ timepoint, lmer.df = lmer_df), adjust = "none"))

    c_2 <- pairs_all[pairs_all$contrast == "timepoint2 - timepoint23", ]
    c_18 <- pairs_all[pairs_all$contrast == "timepoint18 - timepoint23", ]
    c_26 <- pairs_all[pairs_all$contrast == "timepoint23 - timepoint26", ]

    pair_rows[[comp]] <- data.frame(
      Compound = comp,
      N_Obs = nrow(dat),
      LFC_23_vs_2 = -c_2$estimate,
      P_23_vs_2 = c_2$p.value,
      LFC_23_vs_18 = -c_18$estimate,
      P_23_vs_18 = c_18$p.value,
      LFC_23_vs_26 = c_26$estimate,
      P_23_vs_26 = c_26$p.value,
      stringsAsFactors = FALSE
    )
  }
  pair_check <- do.call(rbind, pair_rows)

  ba_idx <- pair_check$Compound %in% ba_cols
  ba_p <- as.vector(t(pair_check[ba_idx, c("P_23_vs_2", "P_23_vs_18", "P_23_vs_26")]))
  ba_q <- p.adjust(ba_p, method = "BH")
  pair_check[ba_idx, c("Q_23_vs_2", "Q_23_vs_18", "Q_23_vs_26")] <- matrix(ba_q, ncol = 3, byrow = TRUE)

  agg_idx <- !ba_idx
  agg_p <- as.vector(t(pair_check[agg_idx, c("P_23_vs_2", "P_23_vs_18", "P_23_vs_26")]))
  agg_q <- p.adjust(agg_p, method = "BH")
  pair_check[agg_idx, c("Q_23_vs_2", "Q_23_vs_18", "Q_23_vs_26")] <- matrix(agg_q, ncol = 3, byrow = TRUE)

  for (nm in c("Q_23_vs_2", "Q_23_vs_18", "Q_23_vs_26")) {
    pair_check[[sub("^Q_", "Sig_", nm)]] <- vapply(pair_check[[nm]], star, character(1))
  }
  write.csv(pair_check, output_path, row.names = FALSE)
}

run_pairwise(median_adj_report_style, "comparison_work/r_pairwise_report_style_check.csv", report_style = TRUE, lmer_df = "kenward-roger")
run_pairwise(median_adj_corrected, "comparison_work/r_pairwise_corrected_missing_check.csv", report_style = FALSE, lmer_df = "satterthwaite")

cat("Wrote comparison_work/lmm_composite_check.csv, comparison_work/r_pairwise_report_style_check.csv, and comparison_work/r_pairwise_corrected_missing_check.csv\n")

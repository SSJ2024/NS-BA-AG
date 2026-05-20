# Complete Bile Acid Analysis, Visualization, and Word Report Generation

library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(officer)
library(flextable)
library(patchwork)

# Create plots folder if it doesn't exist
dir.create("plots", showWarnings = FALSE)

# ---------------------------------------------------------
# 1. Load and Clean Data
# ---------------------------------------------------------
csv_path <- "IRMA23930_MS-BA_Results.csv"
df <- read.csv(csv_path)

# Separate biological samples from reference samples
samples_raw <- df[df$sampleType == "sample", ]
ltrs_raw <- df[df$sampleType == "ltr", ]

# Parse tube label to Subject and Timepoint
tube_split <- strsplit(as.character(samples_raw$tubeLabel), "-")
samples_raw$subject <- sapply(tube_split, function(x) x[1])
samples_raw$timepoint <- sapply(tube_split, function(x) x[2])

# Filter for the relevant timepoints (excluding 10, keeping 2, 18, 23, 26)
# 2 = Baseline, 18 = NS4 Fasting, 23 = NS4 After-Meal, 26 = RTDS Fasting
samples_sub <- samples_raw[samples_raw$timepoint %in% c("2", "18", "23", "26"), ]
samples_sub$timepoint <- factor(samples_sub$timepoint, levels = c("2", "18", "23", "26"))

# The 15 bile acids to analyze
ba_cols <- c("CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA", 
             "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA")

# Define aggregate classes
glycine_cols <- c("GCA", "GCDCA", "GDCA", "GLCA", "GUDCA")
taurine_cols <- c("TCA", "TCDCA", "TDCA", "TLCA", "TUDCA")
unconjugated_cols <- c("CA", "CDCA", "DCA", "LCA", "UDCA")
total_cols <- c(glycine_cols, taurine_cols, unconjugated_cols)
comp_cols <- c("Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated")

# Sort bile acids by class for logical presentation and grid ordering
# Class ordering: Conjugated Primaries, Primaries, Conjugated Secondaries, Secondaries, Tertiaries
ba_order <- c("GCA", "TCA", "GCDCA", "TCDCA", "CA", "CDCA", "GDCA", "TDCA", "GLCA", "TLCA", "DCA", "LCA", "UDCA", "GUDCA", "TUDCA")

# ---------------------------------------------------------
# 2. LTR-based Plate Normalization
# ---------------------------------------------------------
df_adj <- df
normalization_factors <- list()

for (comp in ba_cols) {
  # Calculate plate-specific medians for LTR samples
  plate_medians <- tapply(df[[comp]][df$sampleType == "ltr"], 
                          df$plateID[df$sampleType == "ltr"], 
                          median, na.rm = TRUE)
  
  # Calculate global median of LTR samples
  global_median <- median(df[[comp]][df$sampleType == "ltr"], na.rm = TRUE)
  
  factors <- global_median / plate_medians
  normalization_factors[[comp]] <- factors
  
  # Adjust concentrations for all samples
  for (plate in names(plate_medians)) {
    factor <- factors[plate]
    if (is.na(factor) || is.nan(factor) || is.infinite(factor)) {
      factor <- 1.0
    }
    df_adj[[comp]][df$plateID == plate] <- df[[comp]][df$plateID == plate] * factor
  }
}

# Extract the adjusted biological samples
df_adj$Total_BA <- rowSums(df_adj[, total_cols], na.rm = TRUE)
df_adj$Glycine_Conjugated <- rowSums(df_adj[, glycine_cols], na.rm = TRUE)
df_adj$Taurine_Conjugated <- rowSums(df_adj[, taurine_cols], na.rm = TRUE)
df_adj$Unconjugated <- rowSums(df_adj[, unconjugated_cols], na.rm = TRUE)

samples_adj <- df_adj[df_adj$sampleType == "sample", ]
tube_split_adj <- strsplit(as.character(samples_adj$tubeLabel), "-")
samples_adj$subject <- sapply(tube_split_adj, function(x) x[1])
samples_adj$timepoint <- sapply(tube_split_adj, function(x) x[2])

# Filter for the relevant timepoints in adjusted data
samples_sub_adj <- samples_adj[samples_adj$timepoint %in% c("2", "18", "23", "26"), ]
samples_sub_adj$timepoint <- factor(samples_sub_adj$timepoint, levels = c("2", "18", "23", "26"))

# ---------------------------------------------------------
# Helper function to get significance stars
# ---------------------------------------------------------
get_stars <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# ---------------------------------------------------------
# Helper function to get full name of Bile Acid
# ---------------------------------------------------------
get_ba_full_name <- function(comp) {
  names <- list(
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
    Total_BA = "Total Bile Acids (Total_BA)",
    Glycine_Conjugated = "Glycine-Conjugated Bile Acids",
    Taurine_Conjugated = "Taurine-Conjugated Bile Acids",
    Unconjugated = "Unconjugated Bile Acids"
  )
  return(names[[comp]])
}

# ---------------------------------------------------------
# Helper function to get biological category
# ---------------------------------------------------------
get_ba_category <- function(comp) {
  primaries <- c("CA", "CDCA")
  conjugated_primaries <- c("GCA", "TCA", "GCDCA", "TCDCA")
  secondaries <- c("DCA", "LCA")
  conjugated_secondaries <- c("GDCA", "TDCA", "GLCA", "TLCA")
  tertiaries <- c("UDCA", "GUDCA", "TUDCA")
  
  if (comp == "Total_BA") return("Aggregate Total Bile Acid Pool")
  if (comp == "Glycine_Conjugated") return("Aggregate Conjugated Bile Acids (Glycine Class)")
  if (comp == "Taurine_Conjugated") return("Aggregate Conjugated Bile Acids (Taurine Class)")
  if (comp == "Unconjugated") return("Aggregate Unconjugated Bile Acids")
  
  if (comp %in% primaries) return("Primary Unconjugated Bile Acid")
  if (comp %in% conjugated_primaries) return("Primary Conjugated Bile Acid")
  if (comp %in% secondaries) return("Secondary Unconjugated Bile Acid")
  if (comp %in% conjugated_secondaries) return("Secondary Conjugated Bile Acid")
  if (comp %in% tertiaries) return("Tertiary / Microbial-derived Bile Acid")
  return("Bile Acid")
}

# ---------------------------------------------------------
# Helper function to format p-values cleanly
# ---------------------------------------------------------
format_p_value <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) return("< 0.001")
  return(sprintf("%.3f", p))
}

# ---------------------------------------------------------
# Helper function to get biological category (short form)
# ---------------------------------------------------------
get_ba_category_short <- function(comp) {
  primaries <- c("CA", "CDCA")
  conjugated_primaries <- c("GCA", "TCA", "GCDCA", "TCDCA")
  secondaries <- c("DCA", "LCA")
  conjugated_secondaries <- c("GDCA", "TDCA", "GLCA", "TLCA")
  tertiaries <- c("UDCA", "GUDCA", "TUDCA")
  
  if (comp == "Total_BA") return("Total Pool")
  if (comp == "Glycine_Conjugated") return("Gly Conjugates")
  if (comp == "Taurine_Conjugated") return("Tau Conjugates")
  if (comp == "Unconjugated") return("Unconjugated Pool")
  
  if (comp %in% primaries) return("Primary Unconj")
  if (comp %in% conjugated_primaries) return("Primary Conj")
  if (comp %in% secondaries) return("Secondary Unconj")
  if (comp %in% conjugated_secondaries) return("Secondary Conj")
  if (comp %in% tertiaries) return("Tertiary / Microbial")
  return("Bile Acid")
}

# ---------------------------------------------------------
# Helper function to generate dynamic biological description
# ---------------------------------------------------------
generate_bio_desc <- function(comp, s) {
  # Static biological background for each compound
  bg <- switch(comp,
    GCA = "Glycocholic Acid (GCA) is a major primary conjugated bile acid synthesized in hepatocytes by the conjugation of cholic acid with glycine. Secreted actively into bile, GCA is released into the duodenum following meal ingestion via cholecystokinin (CCK)-mediated gallbladder contraction. GCA's enterohepatic recovery is facilitated by active transport via the ASBT transporter in the ileum.",
    TCA = "Taurocholic Acid (TCA) is a primary taurine-conjugated bile acid synthesized in hepatocytes. Because taurine conjugation provides a highly polar structure, TCA is restricted to the intestinal lumen until it reaches the active transport systems of the terminal ileum. During fasting, systemic TCA concentrations remain extremely low, reflecting quiescent enterohepatic circulation between meals.",
    GCDCA = "Glycochenodeoxycholic Acid (GCDCA) is the glycine conjugate of chenodeoxycholic acid and represents one of the most abundant bile acids in human bile. GCDCA plays a central role in postprandial fat emulsification and is actively reabsorbed via ileal transporters during enterohepatic circulation.",
    TCDCA = "Taurochenodeoxycholic Acid (TCDCA) is the taurine conjugate of chenodeoxycholic acid. TCDCA is a major signaling ligand for the nuclear farnesoid X receptor (FXR). Postprandial surges of TCDCA are crucial for triggering hepatic FXR signaling to suppress de novo bile acid synthesis via FGF15/19 feedback.",
    CA = "Cholic Acid (CA) is an unconjugated primary bile acid synthesized directly from cholesterol via the classical pathway in hepatocytes. Unlike its conjugated forms GCA and TCA, CA is less polar and can undergo passive absorption in the jejunum and colon, though active transport in the ileum is still dominant. A portion of the circulating CA pool arises from microbial deconjugation of GCA and TCA in the distal gut.",
    CDCA = "Chenodeoxycholic Acid (CDCA) is an unconjugated primary bile acid and the most potent endogenous agonist for the farnesoid X receptor (FXR), playing a foundational role in metabolic regulation. Since CDCA is highly lipophilic, it is easily absorbed passively in the intestine.",
    GDCA = "Glycodeoxycholic Acid (GDCA) is a secondary conjugated bile acid formed when gut bacteria dehydroxylate cholic acid derivatives in the colon, which are then conjugated with glycine, reabsorbed, and recirculated through enterohepatic cycling.",
    TDCA = "Taurodeoxycholic Acid (TDCA) is the taurine conjugate of the secondary bile acid deoxycholic acid. Formed via microbial dehydroxylation and subsequent hepatic taurine conjugation, TDCA actively participates in enterohepatic recirculation. As a potent TGR5 receptor agonist, systemic TDCA can stimulate GLP-1 release from intestinal L-cells to optimize glucose-stimulated insulin secretion.",
    GLCA = "Glycolithocholic Acid (GLCA) is a glycine-conjugated secondary bile acid derived from lithocholic acid (LCA), which is highly hepatotoxic. To prevent toxicity, the liver efficiently sulfates and conjugates LCA. Systemic concentrations of GLCA are very low because the LCA pool is kept extremely small via rapid fecal excretion and sulfation.",
    TLCA = "Taurolithocholic Acid (TLCA) is the taurine conjugate of the secondary bile acid lithocholic acid. Like GLCA, it is highly hydrophobic and potentially toxic, requiring tight regulation. Systemic circulation of TLCA is minimal.",
    DCA = "Deoxycholic Acid (DCA) is a major unconjugated secondary bile acid produced via microbial 7-alpha-dehydroxylation of cholic acid in the large intestine. Unconjugated BAs typically require microbial deconjugation and show slower absorption kinetics than conjugated forms.",
    LCA = "Lithocholic Acid (LCA) is an unconjugated secondary bile acid formed by bacterial 7-alpha-dehydroxylation of CDCA in the colon. Extremely hydrophobic and potentially toxic, LCA is mostly excreted in feces or sulfated in hepatocytes. Systemic concentrations are extremely low.",
    UDCA = "Ursodeoxycholic Acid (UDCA) is an unconjugated tertiary bile acid formed via bacterial epimerization of CDCA. UDCA has well-known hydrophilic, cytoprotective properties and its presence helps buffer the cytolytic potential of hydrophobic secondary bile acids.",
    GUDCA = "Glycoursodeoxycholic Acid (GUDCA) is the glycine conjugate of the tertiary bile acid ursodeoxycholic acid. Known for its highly hydrophilic and FXR-antagonistic or cytoprotective properties, GUDCA contributes to the hydrophilic, protective bile acid pool.",
    TUDCA = "Tauroursodeoxycholic Acid (TUDCA) is the taurine conjugate of the tertiary bile acid ursodeoxycholic acid. TUDCA is highly hydrophilic and acts as a molecular chaperone, reducing endoplasmic reticulum (ER) stress in hepatocytes.",
    Total_BA = "The Total Bile Acid pool is a vital physiological biomarker representing the overall capacity of enterohepatic circulation. Meal-induced gallbladder contraction, mediated by cholecystokinin (CCK), triggers release of the total pool into the duodenum, followed by active reabsorption in the terminal ileum.",
    Glycine_Conjugated = "Glycine conjugation is the predominant pathway for bile acid conjugation in humans. Conjugated bile acids are fully ionized at physiological intestinal pH, restricting passive diffusion across cell membranes and confining them to the intestinal lumen until they reach the active transporters in the terminal ileum.",
    Taurine_Conjugated = "Taurine conjugation is a highly polar conjugation pathway. Due to their ultra-low pKa (~1.5), taurine conjugates remain fully negatively charged throughout the entire small intestine, preventing passive absorption and making them dependent on active sodium-coupled transporters (ASBT) in the terminal ileum. As potent agonists of the TGR5 receptor, postprandial taurine conjugate surges serve as systemic signaling cues to trigger GLP-1 release and regulate energy expenditure.",
    Unconjugated = "Unconjugated bile acids (CA, CDCA, DCA, LCA, UDCA) are formed when primary conjugated bile acids undergo deconjugation by gut bacterial bile salt hydrolases (BSH) in the distal ileum and colon. Unconjugated bile acids are highly lipophilic and can be absorbed via passive non-ionic diffusion throughout the intestine, though their appearance in plasma is partially delayed due to the necessity of microbial BSH action.",
    "No description available."
  )
  
  # Dynamic statistical interpretation based on FDR-corrected q-values
  q_2 <- s$q_23_vs_2
  q_18 <- s$q_23_vs_18
  q_26 <- s$q_23_vs_26
  sig_2 <- !is.na(q_2) && q_2 < 0.05
  sig_18 <- !is.na(q_18) && q_18 < 0.05
  sig_26 <- !is.na(q_26) && q_26 < 0.05
  n_sig <- sum(sig_2, sig_18, sig_26)
  
  if (n_sig == 0) {
    stats_text <- paste0(
      "The LMM analysis did not reveal statistically significant postprandial changes ",
      "at Timepoint 23 relative to any of the three fasting states after FDR correction ",
      "(all q > 0.05). This indicates that systemic concentrations of this bile acid remain ",
      "relatively stable across meal and fasting conditions in this cohort (N = 32 subjects)."
    )
  } else {
    sig_parts <- c()
    if (sig_2) {
      dir <- ifelse(s$est_23_vs_2 > 0, "elevated", "decreased")
      sig_parts <- c(sig_parts, paste0(
        "vs. Baseline (FDR q = ", signif(q_2, 3), ", ", get_stars(q_2), ", ", dir, ")"
      ))
    }
    if (sig_18) {
      dir <- ifelse(s$est_23_vs_18 > 0, "elevated", "decreased")
      sig_parts <- c(sig_parts, paste0(
        "vs. NS4 Fasting (FDR q = ", signif(q_18, 3), ", ", get_stars(q_18), ", ", dir, ")"
      ))
    }
    if (sig_26) {
      dir <- ifelse(s$est_23_vs_26 > 0, "elevated", "decreased")
      sig_parts <- c(sig_parts, paste0(
        "vs. RTDS Fasting (FDR q = ", signif(q_26, 3), ", ", get_stars(q_26), ", ", dir, ")"
      ))
    }
    stats_text <- paste0(
      "The LMM analysis reveals statistically significant postprandial changes at Timepoint 23 ",
      paste(sig_parts, collapse = "; "), "."
    )
    
    # Note non-significant contrasts
    nsig_parts <- c()
    if (!sig_2) nsig_parts <- c(nsig_parts, "Baseline")
    if (!sig_18) nsig_parts <- c(nsig_parts, "NS4 Fasting")
    if (!sig_26) nsig_parts <- c(nsig_parts, "RTDS Fasting")
    if (length(nsig_parts) > 0) {
      stats_text <- paste0(
        stats_text,
        " The contrast(s) vs. ", paste(nsig_parts, collapse = " and "),
        " did not reach significance after FDR correction."
      )
    }
  }
  
  paste(bg, stats_text)
}

# ---------------------------------------------------------
# 3. Statistical Analysis & Plotting Loop
# ---------------------------------------------------------
stats_results <- list()
grid_plots <- list()

for (comp in c(ba_cols, comp_cols)) {
  cat("\nProcessing:", comp, "...\n")
  
  # Fit log-transformed LMM
  # Using natural log: log(value). We handle zeroes by adding 1 if necessary (though BA data is mostly >0)
  min_val <- min(samples_sub_adj[[comp]], na.rm = TRUE)
  offset <- if (min_val <= 0) 1 else 0
  
  # Prepare subset data without NA for this compound
  sub_data <- samples_sub_adj[!is.na(samples_sub_adj[[comp]]), ]
  if (offset > 0) {
    sub_data$y_log <- log(sub_data[[comp]] + offset)
  } else {
    sub_data$y_log <- log(sub_data[[comp]])
  }
  
  # Fit LMM
  model <- lmer(y_log ~ timepoint + (1 | subject), data = sub_data)
  
  # Perform pairwise comparisons against after-meal state (23)
  em <- emmeans(model, ~ timepoint)
  pairs_all <- as.data.frame(pairs(em, adjust = "none"))
  
  # Extract specific comparisons of interest: 23 vs 2, 23 vs 18, 23 vs 26
  comp_23_vs_2 <- pairs_all[pairs_all$contrast == "timepoint2 - timepoint23", ]
  comp_23_vs_18 <- pairs_all[pairs_all$contrast == "timepoint18 - timepoint23", ]
  comp_23_vs_26 <- pairs_all[pairs_all$contrast == "timepoint23 - timepoint26", ]
  
  # Save statistical values
  stats_results[[comp]] <- list(
    p_23_vs_2 = comp_23_vs_2$p.value,
    p_23_vs_18 = comp_23_vs_18$p.value,
    p_23_vs_26 = comp_23_vs_26$p.value,
    est_23_vs_2 = -comp_23_vs_2$estimate, # Invert to represent 23 - 2
    est_23_vs_18 = -comp_23_vs_18$estimate, # Invert to represent 23 - 18
    est_23_vs_26 = comp_23_vs_26$estimate, # Already represents 23 - 26 (no inversion needed)
    se_23_vs_2 = comp_23_vs_2$SE,
    se_23_vs_18 = comp_23_vs_18$SE,
    se_23_vs_26 = comp_23_vs_26$SE
  )
  
  # ---------------------------------------------------------
  # Prepare Plot Data (Fasting vs Meal Trajectories)
  # ---------------------------------------------------------
  # We construct two lines:
  # Fasting: 2 -> 18 -> 26
  # After-Meal: 2 -> 23
  # Calculate mean and SEM on the raw concentration scale for plotting
  summary_stats <- sub_data %>%
    group_by(timepoint) %>%
    summarise(
      mean_val = mean(.data[[comp]], na.rm = TRUE),
      sd_val = sd(.data[[comp]], na.rm = TRUE),
      n_val = sum(!is.na(.data[[comp]])),
      sem_val = sd_val / sqrt(n_val),
      .groups = "drop"
    )
  
  # Create Fasting Line (connects Baseline [2] -> NS4 Fasting [18] -> RTDS Fasting [26])
  fasting_line <- summary_stats %>%
    filter(timepoint %in% c("2", "18", "26")) %>%
    mutate(Group = "Fasting")
  
  # Create Meal Line (starts at Baseline [2] and goes directly to After-Meal [23], bypassing NS4 Fasting [18])
  meal_line <- summary_stats %>%
    filter(timepoint %in% c("2", "23")) %>%
    mutate(Group = "After-Meal")
  
  # Create Dotted Orange Recovery Line (connects NS4 After-Meal [23] to RTDS Fasting [26])
  recovery_line <- summary_stats %>%
    filter(timepoint %in% c("23", "26"))
  
  plot_data <- rbind(fasting_line, meal_line)
  plot_data$Group <- factor(plot_data$Group, levels = c("Fasting", "After-Meal"))
  
  # ---------------------------------------------------------
  # Generate ggplot2 Plot
  # ---------------------------------------------------------
  # Find appropriate y-limits for significance brackets
  max_mean_sem <- max(plot_data$mean_val + plot_data$sem_val)
  min_mean_sem <- min(plot_data$mean_val - plot_data$sem_val)
  
  # Use log10 scale
  y_max_log <- log10(max_mean_sem)
  y_min_log <- log10(min_mean_sem)
  y_range_log <- y_max_log - y_min_log
  
  # Calculate heights for cascading brackets on the log10 scale
  bracket_y1 <- 10^(y_max_log + y_range_log * 0.15) # 23 vs 26
  bracket_y2 <- 10^(y_max_log + y_range_log * 0.35) # 23 vs 18
  bracket_y3 <- 10^(y_max_log + y_range_log * 0.55) # 23 vs 2
  y_limit_max <- 10^(y_max_log + y_range_log * 0.75)
  
  # Format p-values/stars for brackets
  p_2_star <- paste0(get_stars(comp_23_vs_2$p.value), " (p = ", format.pval(comp_23_vs_2$p.value, digits=3), ")")
  p_18_star <- paste0(get_stars(comp_23_vs_18$p.value), " (p = ", format.pval(comp_23_vs_18$p.value, digits=3), ")")
  p_26_star <- paste0(get_stars(comp_23_vs_26$p.value), " (p = ", format.pval(comp_23_vs_26$p.value, digits=3), ")")
  
  # Map factor x-coordinates: 2 (1), 18 (2), 23 (3), 26 (4)
  # Color Palette: Deep Slate Navy (#1F3A60) for Fasting, Coral Salmon (#E07A5F) for After-Meal
  p <- ggplot(plot_data, aes(x = timepoint, y = mean_val, group = Group, color = Group)) +
    # Gridlines and Theme
    theme_minimal(base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#4A4A4A", linewidth = 0.6),
      axis.ticks = element_line(color = "#4A4A4A", linewidth = 0.6),
      text = element_text(color = "#2D3748"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b=10)),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    # Error Bars (mean ± SEM)
    geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val), 
                  width = 0.15, linewidth = 0.8, alpha = 0.8, show.legend = FALSE) +
    # Line Geoms
    geom_line(linewidth = 1.2, alpha = 0.9) +
    # Dotted recovery line connecting NS4 After-Meal (23) to RTDS Fasting (26)
    geom_line(data = recovery_line, aes(x = timepoint, y = mean_val, group = 1), 
              linetype = "dotted", color = "#E07A5F", linewidth = 1.2, alpha = 0.9) +
    # Points Geoms
    geom_point(size = 3.5, stroke = 1.0, fill = "white", shape = 21, show.legend = TRUE) +
    # Scales
    scale_x_discrete(labels = c("2" = "Baseline\n(2)", "18" = "NS4 Fasting\n(18)", "23" = "NS4 After-Meal\n(23)", "26" = "RTDS Fasting\n(26)")) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c("Fasting" = "#1F3A60", "After-Meal" = "#E07A5F")) +
    labs(
      title = paste("Plasma", get_ba_full_name(comp), "Trajectory"),
      x = "Experimental Phase (Timepoint ID)",
      y = "Concentration (Scaled Peak Area, Log10 Scale)"
    ) +
    # Bracket 1: 23 (x=3) vs 26 (x=4)
    annotate("segment", x = 3, xend = 4, y = bracket_y1, yend = bracket_y1, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = bracket_y1, yend = bracket_y1 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 4, xend = 4, y = bracket_y1, yend = bracket_y1 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("text", x = 3.5, y = bracket_y1 * 1.1, label = p_26_star, size = 3.2, fontface = "bold", color = "#2D3748") +
    
    # Bracket 2: 23 (x=3) vs 18 (x=2)
    annotate("segment", x = 2, xend = 3, y = bracket_y2, yend = bracket_y2, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y2, yend = bracket_y2 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = bracket_y2, yend = bracket_y2 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("text", x = 2.5, y = bracket_y2 * 1.1, label = p_18_star, size = 3.2, fontface = "bold", color = "#2D3748") +
    
    # Bracket 3: 23 (x=3) vs 2 (x=1)
    annotate("segment", x = 1, xend = 3, y = bracket_y3, yend = bracket_y3, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y3, yend = bracket_y3 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = bracket_y3, yend = bracket_y3 * 0.9, color = "#4A4A4A", linewidth = 0.4) +
    annotate("text", x = 2.0, y = bracket_y3 * 1.1, label = p_2_star, size = 3.2, fontface = "bold", color = "#2D3748") +
    
    coord_cartesian(ylim = c(min_mean_sem * 0.8, y_limit_max))
  
  # Save the plot
  plot_filename <- paste0("plots/", comp, "_trajectory.png")
  ggsave(plot_filename, plot = p, width = 6.5, height = 5.2, dpi = 300)
  
  # Create simplified plot for the comprehensive grid
  p_grid <- ggplot(plot_data, aes(x = timepoint, y = mean_val, group = Group, color = Group)) +
    theme_minimal(base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#4A4A4A", linewidth = 0.4),
      axis.ticks = element_line(color = "#4A4A4A", linewidth = 0.4),
      text = element_text(color = "#2D3748"),
      plot.title = element_text(face = "bold", size = 7, hjust = 0.5, margin = margin(b=2)),
      axis.title = element_blank(),
      axis.text = element_text(size = 5, face = "bold"),
      legend.position = "bottom"
    ) +
    geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val), 
                  width = 0.1, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    # Dotted recovery line connecting NS4 After-Meal (23) to RTDS Fasting (26)
    geom_line(data = recovery_line, aes(x = timepoint, y = mean_val, group = 1), 
              linetype = "dotted", color = "#E07A5F", linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 1.5, stroke = 0.5, fill = "white", shape = 21, show.legend = TRUE) +
    scale_x_discrete(labels = c("2" = "2", "18" = "18", "23" = "23", "26" = "26")) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c("Fasting" = "#1F3A60", "After-Meal" = "#E07A5F")) +
    labs(title = comp)
  
  grid_plots[[comp]] <- p_grid
}

cat("\nStatistical modeling and figure generation completed successfully!\n")

# ---------------------------------------------------------
# 3.1 Apply FDR (Benjamini-Hochberg) Correction
# ---------------------------------------------------------
cat("\nApplying FDR (Benjamini-Hochberg) correction across all pairwise contrasts...\n")

# Collect all raw p-values for individual bile acids (15 compounds x 3 contrasts = 45 tests)
ba_p_values <- c()
for (comp in ba_cols) {
  s <- stats_results[[comp]]
  ba_p_values <- c(ba_p_values, s$p_23_vs_2, s$p_23_vs_18, s$p_23_vs_26)
}
ba_q_values <- p.adjust(ba_p_values, method = "BH")

# Store FDR q-values back into stats_results for individual BAs
idx <- 1
for (comp in ba_cols) {
  stats_results[[comp]]$q_23_vs_2 <- ba_q_values[idx]
  stats_results[[comp]]$q_23_vs_18 <- ba_q_values[idx + 1]
  stats_results[[comp]]$q_23_vs_26 <- ba_q_values[idx + 2]
  idx <- idx + 3
}

# Separately correct aggregate classes (4 classes x 3 contrasts = 12 tests)
agg_p_values <- c()
for (comp in comp_cols) {
  s <- stats_results[[comp]]
  agg_p_values <- c(agg_p_values, s$p_23_vs_2, s$p_23_vs_18, s$p_23_vs_26)
}
agg_q_values <- p.adjust(agg_p_values, method = "BH")

idx <- 1
for (comp in comp_cols) {
  stats_results[[comp]]$q_23_vs_2 <- agg_q_values[idx]
  stats_results[[comp]]$q_23_vs_18 <- agg_q_values[idx + 1]
  stats_results[[comp]]$q_23_vs_26 <- agg_q_values[idx + 2]
  idx <- idx + 3
}

cat("FDR correction applied: ", length(ba_p_values), " tests for individual BAs, ", length(agg_p_values), " tests for aggregate classes.\n")

# ---------------------------------------------------------
# 3.5 Generate Comprehensive Grid of All 15 Plots
# ---------------------------------------------------------
cat("\nGenerating comprehensive 15-plot grid...\n")

# Reorder grid plots to match the biological class presentation order
grid_plots_ordered <- grid_plots[ba_order]

# Wrap plots using patchwork
composite_plot <- wrap_plots(grid_plots_ordered, ncol = 3, nrow = 5) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold")
  )

# Save grid image
ggsave("plots/all_trajectories_grid.png", plot = composite_plot, width = 8.27, height = 11.69, dpi = 300)
cat("Grid plot saved to plots/all_trajectories_grid.png\n")

# ---------------------------------------------------------
# 3.6 Generate Aggregate Classes Grid (2x2 Layout)
# ---------------------------------------------------------
cat("\nGenerating aggregate classes 2x2 grid...\n")
agg_plots_ordered <- grid_plots[comp_cols]
# Adjust font and size slightly for the 2x2 grid to make it look premium
for (nm in comp_cols) {
  agg_plots_ordered[[nm]] <- agg_plots_ordered[[nm]] +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 9, face = "bold")
    )
}

composite_agg_plot <- wrap_plots(agg_plots_ordered, ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold")
  )

ggsave("plots/composite_classes_grid.png", plot = composite_agg_plot, width = 7.5, height = 6.5, dpi = 300)
cat("Aggregate grid plot saved to plots/composite_classes_grid.png\n")

# ---------------------------------------------------------
# 3.8 Generate Master Summary Table
# ---------------------------------------------------------
cat("\nGenerating Master Statistical Summary Table...\n")

ba_order_all <- c("GCA", "TCA", "GCDCA", "TCDCA", "CA", "CDCA", "GDCA", "TDCA", "GLCA", "TLCA", "DCA", "LCA", "UDCA", "GUDCA", "TUDCA",
                  "Total_BA", "Glycine_Conjugated", "Taurine_Conjugated", "Unconjugated")

master_rows <- list()
for (comp in ba_order_all) {
  s <- stats_results[[comp]]
  
  master_rows[[comp]] <- data.frame(
    Compound = comp,
    Class = get_ba_category_short(comp),
    LFC_2 = round(s$est_23_vs_2, 3),
    P_2 = format_p_value(s$p_23_vs_2),
    Q_2 = paste0(format_p_value(s$q_23_vs_2), " (", get_stars(s$q_23_vs_2), ")"),
    LFC_18 = round(s$est_23_vs_18, 3),
    P_18 = format_p_value(s$p_23_vs_18),
    Q_18 = paste0(format_p_value(s$q_23_vs_18), " (", get_stars(s$q_23_vs_18), ")"),
    LFC_26 = round(s$est_23_vs_26, 3),
    P_26 = format_p_value(s$p_23_vs_26),
    Q_26 = paste0(format_p_value(s$q_23_vs_26), " (", get_stars(s$q_23_vs_26), ")"),
    stringsAsFactors = FALSE
  )
}
master_df <- do.call(rbind, master_rows)

master_ft <- flextable(master_df) %>%
  theme_vanilla() %>%
  set_header_labels(
    Compound = "Bile Acid",
    Class = "Class / Category",
    LFC_2 = "LFC",
    P_2 = "p (raw)",
    Q_2 = "FDR q",
    LFC_18 = "LFC",
    P_18 = "p (raw)",
    Q_18 = "FDR q",
    LFC_26 = "LFC",
    P_26 = "p (raw)",
    Q_26 = "FDR q"
  ) %>%
  add_header_row(
    values = c("", "", 
               "Postprandial (23) vs. Baseline (2)", 
               "Postprandial (23) vs. NS4 Fasting (18)", 
               "Postprandial (23) vs. RTDS Fasting (26)"),
    colwidths = c(1, 1, 3, 3, 3)
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Compound", "Class"), align = "left", part = "all") %>%
  fontsize(size = 8, part = "all") %>%
  fontsize(size = 7.5, part = "body") %>%
  padding(padding = 2, part = "all")

for (i in 1:nrow(master_df)) {
  comp <- master_df$Compound[i]
  s <- stats_results[[comp]]
  
  if (!is.na(s$q_23_vs_2) && s$q_23_vs_2 < 0.05) {
    master_ft <- color(master_ft, i = i, j = c("LFC_2", "P_2", "Q_2"), color = "#C0392B", part = "body")
    master_ft <- bold(master_ft, i = i, j = c("LFC_2", "P_2", "Q_2"), part = "body")
  }
  if (!is.na(s$q_23_vs_18) && s$q_23_vs_18 < 0.05) {
    master_ft <- color(master_ft, i = i, j = c("LFC_18", "P_18", "Q_18"), color = "#C0392B", part = "body")
    master_ft <- bold(master_ft, i = i, j = c("LFC_18", "P_18", "Q_18"), part = "body")
  }
  if (!is.na(s$q_23_vs_26) && s$q_23_vs_26 < 0.05) {
    master_ft <- color(master_ft, i = i, j = c("LFC_26", "P_26", "Q_26"), color = "#C0392B", part = "body")
    master_ft <- bold(master_ft, i = i, j = c("LFC_26", "P_26", "Q_26"), part = "body")
  }
}
master_ft <- autofit(master_ft)

# ---------------------------------------------------------
# 4. Generate Microsoft Word Document (.docx)
# ---------------------------------------------------------
cat("\nStarting Word document generation...\n")

doc <- read_docx()

# Style formatting helper
set_doc_styles <- function(doc) {
  # We just use officer API directly
  doc
}

# Header / Title Page
doc <- doc %>%
  body_add_par("ENTEROHEPATIC BILE ACID DYNAMICS REPORT", style = "heading 1") %>%
  body_add_par("Longitudinal Crossover Analysis of Fasting vs. Postprandial States across 15 Bile Acids", style = "Normal") %>%
  body_add_par("Statistical Method: Batch-corrected Linear Mixed-effects Model (LMM)", style = "Normal") %>%
  body_add_par(paste("Date Generated:", Sys.Date()), style = "Normal") %>%
  body_add_break()

# Executive Summary & Methodology Section
doc <- doc %>%
  body_add_par("1. Executive Summary & Research Context", style = "heading 1") %>%
  body_add_par("Bile acids (BAs) are critical amphipathic sterol molecules synthesized from cholesterol in hepatocytes, acting as essential detergents for emulsification and absorption of dietary lipids in the small intestine. Beyond their physical detergent properties, bile acids function as potent signaling ligands that bind to nuclear receptors, such as the farnesoid X receptor (FXR), and G-protein coupled receptors, like TGR5. These actions regulate lipid, glucose, energy, and xenobiotic metabolism throughout the enterohepatic circulation pathway. This report presents a rigorous longitudinal analysis of 15 primary, secondary, and tertiary bile acids to characterize baseline levels, prolonged fasting trajectories, and postprandial (after meal) physiological response.") %>%
  
  body_add_par("2. Detailed Statistical & Preprocessing Methodology", style = "heading 1") %>%
  body_add_par("To ensure chemical and biological validity, a multi-stage statistical pipeline was developed and executed:") %>%
  
  body_add_par("2.1 Technical Plate-to-Plate Batch Effect Normalization", style = "heading 2") %>%
  body_add_par("Because the biological subjects are fully nested within individual plates (e.g., all 5 timepoints for a single subject are processed on one plate), simple biological median centering would confound real biological subject variation with plate-level technical bias. To resolve this, technical batch correction was achieved using identical pooled Long-Term Reference (LTR) samples run repeatedly across all three plates (SRVp008, SRVp009, SRVp010). For each compound, a plate-specific adjustment factor was calculated as the ratio of the global LTR median to the plate-specific LTR median. This technical adjustment was applied to all sample concentrations to normalize technical plate bias while fully preserving individual biological baselines and crossover profiles.") %>%
  
  body_add_par("2.2 Longitudinal Linear Mixed-effects Modeling (LMM)", style = "heading 2") %>%
  body_add_par("A major statistical error in longitudinal or crossover data analysis is treating repeated measurements on the same subject as independent observations (using standard ANOVA or simple linear models). This violates the core independence assumption, inflates type-I error rates, and miscalculates p-values. Here, we implemented a Linear Mixed-effects Model (LMM) in R using the lmerTest package. The model is formulated as:") %>%
  body_add_par("  Value_log ~ Timepoint + (1 | Subject)", style = "Normal") %>%
  body_add_par("where Timepoint is a categorical fixed effect (levels: 2 = Baseline, 18 = NS4 Fasting, 26 = RTDS Fasting, 23 = NS4 After-Meal) and (1 | Subject) is the random effect modeling a subject-specific intercept. The random intercept mathematically models the unique baseline of each individual and successfully accounts for intra-subject covariance due to repeated measures.") %>%
  
  body_add_par("2.3 Log-Transformation and Model Validity", style = "heading 2") %>%
  body_add_par("Plasma concentrations of bile acids are strictly non-negative and highly skewed, naturally spanning multiple orders of magnitude. Analyzing raw concentrations directly in an LMM would severely violate the core assumptions of residual normality and homoscedasticity. Therefore, a natural log-transformation was applied to all adjusted values prior to fitting the mixed models, ensuring variance stabilization and residual normality.") %>%
  
  body_add_par("2.4 Post-hoc Pairwise Comparisons & FDR Correction", style = "heading 2") %>%
  body_add_par("Post-hoc pairwise comparisons were conducted using the emmeans package in R. Pairwise contrasts were specifically targeted to compare the postprandial state (Timepoint 23) against each of the three fasting states: Baseline (Timepoint 2), NS4 Fasting (Timepoint 18), and RTDS Fasting (Timepoint 26). P-values were calculated using Satterthwaite's approximation for degrees of freedom. To control the false discovery rate (FDR) arising from multiple testing across all 15 bile acids (45 pairwise tests), Benjamini-Hochberg (BH) correction was applied. Statistical significance is determined based on FDR-corrected q-values (q < 0.05) and annotated using standard stars (*q < 0.05, **q < 0.01, ***q < 0.001, ns = non-significant). Individual trajectory plots display unadjusted p-values for reference; all summary tables and interpretive text use FDR-corrected q-values.") %>%
  body_add_break()

# ---------------------------------------------------------
# 3. Comprehensive Statistical Summary Table (Master Table)
# ---------------------------------------------------------
doc <- doc %>%
  body_add_par("3. Comprehensive Statistical Summary Table (Master Table)", style = "heading 1") %>%
  body_add_par("The following table provides a complete, publication-ready overview of the statistical comparisons for all 15 individual bile acids and 4 aggregate classes (Total Bile Acids, Glycine-Conjugated, Taurine-Conjugated, and Unconjugated pools). Log Fold Changes (LFC) are reported on the natural log scale as estimated by the Linear Mixed-effects Models (LMM). Both raw p-values (Satterthwaite) and FDR-corrected q-values (Benjamini-Hochberg) are reported. Statistical significance is determined by FDR q-values and highlighted in bold dark red for all contrasts yielding q < 0.05 (*q < 0.05, **q < 0.01, ***q < 0.001, ns = non-significant).") %>%
  body_add_par("") %>%
  body_add_flextable(master_ft) %>%
  body_add_break()

# Detailed Results per Compound Section
doc <- doc %>%
  body_add_par("4. Detailed Bile Acid Trajectory Analysis & Interpretations", style = "heading 1")

# Sort bile acids by class for logical presentation (already defined above)

for (comp in ba_order) {
  cat("Adding", comp, "to Word document...\n")
  
  full_name <- get_ba_full_name(comp)
  category <- get_ba_category(comp)
  s <- stats_results[[comp]]
  
  # Create a summary table for LMM contrasts
  raw_p <- c(s$p_23_vs_2, s$p_23_vs_18, s$p_23_vs_26)
  raw_q <- c(s$q_23_vs_2, s$q_23_vs_18, s$q_23_vs_26)
  
  tbl_df <- data.frame(
    Contrast = c("Postprandial (23) vs. Baseline (2)", 
                 "Postprandial (23) vs. NS4 Fasting (18)", 
                 "Postprandial (23) vs. RTDS Fasting (26)"),
    `Log_Fold_Change` = round(c(s$est_23_vs_2, s$est_23_vs_18, s$est_23_vs_26), 4),
    `Std_Error` = round(c(s$se_23_vs_2, s$se_23_vs_18, s$se_23_vs_26), 4),
    `P_Value` = format.pval(raw_p, digits = 4),
    `FDR_q` = format.pval(raw_q, digits = 4),
    Significance = c(get_stars(s$q_23_vs_2), get_stars(s$q_23_vs_18), get_stars(s$q_23_vs_26)),
    q_num = raw_q,
    stringsAsFactors = FALSE
  )
  
  ft <- flextable(tbl_df[, c("Contrast", "Log_Fold_Change", "Std_Error", "P_Value", "FDR_q", "Significance")]) %>%
    theme_vanilla() %>%
    autofit() %>%
    set_header_labels(
      Contrast = "Contrast Comparison",
      Log_Fold_Change = "Log Fold Change (LFC)",
      Std_Error = "Std Error (SE)",
      P_Value = "p-value (raw)",
      FDR_q = "FDR q-value",
      Significance = "Signif"
    ) %>%
    color(i = which(tbl_df$q_num < 0.05), color = "#C0392B") %>%
    bold(i = which(tbl_df$q_num < 0.05), part = "body") %>%
    align(align = "center", part = "all")
  
  # Dynamic biological interpretation text based on actual statistical results
  bio_desc <- generate_bio_desc(comp, s)
  
  # Append to Word Document
  plot_filename <- paste0("plots/", comp, "_trajectory.png")
  
  doc <- doc %>%
    body_add_par(paste("4.", match(comp, ba_order), " - ", full_name), style = "heading 2") %>%
    body_add_par(paste("Bile Acid Classification:", category), style = "Normal") %>%
    body_add_par("") %>%
    body_add_img(src = plot_filename, width = 5.2, height = 4.16, style = "centered") %>%
    body_add_par("Figure Caption: Mean concentration trajectory (Scaled Peak Area, log10 scale) showing the Fasting trajectory (Baseline [2] -> NS4 Fasting [18] -> RTDS Fasting [26], blue line) and the Meal trajectory (overlapping fasting, and extending to NS4 After-Meal [23], orange line). Error bars represent mean ± SEM. Statistical brackets show unadjusted p-values from longitudinal crossover Linear Mixed-effects Models (LMM). FDR-corrected q-values are reported in the statistical table and interpretive text. *p < 0.05, **p < 0.01, ***p < 0.001, ns = non-significant.", style = "Normal") %>%
    body_add_par("") %>%
    body_add_par("Statistical Contrast Results Table:", style = "Normal") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par("Statistical Relevance & Biological Interpretation:", style = "Normal") %>%
    body_add_par(bio_desc, style = "Normal") %>%
    body_add_break()
}

# ---------------------------------------------------------
# 4. Generate Aggregate Classes Section
# ---------------------------------------------------------
cat("\nAdding Aggregate Classes Section to Word document...\n")

doc <- doc %>%
  body_add_par("5. Aggregate Bile Acid Class Trajectory Analysis & Interpretations", style = "heading 1") %>%
  body_add_par("This section provides a high-level physiological analysis of the four major bile acid pool classes: Total Bile Acids, Glycine-Conjugated Bile Acids, Taurine-Conjugated Bile Acids, and Unconjugated Bile Acids. Evaluating these broad chemical classes offers a comprehensive overview of enterohepatic circulation capacity, conjugation preferences, active transport efficacy, and the metabolic influence of gut microbial biotransformations.") %>%
  body_add_par("") %>%
  body_add_par("5.1 Aggregate Classes 2x2 Trajectory Grid", style = "heading 2") %>%
  body_add_img(src = "plots/composite_classes_grid.png", width = 6.0, height = 5.2, style = "centered") %>%
  body_add_par("Figure Caption: 2x2 composite trajectory grid of the four aggregate bile acid classes showing the Fasting trajectory (Baseline [2] -> NS4 Fasting [18] -> RTDS Fasting [26], blue line) and the Meal trajectory (starts at Baseline [2] and connects directly to NS4 After-Meal [23], orange line, with a dotted line connecting to RTDS Fasting [26]). Error bars represent mean ± SEM. Values are in Scaled Peak Area, log10 scale.", style = "Normal") %>%
  body_add_break()

for (comp in comp_cols) {
  cat("Adding aggregate:", comp, "to Word document...\n")
  
  full_name <- get_ba_full_name(comp)
  category <- get_ba_category(comp)
  s <- stats_results[[comp]]
  
  # Create a summary table for LMM contrasts
  raw_p <- c(s$p_23_vs_2, s$p_23_vs_18, s$p_23_vs_26)
  raw_q <- c(s$q_23_vs_2, s$q_23_vs_18, s$q_23_vs_26)
  
  tbl_df <- data.frame(
    Contrast = c("Postprandial (23) vs. Baseline (2)", 
                 "Postprandial (23) vs. NS4 Fasting (18)", 
                 "Postprandial (23) vs. RTDS Fasting (26)"),
    `Log_Fold_Change` = round(c(s$est_23_vs_2, s$est_23_vs_18, s$est_23_vs_26), 4),
    `Std_Error` = round(c(s$se_23_vs_2, s$se_23_vs_18, s$se_23_vs_26), 4),
    `P_Value` = format.pval(raw_p, digits = 4),
    `FDR_q` = format.pval(raw_q, digits = 4),
    Significance = c(get_stars(s$q_23_vs_2), get_stars(s$q_23_vs_18), get_stars(s$q_23_vs_26)),
    q_num = raw_q,
    stringsAsFactors = FALSE
  )
  
  ft <- flextable(tbl_df[, c("Contrast", "Log_Fold_Change", "Std_Error", "P_Value", "FDR_q", "Significance")]) %>%
    theme_vanilla() %>%
    autofit() %>%
    set_header_labels(
      Contrast = "Contrast Comparison",
      Log_Fold_Change = "Log Fold Change (LFC)",
      Std_Error = "Std Error (SE)",
      P_Value = "p-value (raw)",
      FDR_q = "FDR q-value",
      Significance = "Signif"
    ) %>%
    color(i = which(tbl_df$q_num < 0.05), color = "#C0392B") %>%
    bold(i = which(tbl_df$q_num < 0.05), part = "body") %>%
    align(align = "center", part = "all")
  
  # Dynamic biological interpretation text based on actual statistical results
  bio_desc <- generate_bio_desc(comp, s)
  
  # Append to Word Document
  plot_filename <- paste0("plots/", comp, "_trajectory.png")
  
  doc <- doc %>%
    body_add_par(paste("5.", 1 + match(comp, comp_cols), " - ", full_name), style = "heading 2") %>%
    body_add_par(paste("Bile Acid Classification:", category), style = "Normal") %>%
    body_add_par("") %>%
    body_add_img(src = plot_filename, width = 5.2, height = 4.16, style = "centered") %>%
    body_add_par("Figure Caption: Mean class concentration trajectory (Scaled Peak Area, log10 scale) showing the Fasting trajectory (Baseline [2] -> NS4 Fasting [18] -> RTDS Fasting [26], blue line) and the Meal trajectory (starts at Baseline [2] and connects directly to NS4 After-Meal [23], orange line, with a dotted recovery line connecting to RTDS Fasting [26]). Error bars represent mean ± SEM. Statistical brackets show unadjusted p-values from longitudinal crossover Linear Mixed-effects Models (LMM). FDR-corrected q-values are reported in the statistical table and interpretive text. *p < 0.05, **p < 0.01, ***p < 0.001, ns = non-significant.", style = "Normal") %>%
    body_add_par("") %>%
    body_add_par("Statistical Contrast Results Table:", style = "Normal") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par("Statistical Relevance & Biological Interpretation:", style = "Normal") %>%
    body_add_par(bio_desc, style = "Normal") %>%
    body_add_break()
}

# ---------------------------------------------------------
# 5. Append Comprehensive Grid Page to Report
# ---------------------------------------------------------
cat("\nAppending comprehensive grid page to report...\n")

doc <- doc %>%
  body_add_par("6. Comprehensive Bile Acid Trajectory Grid (All 15 Compounds)", style = "heading 1") %>%
  body_add_par("The following grid combines the longitudinal trajectories of all 15 bile acids in a single comprehensive page. This overview facilitates comparative analysis of baseline levels, prolonged fasting trajectories, and postprandial responses across primary, secondary, and tertiary bile acid families. The X-axis timepoint IDs represent Baseline (2), NS4 Fasting (18), NS4 After-Meal (23), and RTDS Fasting (26).") %>%
  body_add_par("") %>%
  body_add_img(src = "plots/all_trajectories_grid.png", width = 6.5, height = 9.19, style = "centered")

# Save Word document
print(doc, target = "Bile_Acid_Analysis_Report.docx")
cat("\nWord document Bile_Acid_Analysis_Report.docx generated successfully!\n")

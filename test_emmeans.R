library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

df <- read.csv("IRMA23930_MS-BA_Results.csv")
samples_raw <- df[df$sampleType == "sample", ]
tube_split <- strsplit(as.character(samples_raw$tubeLabel), "-")
samples_raw$subject <- sapply(tube_split, function(x) x[1])
samples_raw$timepoint <- sapply(tube_split, function(x) x[2])

# Filter for the relevant timepoints
samples_sub <- samples_raw[samples_raw$timepoint %in% c("2", "18", "23", "26"), ]
samples_sub$timepoint <- factor(samples_sub$timepoint, levels = c("2", "18", "26", "23"))

samples_sub$y_log <- log(samples_sub$CDCA)

model <- lmer(y_log ~ timepoint + (1 | subject), data = samples_sub)
em <- emmeans(model, ~ timepoint)
pairs_all <- as.data.frame(pairs(em, adjust = "none"))

cat("Contrast Names in emmeans output:\n")
print(pairs_all$contrast)
print(pairs_all)

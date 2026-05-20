df <- read.csv("IRMA23930_MS-BA_Results.csv")
samples <- df[df$sampleType == "sample", ]
ba_cols <- c("CDCA", "GCDCA", "TCDCA", "UDCA", "GUDCA", "TUDCA", 
             "DCA", "GDCA", "LCA", "GLCA", "TLCA", "CA", "GCA", "TCA", "TDCA")

for (comp in ba_cols) {
  vals <- samples[[comp]]
  num_na <- sum(is.na(vals))
  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  num_zero <- sum(vals == 0, na.rm = TRUE)
  num_neg <- sum(vals < 0, na.rm = TRUE)
  cat(sprintf("%-10s | Min: %10.4f | Max: %12.4f | NA: %3d | Zeros: %2d | Neg: %2d\n", 
              comp, min_val, max_val, num_na, num_zero, num_neg))
}

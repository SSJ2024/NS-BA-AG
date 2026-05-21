df <- read.csv("IRMA23930_MS-BA_Results.csv")
cat("Dimensions: ", dim(df), "\n")
cat("Sample Type table:\n")
print(table(df$sampleType, useNA = "always"))

samples <- df[df$sampleType == "sample", ]
cat("\nFirst 10 tube labels:\n")
print(head(samples$tubeLabel, 10))

# Check how subject and timepoints split from tube labels
tube_split <- strsplit(as.character(samples$tubeLabel), "-")
subjects <- sapply(tube_split, function(x) x[1])
timepoints <- sapply(tube_split, function(x) x[2])
cat("\nUnique Subjects:\n")
print(unique(subjects))
cat("\nUnique Timepoints:\n")
print(table(timepoints, useNA = "always"))

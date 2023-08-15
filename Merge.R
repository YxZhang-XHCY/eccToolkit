# Load necessary library
library(dplyr)

# Get command line arguments
args <- commandArgs()

# Load data
sample1 <- read.csv(file = paste(args[6], "merge.csv", sep = "."), header = FALSE, sep = "\t")

# Select the relevant columns
sample2 <- select(sample1, V1, V2, V3, V12, V13, V14, V5, V6)

# Initialize new columns
sample2$df1 <- numeric(nrow(sample2))
sample2$df2 <- numeric(nrow(sample2))
sample2$class <- character(nrow(sample2))

# Calculate and classify differences
for (i in 1:nrow(sample2)) {
  sample2$df1[i] <- sample2$V13[i] - sample2$V2[i]
  sample2$df2[i] <- sample2$V14[i] - sample2$V3[i]

  sample2$class[i] <- ifelse(abs(sample2$df1[i]) <= 10 & abs(sample2$df2[i]) <= 10, "HiConf", "LoConf")
}

# Select the required columns and add sample and tissue info
sample3 <- select(sample2, V1, V2, V3, class, V5, V6)
sample3$sample <- args[6]
sample3$tissue <- args[7]
sample3$length <- sample3$V3 - sample3$V2 + 1
sample3$df <- paste(sample3$V1, paste(sample3$V2, sample3$V3, sep = "-"), sep = ": ")

# Reorder and rename columns for clarity
sample4 <- select(sample3, tissue, sample, class, df, V1, V2, V3, length, V5, V6)
names(sample4) <- c("Tissue", "Sample", "Class", "Name", "Chr", "Start", "End", "Length", "Split reads", "Circle score")

# Write the result to a CSV file
write.csv(sample4, file = paste(args[6], "eccDNA.obtained.csv", sep = "."))

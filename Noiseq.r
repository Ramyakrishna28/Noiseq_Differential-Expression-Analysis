# =======================================================
# Differential Expression Analysis with NOISeq
# =======================================================

# Description:
#   - Reads transcript quantification files from control and patient samples
#   - Summarizes counts at the gene/transcript level
#   - Performs differential expression analysis using NOISeq
#   - Outputs DEG tables (full + filtered)
# =======================================================

# ---- 0. Load libraries ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("NOISeq", character.only = TRUE)) {
  BiocManager::install("NOISeq")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(NOISeq)
})

# ---- 1. Define file paths (generalized) ----
# Place your files in the "data" directory
control_file <- "data/control_quant.csv"
patient_file <- "data/patient_quant.csv"

# ---- 2. Load input data ----
control_df <- read.csv(control_file, stringsAsFactors = FALSE)
patient_df <- read.csv(patient_file, stringsAsFactors = FALSE)

# Ensure column naming consistency
colnames(control_df)[colnames(control_df) == "Name"] <- "Gene_Symbol"
colnames(patient_df)[colnames(patient_df) == "Name"] <- "Gene_Symbol"

# ---- 3. Summarize counts at the gene/transcript level ----
control_df_unique <- control_df %>%
  group_by(Gene_Symbol) %>%
  summarise(NumReads = sum(NumReads, na.rm = TRUE), .groups = "drop")

patient_df_unique <- patient_df %>%
  group_by(Gene_Symbol) %>%
  summarise(NumReads = sum(NumReads, na.rm = TRUE), .groups = "drop")

# ---- 4. Merge control & patient counts ----
merged_df <- inner_join(control_df_unique, patient_df_unique,
                        by = "Gene_Symbol",
                        suffix = c("_Control", "_Patient"))

count_matrix <- data.frame(
  row.names = merged_df$Gene_Symbol,
  Control = merged_df$NumReads_Control,
  Patient = merged_df$NumReads_Patient
)

# ---- 5. Experimental design ----
factors <- data.frame(
  Sample = c("Control", "Patient"),
  Condition = c("Control", "Patient")
)

# ---- 6. Prepare data for NOISeq ----
mydata <- readData(data = count_matrix, factors = factors)

# ---- 7. Run NOISeq ----
mynoiseq <- noiseq(
  mydata,
  factor = "Condition",
  norm = "tmm",
  replicates = "no",  # only one control & one patient
  pnr = 0.2,
  nss = 5,
  v = 0.02,
  lc = 1
)

# Extract results table
mynoiseq_results <- mynoiseq@results[[1]]

# ---- 8. Create results directory ----
if (!dir.exists("results")) {
  dir.create("results")
}

# Save unfiltered results
write.csv(mynoiseq_results, "results/NOISeq_DEG_results.csv", row.names = TRUE)

# ---- 9. Filter DEGs ----
deg_filtered <- mynoiseq_results[
  abs(mynoiseq_results$M) > 1 &    # log2 fold-change > 1
    mynoiseq_results$prob > 0.95,  # probability > 0.95
]

# Add transcript/gene IDs as column
deg_filtered$Transcript_ID <- rownames(deg_filtered)

# Reorder columns
deg_filtered <- deg_filtered[, c("Transcript_ID", setdiff(names(deg_filtered), "Transcript_ID"))]

# Save filtered DEGs
write.csv(deg_filtered, "results/filtered_DEG_m-1_prob-0.95.csv", row.names = FALSE)

cat("Analysis complete.\n",
    "Full results: results/NOISeq_DEG_results.csv\n",
    "Filtered DEGs: results/filtered_DEG_m-1_prob-0.95.csv\n")

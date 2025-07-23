# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(gprofiler2)
BiocManager::install("gprofiler2")
# 1. Load your Excel file (adjust the file path as needed)
df <- read_excel("/Users/brandon.oswald/Desktop/brandon_lab_data/Experiments/BMO-170/proteomics_results.xlsx")

# 2. Identify replicate columns by condition (same as before)
all_cols <- names(df)
al_wt_cols  <- grep("AL-WT", all_cols, value = TRUE)
al_dko_cols <- grep("AL-DKO", all_cols, value = TRUE)
dr_wt_cols  <- grep("DR-WT", all_cols, value = TRUE)
dr_dko_cols <- grep("DR-DKO", all_cols, value = TRUE)

# Combine all replicate columns into a single vector for easy processing
replicate_cols <- c(al_wt_cols, al_dko_cols, dr_wt_cols, dr_dko_cols)
replicate_cols

# 3. Force these replicate columns to be numeric.
#    "NaN" strings will convert to numeric NaN, and any non-numeric strings become NA.
df[replicate_cols] <- lapply(df[replicate_cols], function(x) {
  as.numeric(x)  # R will parse "NaN" as NaN, otherwise parse failures become NA
})

# Now you can safely run rowMeans() without the "must be numeric" error
df <- df %>%
  mutate(
    AL_WT_mean  = rowMeans(select(., all_of(al_wt_cols)), na.rm = TRUE),
    AL_DKO_mean = rowMeans(select(., all_of(al_dko_cols)), na.rm = TRUE),
    DR_WT_mean  = rowMeans(select(., all_of(dr_wt_cols)), na.rm = TRUE),
    DR_DKO_mean = rowMeans(select(., all_of(dr_dko_cols)), na.rm = TRUE)
  )
df

# Continue with your fold-change, log2-FC, t-tests, and plots as before.
df <- df %>%
  mutate(
    AL_fold_change = AL_DKO_mean / AL_WT_mean,
    DR_fold_change = DR_DKO_mean / DR_WT_mean,
    AL_log2_fc = log2(AL_fold_change),
    DR_log2_fc = log2(DR_fold_change)
  )

# Example: Simple t-test (same code as before)
compute_p_value <- function(row_values, group1_cols, group2_cols) {
  values1 <- as.numeric(row_values[group1_cols])
  values2 <- as.numeric(row_values[group2_cols])
  if (sum(!is.na(values1)) > 1 && sum(!is.na(values2)) > 1) {
    return(t.test(values1, values2, var.equal = FALSE)$p.value)
  } else {
    return(NA)
  }
}

df$AL_pval <- apply(df, 1, function(row) compute_p_value(row, al_wt_cols, al_dko_cols))
df$DR_pval <- apply(df, 1, function(row) compute_p_value(row, dr_wt_cols, dr_dko_cols))

# ------------------------------
# 0. Setup & Required Libraries
# ------------------------------
# (Install packages if needed)
# BiocManager::install(c("gprofiler2", "EnhancedVolcano", "fgsea", "msigdbr"))

library(readxl)
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)

# ------------------------------
# 1. Load and Process Data
# ------------------------------
# Load your Excel file (adjust path as needed)
df <- read_excel("/Users/brandon.oswald/Desktop/brandon_lab_data/Experiments/BMO-170/proteomics_results.xlsx")

# Identify replicate columns by condition
all_cols <- names(df)
al_wt_cols  <- grep("AL-WT",  all_cols, value = TRUE)
al_dko_cols <- grep("AL-DKO", all_cols, value = TRUE)
dr_wt_cols  <- grep("DR-WT",  all_cols, value = TRUE)
dr_dko_cols <- grep("DR-DKO", all_cols, value = TRUE)

# Combine replicate columns for conversion
replicate_cols <- c(al_wt_cols, al_dko_cols, dr_wt_cols, dr_dko_cols)

# Force replicate columns to be numeric (converting "NaN" strings appropriately)
df[replicate_cols] <- lapply(df[replicate_cols], function(x) as.numeric(x))

# Compute per-condition means
df <- df %>%
  mutate(
    AL_WT_mean  = rowMeans(select(., all_of(al_wt_cols)), na.rm = TRUE),
    AL_DKO_mean = rowMeans(select(., all_of(al_dko_cols)), na.rm = TRUE),
    DR_WT_mean  = rowMeans(select(., all_of(dr_wt_cols)), na.rm = TRUE),
    DR_DKO_mean = rowMeans(select(., all_of(dr_dko_cols)), na.rm = TRUE)
  )

# Compute fold changes & log2 fold changes for our two comparisons:
# (a) AL WT vs DR WT and (b) AL WT vs DR DKO.
df <- df %>%
  mutate(
    WT_fold_change         = DR_WT_mean / AL_WT_mean,
    WT_log2_fc             = log2(WT_fold_change),
    ALWT_DRDKO_fold_change = DR_DKO_mean / AL_WT_mean,
    ALWT_DRDKO_log2_fc     = log2(ALWT_DRDKO_fold_change)
  )

df

# Define a function to compute a simple t-test between two sets of replicates
compute_p_value <- function(row_values, group1_cols, group2_cols) {
  values1 <- as.numeric(row_values[group1_cols])
  values2 <- as.numeric(row_values[group2_cols])
  if (sum(!is.na(values1)) > 1 && sum(!is.na(values2)) > 1) {
    return(t.test(values1, values2, var.equal = FALSE)$p.value)
  } else {
    return(NA)
  }
}

# Compute p-values for:
# (a) AL WT vs DR WT and (b) AL WT vs DR DKO.
df$WT_pval <- apply(df, 1, function(row) compute_p_value(row, al_wt_cols, dr_wt_cols))
df$ALWT_DRDKO_pval <- apply(df, 1, function(row) compute_p_value(row, al_wt_cols, dr_dko_cols))

# ------------------------------
# 2. Volcano Plots (EnhancedVolcano)
# ------------------------------

# (a) AL WT vs DR WT Volcano Plot
res_WT <- data.frame(
  log2FC = df$WT_log2_fc,
  pvalue = df$WT_pval,
  Genes  = df$`PG.Genes`
)
res_WT
pdf("Volcano_ALWT_vs_DRWT.pdf", width = 8, height = 6)
EnhancedVolcano(
  res_WT,
  lab = res_WT$Genes,
  x = 'log2FC',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = 'DR-WT vs AL-WT',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'p-value'),
  pointSize = 2.0,
  labSize = 3.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  ylim = c(0, 5),
  xlim = c(-4, 5)
)
dev.off()
library(dplyr)
up_proteins <- res_WT %>% filter(log2FC <= -1, pvalue < 0.05)
up_proteins
EnhancedVolcano(
  res_WT,
  lab = res_WT$Genes,
  x = 'log2FC',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = 'DR-WT vs AL-WT',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'p-value'),
  pointSize = 2.0,
  labSize = 3.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  drawConnectors = TRUE,           # Enables connector lines
  widthConnectors = 0.5,           # Adjust connector line width as needed
  colConnectors = "grey50",        # Color for the connector lines
  ylim = c(0, 10),
)

# Install writexl if you haven't already:
install.packages("writexl")
library(writexl)

# Export res_WT to an Excel file named "res_WT.xlsx"
write_xlsx(res_WT, "res_WT.xlsx")



# (b) AL WT vs DR DKO Volcano Plot
res_DRDKO <- data.frame(
  log2FC = df$ALWT_DRDKO_log2_fc,
  pvalue = df$ALWT_DRDKO_pval,
  Genes  = df$`PG.Genes`
)

pdf("Volcano_ALWT_vs_DRDKO.pdf", width = 3, height = 2)
EnhancedVolcano(
  res_DRDKO,
  lab = res_DRDKO$Genes,
  x = 'log2FC',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = 'DR-DKO vs AL-WT',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'p-value'),
  pointSize = 2.0,
  labSize = 3.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2')
)
dev.off()

# ------------------------------
# 3. GSEA Analysis on AL WT vs DR WT Using C7 Gene Sets
# ------------------------------

# Retrieve C7 gene sets for mouse using msigdbr and filter out type of results
# Get the C2 (curated) gene sets for mouse and filter for those with "metab" in the gene set name.
# Retrieve MSigDB gene sets for mouse
msig_mouse <- msigdbr(species = "Mus musculus")
# Retrieve gene sets from C2, Hallmark (H), and C5 separately
msig_C2 <- msigdbr(species = "Mus musculus", category = "C2")
msig_H  <- msigdbr(species = "Mus musculus", category = "H")
msig_C5 <- msigdbr(species = "Mus musculus", category = "C5")


# Load required libraries
library(readxl)
library(dplyr)

# 1. Read the Excel file.
# Adjust "path_to_your_excel_file.xlsx" to your file's path.
gene_data <- read_excel("/Users/brandon.oswald/Desktop/brandon_lab_data/Experiments/BMO-170/gene_list_danials.xlsx", col_names = TRUE)

# 2. Make sure the columns are correctly named. For example, if your file has headers:
# "Group" and "Gene". If not, rename them accordingly:
colnames(gene_data) <- c("Group", "Gene")

# 3. Ensure the Gene column is of character type.
gene_data <- gene_data %>%
  mutate(Gene = as.character(Gene))

# 4. Split the gene names by group.
# This creates a named list where each element is a character vector of gene names.
gene_sets <- split(gene_data$Gene, gene_data$Group)

# 5. (Optional) To check the structure:
str(gene_sets)
print(gene_sets)

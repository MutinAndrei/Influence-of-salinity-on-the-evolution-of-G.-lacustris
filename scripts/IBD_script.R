# ==========================================
# FULL ANALYSIS: IBD/IBE with 8 genes
# and multiple testing correction (FDR)
# ==========================================

# Load packages
library(ape)
library(stringr)
library(geosphere)
library(vegan)
library(ecodist)
library(dplyr)
library(openxlsx)

# ==========================================
# 1. PREPARE DATA TABLE
# ==========================================

info <- read.xlsx("salt_table.xlsx", 4)

colnames(info)[1] <- "salinity_code"
colnames(info)[3] <- "MOTU"
colnames(info)[4] <- "Latitude"
colnames(info)[5] <- "Longitude"

clean_coord <- function(x) {
  x <- gsub("[^0-9.,-]", "", x)
  x <- gsub(",", ".", x)
  as.numeric(x)
}
info$Latitude <- clean_coord(info$Latitude)
info$Longitude <- clean_coord(info$Longitude)

fix_coords <- data.frame(
  MOTU = c("MOTU21", "MOTU75", "MOTU77", "MOTU85", "MOTU86", "MOTU87"),
  lon_correct = c(9.6469, NA, NA, 3.9865, 3.9865, 12.4832),
  lat_correct = c(NA, 26.02, 26.4, NA, NA, NA)
)
for (i in seq_len(nrow(fix_coords))) {
  idx <- which(info$MOTU == fix_coords$MOTU[i])
  if (length(idx) > 0) {
    if (!is.na(fix_coords$lon_correct[i])) info$Longitude[idx] <- fix_coords$lon_correct[i]
    if (!is.na(fix_coords$lat_correct[i])) info$Latitude[idx] <- fix_coords$lat_correct[i]
  }
}

info$salinity_num <- case_when(
  tolower(trimws(info$salinity_code)) == "fresh" ~ 0,
  tolower(trimws(info$salinity_code)) == "brackish" ~ 1,
  tolower(trimws(info$salinity_code)) == "salty" ~ 2,
  TRUE ~ NA_real_
)

# Convert MOTU to "MOTU###" format
info$MOTU <- ifelse(grepl("^MOTU", info$MOTU), info$MOTU, paste0("MOTU", info$MOTU))

# Remove rows with NA and invalid coordinates
info_clean <- info %>%
  filter(!is.na(MOTU), !is.na(Latitude), !is.na(Longitude), !is.na(salinity_num)) %>%
  filter(Longitude >= -180, Longitude <= 180, Latitude >= -90, Latitude <= 90)

# Remove duplicates
info_clean <- info_clean[!duplicated(info_clean$MOTU), ]
rownames(info_clean) <- info_clean$MOTU

cat("After table cleaning, remaining MOTUs:", nrow(info_clean), "\n")

# ==========================================
# 2. READ ALL 8 GENES AND EXTRACT MOTU
# ==========================================

gene_files <- c(
  "AK_lacustris_HOU_table_salt.fa.aln.trim",
  "COIII_lacustris_HOU_table_salt.fa.aln.trim",
  "COI_lacustris_HOU_table_salt.fa.aln.trim",
  "EF1a_lacustris_HOU_table_salt.fa.aln.trim",
  "HSP70_lacustris_HOU_table_salt.fa.aln.trim",
  "PEPCK_lacustris_HOU_table_salt.fa.aln.trim",
  "S12_16_lacustris_HOU_table_salt.fa.aln.trim",
  "S28_lacustris_HOU_table_salt.fa.aln.trim"
)

# Function to read a gene: extract MOTU, select sequence with maximum coverage
read_gene <- function(file) {
  dna <- read.dna(file, format = "fasta")
  motu <- str_extract(rownames(dna), "MOTU\\d+")
  keep <- !is.na(motu)
  dna <- dna[keep, ]
  motu <- motu[keep]
  
  coverage <- apply(as.character(dna), 1, function(x) mean(x != "-"))
  
  df_temp <- data.frame(motu = motu, cov = coverage, index = 1:length(motu))
  df_temp <- df_temp %>%
    group_by(motu) %>%
    filter(cov == max(cov)) %>%
    slice(1) %>%
    ungroup()
  
  dna <- dna[df_temp$index, ]
  motu <- df_temp$motu
  rownames(dna) <- motu
  return(dna)
}

gene_alignments <- list()
for (file in gene_files) {
  cat("Reading", file, "...\n")
  gene_alignments[[file]] <- read_gene(file)
}

motus_in_genes <- unique(unlist(lapply(gene_alignments, rownames)))
cat("Unique MOTUs in genes:", length(motus_in_genes), "\n")

# ==========================================
# 3. COMBINE FOR MOTUs PRESENT IN TABLE
# ==========================================

common_motus <- intersect(rownames(info_clean), motus_in_genes)
cat("Common MOTUs for analysis (before coverage filtering):", length(common_motus), "\n")
if (length(common_motus) == 0) stop("No common MOTUs — check MOTU format")

info_sel <- info_clean[common_motus, ]

gene_lengths <- sapply(gene_alignments, ncol)
total_length <- sum(gene_lengths)

concat_seqs <- matrix("-", nrow = length(common_motus), ncol = total_length,
                      dimnames = list(common_motus, NULL))

col_start <- 1
for (i in seq_along(gene_alignments)) {
  gene <- gene_alignments[[i]]
  gene_len <- ncol(gene)
  col_end <- col_start + gene_len - 1
  
  for (motu in common_motus) {
    if (motu %in% rownames(gene)) {
      concat_seqs[motu, col_start:col_end] <- as.character(gene)[motu, ]
    }
  }
  col_start <- col_end + 1
}

concat_dna <- as.DNAbin(concat_seqs)
cat("Concatenated matrix size:", dim(concat_dna), "\n")

# ==========================================
# 4. FILTER BY COVERAGE (GAP FRACTION <= 50%)
# ==========================================

gap_fraction <- apply(as.character(concat_dna), 1, function(x) mean(x == "-"))
keep_motu <- rownames(concat_dna)[gap_fraction <= 0.5]
cat("MOTUs after coverage filtering (<=50% gaps):", length(keep_motu), "\n")

concat_dna <- concat_dna[keep_motu, ]
info_sel <- info_sel[keep_motu, ]

# ==========================================
# 5. COMPUTE DISTANCE MATRICES
# ==========================================

gen_dist <- dist.dna(concat_dna, model = "raw")
geo_dist <- as.dist(geosphere::distm(cbind(info_sel$Longitude, info_sel$Latitude),
                                     fun = geosphere::distHaversine) / 1000)
sal_diff <- dist(info_sel$salinity_num, method = "manhattan")

# ==========================================
# 6. REMOVE MOTUs CAUSING NaN
# ==========================================

while (any(is.nan(gen_dist))) {
  gen_mat <- as.matrix(gen_dist)
  na_rows <- apply(gen_mat, 1, function(x) any(is.nan(x)))
  bad_motus <- rownames(gen_mat)[na_rows]
  cat("Removing MOTUs with NaN:", length(bad_motus), "\n")
  
  good_motus <- setdiff(rownames(gen_mat), bad_motus)
  if (length(good_motus) < 10) stop("Too few MOTUs after NaN filtering")
  
  concat_dna <- concat_dna[good_motus, ]
  info_sel <- info_sel[good_motus, ]
  
  gen_dist <- dist.dna(concat_dna, model = "raw")
  geo_dist <- as.dist(geosphere::distm(cbind(info_sel$Longitude, info_sel$Latitude),
                                       fun = geosphere::distHaversine) / 1000)
  sal_diff <- dist(info_sel$salinity_num, method = "manhattan")
}

cat("Final number of MOTUs:", nrow(info_sel), "\n")

# ==========================================
# 7. MANTEL AND MRM TESTS
# ==========================================

set.seed(123)

mantel_geo <- vegan::mantel(gen_dist, geo_dist, permutations = 999)
mantel_sal <- vegan::mantel(gen_dist, sal_diff, permutations = 999)
mantel_partial_sal <- vegan::mantel.partial(gen_dist, sal_diff, geo_dist, permutations = 999)
mantel_partial_geo <- vegan::mantel.partial(gen_dist, geo_dist, sal_diff, permutations = 999)
mrm_result <- ecodist::MRM(gen_dist ~ geo_dist + sal_diff, nperm = 999)

# ==========================================
# 8. EXTRACT P-VALUES AND APPLY FDR CORRECTION
# ==========================================

p_raw <- c(
  "Mantel: gen~geo" = mantel_geo$signif,
  "Mantel: gen~sal" = mantel_sal$signif,
  "Partial Mantel: sal|geo" = mantel_partial_sal$signif,
  "Partial Mantel: geo|sal" = mantel_partial_geo$signif,
  "MRM: geo" = mrm_result$coef["geo_dist", "pval"],
  "MRM: sal" = mrm_result$coef["sal_diff", "pval"]
)

# Apply FDR correction (Benjamini-Hochberg)
p_fdr <- p.adjust(p_raw, method = "BH")

# Summary table
results_df <- data.frame(
  Test = names(p_raw),
  p_raw = round(p_raw, 4),
  p_FDR = round(p_fdr, 4)
)

# Print test results
cat("\n=== TEST RESULTS ===\n")
print(mantel_geo)
print(mantel_sal)
print(mantel_partial_sal)
print(mantel_partial_geo)
print(mrm_result)

# Print FDR-adjusted p-values
cat("\n=== FDR CORRECTION ===\n")
print(results_df)

cat("\nSignificance after FDR (p < 0.05):\n")
print(results_df$p_FDR < 0.05)




# Load ggplot2 for plotting
library(ggplot2)

# Data: tests and their FDR p-values
results_fdr <- data.frame(
  Test = c("Mantel: genetic ~ geographic", 
           "Mantel: genetic ~ salinity", 
           "Partial Mantel: salinity | geography", 
           "Partial Mantel: geography | salinity", 
           "MRM: geography", 
           "MRM: salinity"),
  p_FDR = c(0.540, 0.066, 0.066, 0.540, 0.906, 0.066)
)

# Reorder for readability (top to bottom)
results_fdr$Test <- factor(results_fdr$Test, levels = rev(results_fdr$Test))

# Plot
ggplot(results_fdr, aes(x = Test, y = p_FDR)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "FDR-adjusted p-value",
    title = "Multiple testing corrected p-values (FDR)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  coord_flip()


ggsave("fdr_pvalues_plot.png", width = 10, height = 5, dpi = 300)
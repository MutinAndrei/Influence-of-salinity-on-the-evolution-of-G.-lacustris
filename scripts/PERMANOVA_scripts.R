# ============================================================
# PERMANOVA: salinity vs phylogenetic distances
# ============================================================

# Load packages
library(ape)
library(vegan)

# ------------------------------------------------------------
# 1. Read tree and prepare data
# ------------------------------------------------------------
tree <- read.tree("example.nex.treefile")
tree$tip.label <- gsub("brakish|brekish", "brackish", tree$tip.label)

# Ensure positive branch lengths
tree$edge.length[tree$edge.length <= 0] <- .Machine$double.eps
tree$edge.length[tree$edge.length < 1e-8] <- 1e-8

# Extract salinity from tip labels
extract_salinity <- function(x) {
  if (grepl("fresh", x, ignore.case = TRUE)) return("fresh")
  if (grepl("brackish|brakish|brekish", x, ignore.case = TRUE)) return("brackish")
  if (grepl("salty", x, ignore.case = TRUE)) return("salty")
  if (grepl("Lake_14", x, ignore.case = TRUE)) return("fresh")
  if (grepl("NuhuNur", x, ignore.case = TRUE)) return("fresh")
  if (grepl("NamishNur", x, ignore.case = TRUE)) return("brackish")
  return(NA)
}

tip_salinity <- sapply(tree$tip.label, extract_salinity)
names(tip_salinity) <- tree$tip.label

# Remove tips with unknown salinity
keep <- !is.na(tip_salinity)
tree_pruned <- drop.tip(tree, tree$tip.label[!keep])
tip_salinity <- tip_salinity[keep]
names(tip_salinity) <- tree_pruned$tip.label

cat("Salinity distribution:\n")
print(table(tip_salinity))

# ------------------------------------------------------------
# 2. Phylogenetic distance matrix
# ------------------------------------------------------------
phylo_dist <- cophenetic(tree_pruned)

# ------------------------------------------------------------
# 3. Salinity factor
# ------------------------------------------------------------
salinity_factor <- factor(tip_salinity, levels = c("fresh", "brackish", "salty"))
names(salinity_factor) <- tree_pruned$tip.label

# ------------------------------------------------------------
# 4. PERMANOVA
# ------------------------------------------------------------
set.seed(123)
permanova_result <- adonis2(phylo_dist ~ salinity_factor, permutations = 999)

cat("\n===== PERMANOVA =====\n")
print(permanova_result)

# ------------------------------------------------------------
# 5. (Optional) Pairwise PERMANOVA with Holm correction
# ------------------------------------------------------------
pairwise_permanova <- function(dist_mat, groups, p.adjust.m = "holm") {
  groups <- as.factor(groups)
  lev <- levels(groups)
  comb <- combn(lev, 2)
  results <- list()
  for (i in 1:ncol(comb)) {
    idx <- which(groups %in% comb[, i])
    sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
    sub_groups <- factor(groups[idx], levels = comb[, i])
    set.seed(123)
    ad <- adonis2(sub_dist ~ sub_groups, permutations = 999)
    results[[paste(comb[1, i], "vs", comb[2, i])]] <- ad
  }
  pvals <- sapply(results, function(x) x$`Pr(>F)`[1])
  padj <- p.adjust(pvals, method = p.adjust.m)
  data.frame(Comparison = names(pvals), p_raw = pvals, p_adjusted = padj)
}

pairwise_result <- pairwise_permanova(phylo_dist, salinity_factor, p.adjust.m = "holm")
cat("\n===== Pairwise PERMANOVA (Holm correction) =====\n")
print(pairwise_result)

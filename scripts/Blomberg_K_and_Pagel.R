# ============================================================
# Blomberg's K for phylogenetic tree
# ============================================================

# Load required packages
library(ape)
library(phytools)

# 1. Read tree
tree <- read.tree("example.nex.treefile")

# Ensure positive branch lengths (some methods require >0)
tree$edge.length[tree$edge.length <= 0] <- .Machine$double.eps
tree$edge.length[tree$edge.length < 1e-8] <- 1e-8

# 2. Extract salinity category from tip labels
extract_salinity <- function(x) {
  if (grepl("fresh", x, ignore.case = TRUE)) {
    return("fresh")
  } else if (grepl("brackish|brakish|brekish", x, ignore.case = TRUE)) {
    return("brackish")
  } else if (grepl("salty", x, ignore.case = TRUE)) {
    return("salty")
  } else if (grepl("Lake_14", x, ignore.case = TRUE)) {
    return("fresh")       # Lake_14 is freshwater
  } else if (grepl("NuhuNur", x, ignore.case = TRUE)) {
    return("fresh")       # NuhuNur is freshwater
  } else if (grepl("NamishNur", x, ignore.case = TRUE)) {
    return("brackish")    # NamishNur is brackish
  } else {
    return(NA)            # unknown, will be excluded
  }
}

tip_salinity <- sapply(tree$tip.label, extract_salinity)
names(tip_salinity) <- tree$tip.label

# 3. Remove tips with unknown salinity
keep <- !is.na(tip_salinity)
tree_pruned <- drop.tip(tree, tree$tip.label[!keep])
tip_salinity <- tip_salinity[keep]
names(tip_salinity) <- tree_pruned$tip.label

# Check distribution
cat("Salinity distribution:\n")
print(table(tip_salinity))

# 4. Convert salinity to numeric: fresh=0, brackish=1, salty=2
salinity_num <- as.numeric(factor(tip_salinity, levels = c("fresh","brackish","salty"))) - 1
names(salinity_num) <- tree_pruned$tip.label

# 5. Blomberg's K
set.seed(123)   # for reproducibility
K_test <- phylosig(tree_pruned, salinity_num, method = "K", test = TRUE, nsim = 100000)

# 6. Print result
cat("\n===== Blomberg's K =====\n")
print(K_test)




# ============================================================
# Pagel's lambda for phylogenetic tree
# ============================================================

# Load required packages
library(ape)
library(phytools)

# 1. Read tree
tree <- read.tree("example.nex.treefile")

# Ensure positive branch lengths
tree$edge.length[tree$edge.length <= 0] <- .Machine$double.eps
tree$edge.length[tree$edge.length < 1e-8] <- 1e-8

# 2. Extract salinity category from tip labels
extract_salinity <- function(x) {
  if (grepl("fresh", x, ignore.case = TRUE)) {
    return("fresh")
  } else if (grepl("brackish|brakish|brekish", x, ignore.case = TRUE)) {
    return("brackish")
  } else if (grepl("salty", x, ignore.case = TRUE)) {
    return("salty")
  } else if (grepl("Lake_14", x, ignore.case = TRUE)) {
    return("fresh")       # Lake_14 is freshwater
  } else if (grepl("NuhuNur", x, ignore.case = TRUE)) {
    return("fresh")       # NuhuNur is freshwater
  } else if (grepl("NamishNur", x, ignore.case = TRUE)) {
    return("brackish")    # NamishNur is brackish
  } else {
    return(NA)            # unknown, will be excluded
  }
}

tip_salinity <- sapply(tree$tip.label, extract_salinity)
names(tip_salinity) <- tree$tip.label

# 3. Remove tips with unknown salinity
keep <- !is.na(tip_salinity)
tree_pruned <- drop.tip(tree, tree$tip.label[!keep])
tip_salinity <- tip_salinity[keep]
names(tip_salinity) <- tree_pruned$tip.label

# Check distribution
cat("Salinity distribution:\n")
print(table(tip_salinity))

# 4. Convert salinity to numeric: fresh=0, brackish=1, salty=2
salinity_num <- as.numeric(factor(tip_salinity, levels = c("fresh","brackish","salty"))) - 1
names(salinity_num) <- tree_pruned$tip.label

# 5. Pagel's λ
lambda_test <- phylosig(tree_pruned, salinity_num, method = "lambda", test = TRUE)

# 6. Print result
cat("\n===== Pagel's lambda =====\n")
print(lambda_test)

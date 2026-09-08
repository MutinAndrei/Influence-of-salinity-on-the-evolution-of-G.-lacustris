# ============================================================
# Regions
# ============================================================

# Load packages
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(dplyr)

# Read and prepare tree
tree <- read.tree("example.nex.treefile")
tree$tip.label <- gsub("brakish|brekish", "brackish", tree$tip.label)

# Re-root by node 251
new_root_node <- 251
rerooted_tree_ape <- root(tree, node = new_root_node, resolve.root = TRUE)
new_root_node <- tree
# Define groups and colors
patterns <- c(
  "TSG-1" = "TSG-1",
  "TSG-2" = "TSG-2",
  "TSG-3" = "TSG-3",
  "TSG-4" = "TSG-4",
  "TSG-5" = "TSG-5",
  "EUG-1" = "EUG-1",
  "EUG-2" = "EUG-2",
  "HIG"   = "HIG",
  "NAG"   = "NAG",
  "EAG"   = "EAG",
  "NTG"   = "NTG",
  "STG"   = "STG",
  "Angara" = "Angara",
  "Maloe" = "Maloe",
  "Bormashovoe" = "Bormashovoe",
  "NamishNur" = "NamishNur",
  "NuhuNur" = "NuhuNur",
  "Lake_14" = "Lake_14",
  "Kirenga" = "Kirenga"
)

colors <- c(
  "TSG-1" = "#0521fb",
  "TSG-2" = "#0521fb",
  "TSG-3" = "#0521fb",
  "TSG-4" = "#0521fb",
  "TSG-5" = "#0521fb",
  "EUG-1" = "#fbad07",
  "EUG-2" = "#fbad07",
  "HIG"   = "#d006f3",
  "NAG"   = "#67dad8",
  "EAG"   = "#5b8c10",
  "NTG"   = "#f90312",
  "STG"   = "#7f1019",
  "Angara" = "#5b8c10",
  "Maloe" = "#5b8c10",
  "Bormashovoe" = "#5b8c10",
  "NamishNur" = "#5b8c10",
  "NuhuNur" = "#5b8c10",
  "Lake_14" = "#5b8c10",
  "Kirenga" = "#5b8c10"
)

# Assign groups to each tip
tip_groups <- sapply(rerooted_tree_ape$tip.label, function(label) {
  found <- names(patterns)[sapply(patterns, function(pat) grepl(pat, label, ignore.case = TRUE))]
  if (length(found) == 0) return(NA) else return(found[1])
})

tip_data <- data.frame(
  label = rerooted_tree_ape$tip.label,
  group = tip_groups,
  stringsAsFactors = FALSE
)

# Build base tree
p <- ggtree(rerooted_tree_ape, layout = "circular", branch.length = "none") %<+% tip_data +
  geom_tree(color = "gray45", size = 2) +
  geom_tippoint(aes(color = group, x = x + 1), size = 9, alpha = 1) +
  geom_nodepoint(
    aes(subset = !isTip & as.numeric(label) >= 95),
    color = "black",
    size = 5,
    shape = 15
  ) +
  scale_color_manual(
    values = colors,
    na.value = "gray45",
    guide = "none"
  ) +
  theme_tree() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(-1.2, -1.2, -1.2, -1.2), "cm")
  )

p
# Collapse clades
nodes_to_collapse <- c(301, 271, 249, 272, 245, 236, 276, 229, 279, 222, 185)
for (n in nodes_to_collapse) {
  p <- collapse(p, node = n)
}

# Add black triangles at collapsed nodes
collapsed_nodes_data <- p$data[p$data$node %in% nodes_to_collapse, ]
p_final <- p + 
  geom_point(data = collapsed_nodes_data, aes(x = x, y = y),
             shape = 17, size = 9, color = "black", inherit.aes = FALSE)
print(p_final)

# Print and save
ggsave("tree_final_clen_08092926.png", p, width = 14, height = 14, dpi = 300)
ggsave("tree_final_clean_08092926.svg", p_final, width = 14, height = 14, dpi = 300)
print(p_final)




# ============================================================
# Salinity
# ============================================================

# Load packages
library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(stringr)

# Read and prepare tree
tree <- read.tree("example.nex.treefile")
tree$tip.label <- gsub("brakish|brekish", "brackish", tree$tip.label)

# Re-root
new_root_node <- 251
rerooted_tree_ape <- root(tree, node = new_root_node, resolve.root = TRUE)

# Define salinity categories for tip coloring
salinity_cat <- function(label) {
  if (grepl("fresh", label, ignore.case = TRUE)) return("Fresh")
  if (grepl("brackish|brakish|brekish", label, ignore.case = TRUE)) return("Brackish")
  if (grepl("salty", label, ignore.case = TRUE)) return("Salty")
  if (grepl("Angara|Maloe|Kirenga|Lake_14|NuhuNur", label, ignore.case = TRUE)) return("Fresh")
  if (grepl("Bormashovoe", label, ignore.case = TRUE)) return("Salty")
  if (grepl("NamishNur", label, ignore.case = TRUE)) return("Brackish")
  return(NA)
}

tip_data <- data.frame(
  label = rerooted_tree_ape$tip.label,
  group = sapply(rerooted_tree_ape$tip.label, salinity_cat),
  stringsAsFactors = FALSE
)

color_map <- c("Fresh" = "#04859D", "Brackish" = "#FFA100", "Salty" = "#A101A6")

# Define substrings for sample highlighting
highlight_patterns <- c("NuhuNur", "Kirenga", "Bormashovoe", "NamishNur", 
                        "Lake_14", "Maloe", "Angara")

tip_data <- tip_data %>%
  mutate(is_sample = sapply(label, function(l) 
    any(sapply(highlight_patterns, function(pat) grepl(pat, l, ignore.case = TRUE)))))

# Build base tree
p <- ggtree(rerooted_tree_ape, layout = "circular", branch.length = "none") %<+% tip_data

# Draw all branches in gray
p <- p + geom_tree(color = "gray45", size = 2)

# Get tip nodes corresponding to samples
sample_tip_nodes <- p$data$node[p$data$isTip & p$data$is_sample]

# Function to get all ancestors of a node (including itself)
get_all_ancestors <- function(tree, node) {
  anc <- c()
  while (node != Ntip(tree) + 1) {
    anc <- c(anc, node)
    # Находим родителя
    edge <- which(tree$edge[,2] == node)
    if (length(edge) == 0) break
    node <- tree$edge[edge,1]
  }
  anc <- c(anc, node) 
  unique(anc)
}

# Collect all nodes along paths to samples
highlight_nodes <- unique(unlist(lapply(sample_tip_nodes, function(n) get_all_ancestors(rerooted_tree_ape, n))))

# Add colored branches for those nodes
p <- p + geom_tree(
  data = subset(p$data, node %in% highlight_nodes),
  color = "red3", size = 2
)

# Add tip points colored by salinity
p <- p + geom_tippoint(
  aes(color = group, x = x + 1.3),
  size = 7,
  alpha = 1,
  na.rm = TRUE
)

# Add bootstrap support points (>=95)
p <- p + geom_nodepoint(
  aes(subset = !isTip & as.numeric(label) >= 95),
  color = "black",
  size = 5,
  shape = 15
)

# Apply color scale and theme
p <- p + scale_color_manual(values = color_map, na.value = "gray50", guide = "none") +
  theme_tree() +
  theme(legend.position = "none", plot.margin = unit(c(-1.2, -1.2, -1.2, -1.2), "cm"))

# Save and print
ggsave("tree_full_path_highlight_08092926.png", p, width = 14, height = 14, dpi = 300)
ggsave("tree_full_path_highlight_08092926.svg", p, width = 14, height = 14, dpi = 300)
print(p)

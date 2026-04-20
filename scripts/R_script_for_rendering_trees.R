#Packages
library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(stringr)


#Uploading a tree file
tree <- read.tree("example.nex.treefile")
tree$tip.label <- gsub("brakish", "brackish", tree$tip.label)
tree$tip.label <- gsub("brekish", "brackish", tree$tip.label)

#Re-rooting by 219 node
new_root_node <- 248
rerooted_tree_ape <- root(
  tree,
  node = new_root_node,
  resolve.root = TRUE  # Fixing the root structure
)
rerooted_tree_ape <- tree
#Function for rendering according to a custom pattern
create_color_groups <- function(tree, patterns, colors) {
  groups <- setNames(rep(NA, Ntip(tree)), tree$tip.label)
  
  # Assign groups according to patterns (with priority of the first match)
  for(i in seq_along(patterns)) {
    groups[str_detect(names(groups), patterns[i]) & is.na(groups)] <- i
  }
  
  # Creating data.frame for ggtree
  data.frame(
    label = names(groups),
    group = factor(groups, levels = seq_along(patterns), labels = names(patterns)),
    color = colors[groups]
  )
}

#Pattern/legend
patterns <- c(
  Fresh = "fresh", 
  Brackish = "brackish", 
  Salty = "salty",
  Fresh = "Angara",
  Fresh = "Maloe",
  Salty = "Bormashovoe",
  Brackish = "NamishNur",
  Brackish = "NuhuNur"
)
#Pattern/Color
colors <- c(
  "#04859D",
  "#FFA100", 
  "#A101A6",
  "#04859D",
  "#04859D", 
  "#A101A6",
  "#FFA100",
  "#FFA100"
)
#Use the function for rendering according to a custom pattern
tip_data <- create_color_groups(rerooted_tree_ape, patterns, colors)

#A script for visualizing a tree by salinity
ggtree(rerooted_tree_ape, 
       layout = "circular", 
       branch.length = "none") %<+% 
  filter(tip_data, !is.na(group)) + 
  geom_tree(aes(color = "gray45"), size = 3) + 
  geom_tippoint(
    aes(
      color = group, 
      x = x + 1.3 
    ), 
    size = 20,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  xlim(-5, NA) +
  geom_tiplab(
    size = 0,
    color = "transparent",
    offset = 0
  ) +
  geom_nodepoint(
    aes(
      subset = !isTip & as.numeric(label) >= 95,
      fill = "UF Bootstrap ≥95"
    ),
    color = "black",
    size = 8,
    alpha = 1,
    shape = 15
  ) +
  # Шкала цветов (Salinity)
  scale_color_manual(
    name = "Salinity:",
    values = setNames(colors, names(patterns)),
    na.value = "gray45",
    guide = guide_legend(
      override.aes = list(
        size = 5, 
        shape = 16
      ),
      order = 1
    )
  ) +
  scale_shape_manual(
    name = "Node support:",
    values = c("UF Bootstrap ≥95" = 16),
    guide = guide_legend(
      override.aes = list(
        size = 7, 
        color = "black",
        alpha = 0.9
      ),
      order = 2
    )
  ) +
  theme_tree2() +
  theme(
    legend.position = c(0.13, 0.97),
    legend.justification = "top",
    legend.title = element_text(face = "bold", size = 35),
    legend.text = element_text(size = 30),
    legend.spacing.y = unit(1, "cm"),
    legend.box.spacing = unit(2, "cm"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 2, face = "bold"),
    plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "cm")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )


# Save .png  
ggsave("phylogenetic_tree_salt_table_lake14_17122025_1.png", width = 35, height = 35, dpi = 300)
ggsave("phylogenetic_tree_salt_table_lake14_08122025_tips.png", width = 22, height = 22, dpi = 300)
# Save .svg
ggsave("phylogenetic_tree_salt_table_lake14_17122025.svg", width = 24, height = 24, dpi = 300)



###############Biogeography

#packages
library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(stringr)
library(phytools)

# Uploading a tree file
tree <- read.tree("example.nex.treefile")
tree$tip.label <- gsub("brakish", "brackish", tree$tip.label)
tree$tip.label <- gsub("brekish", "brackish", tree$tip.label)
#Re-rooting by 211 node
new_root_node <- 248
rerooted_tree_ape <- root(
  tree,
  node = new_root_node,
  resolve.root = TRUE  # Fixing the root structure
)
rerooted_tree_ape <- tree
#The biogeography pattern
patterns <- c(
  TSG = "TSG-5", 
  EUG = "EUG-1", 
  TSG = "TSG-1",
  HIG = "HIG",
  TSG = "TSG-2",
  NAG = "NAG",
  TSG = "TSG-3",
  EAG = "EAG",
  NTG = "NTG",
  TSG = "TSG-4",
  STG = "STG",
  EUG = "EUG-2",
  EAG = "Angara",
  EAG = "Maloe",
  EAG = "Bormashovoe",
  EAG = "NamishNur",
  EAG = "NuhuNur",
  EAG = "Lake_14",
  EAG = "Kirenga"
)
#Color by biogeography
colors <- c(
  "#0521fbff",
  "#fbad07ff", 
  "#0521fbff", 
  "#d006f3ff", 
  "#0521fbff", 
  "#67dad8ff",
  "#0521fbff",
  "#5b8c10ff",
  "#f90312ff",
  "#0521fbff",
  "#7f1019ff",
  "#fbad07ff",
  "#5b8c10ff",
  "#5b8c10ff",
  "#5b8c10ff",
  "#5b8c10ff",
  "#5b8c10ff",
  "#5b8c10ff",
  "#5b8c10ff"
)
#Use the function for rendering according to a custom pattern
tip_data <- create_color_groups(rerooted_tree_ape, patterns, colors)

#A script for visualizing a tree by biogeography
p <- ggtree(rerooted_tree_ape, 
            layout = "circular", 
       branch.length = "none") %<+% 
  filter(tip_data, !is.na(group)) + 
  geom_tree(aes(color = "gray45"), size = 3) + 
  geom_tippoint(
    aes(
      color = group, 
      x = x + 1.2
    ), 
    size = 20,
    na.rm = TRUE
  ) +
  geom_tiplab(
    size = 0,
    color = "transparent",
    offset = 0
  ) +
  geom_nodepoint(
    aes(
      subset = !isTip & as.numeric(label) >= 95,
      fill = "UF Bootstrap ≥95"
    ),
    color = "black",
    size = 8,
    alpha = 1,
    shape = 15
  ) +
  # Шкала цветов (Regions)
  scale_color_manual(
    name = "Regions:",
    values = setNames(colors, names(patterns)),
    na.value = "gray45",
    guide = guide_legend(
      override.aes = list(
        size = 5, 
        shape = 16
      ),
      order = 1
    )
  ) +
  scale_shape_manual(
    name = "Node support:",
    values = c("UF Bootstrap ≥95" = 16),
    guide = guide_legend(
      override.aes = list(
        size = 7, 
        color = "black",
        alpha = 0.9
      ),
      order = 2
    )
  ) +
  theme_tree2() +
  theme(
    legend.position = c(0.13, 0.97),
    legend.justification = "top",
    legend.title = element_text(face = "bold", size = 35),
    legend.text = element_text(size = 30),
    legend.spacing.y = unit(1, "cm"),
    legend.box.spacing = unit(2, "cm"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 2, face = "bold"),
    plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "cm")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) 

p


nodes_to_collapse <- c(184, 277, 225, 274, 231, 246, 270)

p_collapsed <- Reduce(function(p, node) collapse(p, node), 
                      nodes_to_collapse, 
                      init = p)

p_final <- p_collapsed + 
  # Добавляем маркеры на места схлопнутых узлов
  geom_point2(
    aes(subset = (node %in% nodes_to_collapse)), 
    shape = 17,   
    size = 10,
    color = "black",
    stroke = 1
  ) + 
  theme(plot.margin = unit(c(1, -.4, 1, 1), "cm"))

p_final
# Save .png
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_1.png", width = 24, height = 24, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_colaps_cir_1.png", width = 24, height = 24, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_colaps_ver.png", width = 18, height = 32, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_08122025.png", width = 24, height = 24, dpi = 300)
# Save .svg
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_tips.svg", width = 24, height = 24, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_colaps_cir_1.svg", width = 24, height = 24, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_17122025_colaps_ver.svg", width = 18, height = 32, dpi = 300)


ggtree(rerooted_tree_ape, 
       layout = "circular", 
       branch.length = "none", 
       aes(color = group), 
       size = 4) %<+% 
  filter(tip_data, !is.na(group)) +
  geom_tiplab(
    size = 0,
    color = "transparent",
    offset = 0
  ) + 
  geom_nodepoint(
    aes(
      subset = !isTip & as.numeric(label) >= 95,
      shape = "UF Bootstrap ≥95"
    ),
    color = "black",
    size = 7,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Regions:",
    values = setNames(colors, names(patterns)),
    na.value = "gray45",
    guide = guide_legend(
      override.aes = list(size = 5, shape = 15),
      order = 1
    )
  ) +
  scale_shape_manual(
    name = "Node support:",
    values = c("UF Bootstrap ≥95" = 16),
    guide = guide_legend(
      override.aes = list(
        size = 7, 
        color = "black",
        alpha = 0.9
      ),
      order = 2
    )
  ) +
  theme_tree2() +
  theme(
    legend.position = c(0.11, 0.97),
    legend.justification = "top",
    legend.title = element_text(face = "bold", size = 35),
    legend.text = element_text(size = 30),
    legend.spacing.y = unit(1, "cm"),
    legend.box.spacing = unit(2, "cm"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 3, face = "bold"),
    plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "cm")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) 
# Save .png
ggsave("phylogenetic_tree_salt_table_color_pattern_14_08122025_tips.png", width = 24, height = 24, dpi = 300)
ggsave("phylogenetic_tree_salt_table_color_pattern_14_08122025.png", width = 24, height = 24, dpi = 300)
# Save .svg
ggsave("phylogenetic_tree_salt_table_color_pattern_14.svg", width = 24, height = 24, dpi = 300)

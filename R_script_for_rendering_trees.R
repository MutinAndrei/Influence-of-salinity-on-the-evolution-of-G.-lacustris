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

#Re-rooting by 211 node
new_root_node <- 211
rerooted_tree_ape <- root(
  tree,
  node = new_root_node,
  resolve.root = TRUE  # Fixing the root structure
)

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
  Brackish = "NamishNur",
  Brackish = "NuhuNur",
  Fresh = "Kirenga"
)
#Pattern/Color
colors <- c(
  "#04859D",
  "#FFA100", 
  "#A101A6", 
  "#FFA100", 
  "#FFA100", 
  "#04859D"
)
#Use the function for rendering according to a custom pattern
tip_data <- create_color_groups(rerooted_tree_ape, patterns, colors)

#A script for visualizing a tree by salinity
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
      shape = "UF Bootstrap ≥95"  # Separate aesthetics for the form
    ),
    color = "black",
    size = 7,
    alpha = 0.9
  ) +
  # Scale for colors (salinity)
  scale_color_manual(
    name = "Salinity:",
    values = setNames(colors, names(patterns)),
    na.value = "gray45",
    guide = guide_legend(
      override.aes = list(
        size = 5, 
        shape = 15  # Squares for the color legend
      ),
      order = 1  # Display order
    )
  ) +
  # Scale for shape (node support)
  scale_shape_manual(
    name = "Node support:",
    values = c("UF Bootstrap ≥95" = 16),
    guide = guide_legend(
      override.aes = list(
        size = 7, 
        color = "black",  # Black color for dots
        alpha = 0.9
      ),
      order = 2  # Display order
    )
  ) +
  theme_tree2() +
  theme(
    legend.position = c(0.13, 0.97),
    legend.justification = "top",
    legend.title = element_text(face = "bold", size = 35),
    legend.text = element_text(size = 30),
    legend.spacing.y = unit(1, "cm"),
    legend.box.spacing = unit(2, "cm"),  # The space between the legends
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 2, face = "bold"),
    plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "cm")
  ) +
  # Legend order management
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

# Save .png  
ggsave("phylogenetic_tree_salt_table.png", width = 22, height = 22, dpi = 300)
# Save .svg
ggsave("phylogenetic_tree_salt_table.svg", width = 22, height = 22, dpi = 300)




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
rerooted_tree_ape$tip.label <- gsub("brakish", "brackish", rerooted_tree_ape$tip.label)
rerooted_tree_ape$tip.label <- gsub("brekish", "brackish", rerooted_tree_ape$tip.label)
#Re-rooting by 211 node
new_root_node <- 211
rerooted_tree_ape <- root(
  tree,
  node = new_root_node,
  resolve.root = TRUE  # Fixing the root structure
)

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
  EUG = "EUG-2"
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
  "#fbad07ff"
)
#Use the function for rendering according to a custom pattern
tip_data <- create_color_groups(rerooted_tree_ape, patterns, colors)

#A script for visualizing a tree by biogeography
ggtree(rerooted_tree_ape, 
       layout = "circular", 
       branch.length = "none", 
       aes(color = group), 
       size = 4) %<+% 
  filter(tip_data, !is.na(group)) + 
  
  # Выделение клад цветом
  geom_hilight(
    node = 174,
    fill = "#FF6B6B",
    alpha = 0.3,
    extend = 1
  ) +
  geom_hilight(
    node = 238,
    fill = "#FF6B6B",
    alpha = 0.3,
    extend = 1
  ) +
  geom_hilight(
    node = 166,
    fill = "#FF6B6B",
    alpha = 0.3,
    extend = 1
  ) +
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
ggsave("phylogenetic_tree_salt_table_color_pattern.png", width = 24, height = 24, dpi = 300)
# Save .svg
ggsave("phylogenetic_tree_salt_table_color_pattern.svg", width = 24, height = 24, dpi = 300)

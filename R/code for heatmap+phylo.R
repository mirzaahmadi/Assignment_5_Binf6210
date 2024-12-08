# Load necessary libraries
library(ape)
library(ggtree)
library(ComplexHeatmap)
library(cowplot)
library(circlize)

# Step 1: Create a toy phylogenetic tree
newick_tree <- "(Panthera_uncia_-Snow_Leopard:0.08096672328,(Panthera_leo_-Lion:0.01258227973,Panthera_onca_-Jaguar:0.02922096954):0.01198357028,(Panthera_pardus_-Leopard:0.03103365556,Panthera_tigris_-Tiger:0.02999752708):0.01164645155);"
tree_ <- read.tree(text = newick_tree)

# Step 2: Create a toy similarity matrix
sim <- matrix(
  c(
    1.0000000, 0.3333333, 0.0625, 0, 0.1666667,
    0.3333333, 1.0000000, 0.1875, 0, 0.3333333,
    0.0625000, 0.1875000, 1.0000, 0, 0.2500000,
    0.0000000, 0.0000000, 0.0000, 1, 0.0000000,
    0.1666667, 0.3333333, 0.2500, 0, 1.0000000
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(c("Panthera uncia", "Panthera tigris", "Panthera pardus", "Panthera onca", "Panthera leo"),
                  c("Panthera uncia", "Panthera tigris", "Panthera pardus", "Panthera onca", "Panthera leo"))
)

# Step 3: Plot the tree using ggtree with tip labels
tree_plot <- ggtree(tree_, layout = "rectangular") +
  geom_tiplab(size = 4, fontface = "italic") +  # Add tip labels
  theme_tree()

# Step 4: Create the heatmap using ComplexHeatmap
heatmap <- Heatmap(sim, name = "Similarity",
                   col = colorRamp2(c(0.3, 1), c("white", "blue")),
                   show_row_names = TRUE,  # Show row names (species names)
                   show_column_names = TRUE)  # Show column names (species names)

# Step 5: Convert the ggtree plot to a grob (graphical object)
tree_grob <- ggplotGrob(tree_plot)

# Step 6: Draw the heatmap as a grid object
heatmap_grob <- draw(heatmap, newpage = FALSE)

# Step 7: Combine tree and heatmap using cowplot's plot_grid
combined_plot <- plot_grid(
  tree_grob, # Tree plot on the left
  heatmap_grob, # Heatmap on the right
  ncol = 2, # Arrange in two columns
  rel_widths = c(1, 1.5) # Adjust relative widths to make heatmap wider
)

# Step 8: Display the combined plot
print(combined_plot)

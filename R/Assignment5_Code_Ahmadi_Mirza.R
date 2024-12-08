library(ape)
library(Biostrings)
library(bold)
library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(ggtree)
library(grDevices)
library(gridExtra)
library(knitr)
library(maps)
library(msa)
library(rentrez)
library(skimr)
library(tidyverse)
library(tools)
library(vioplot)



# 1. DATA ACQUISITION ---- 
# Acquire geographic, species, and sequence data from BOLD and NCBI

#panthera_df <- as_tibble(bold_specimens(taxon='panthera'))
BOLD_data <- read_tsv("./data/panthera_df.tsv")


#The commented-out code below details how COI sequence data was obtained from NCBI 

# Set 'retmax' to 0 to retrieve total hit counts for all COI sequences for each species, allowing for subsequent fetching of all sequences
# leopard_COI_seqs <- entrez_search(db="nucleotide", term='"Panthera pardus"[Organism] AND COI[Gene]', retmax = 0)
# max_leopard_COI_hits <- leopard_COI_seqs$count #21
# 
# snow_COI_seqs <- entrez_search(db="nucleotide", term='"Panthera uncia"[Organism] AND COI[Gene]', retmax = 0)
# max_snow_COI_hits <- snow_COI_seqs$count #5
# 
# lion_COI_seqs <- entrez_search(db="nucleotide", term='"Panthera leo"[Organism] AND COI[Gene]', retmax = 0)
# max_lion_COI_hits <- lion_COI_seqs$count #12
# 
# jaguar_COI_seqs <- entrez_search(db="nucleotide", term='"Panthera onca"[Organism] AND COI[Gene]', retmax = 0)
# max_jaguar_COI_hits <- jaguar_COI_seqs$count #16
# 
# tiger_COI_seqs <- entrez_search(db="nucleotide", term='"Panthera tigris"[Organism] AND COI[Gene]', retmax = 0)
# max_leopard_COI_hits <- tiger_COI_seqs$count #67
# 
# 
# leopard_COI_seqs <- entrez_search(db="nuccore", term='"Panthera pardus"[Organism] AND COI[Gene]', retmax = max_leopard_COI_hits, use_history = TRUE)
# snow_COI_seqs <- entrez_search(db="nuccore", term='"Panthera uncia"[Organism] AND COI[Gene]', retmax = max_snow_COI_hits, use_history = TRUE)
# lion_COI_seqs <- entrez_search(db="nuccore", term='"Panthera leo"[Organism] AND COI[Gene]', retmax = max_lion_COI_hits, use_history = TRUE)
# jaguar_COI_seqs <- entrez_search(db="nuccore", term='"Panthera onca"[Organism] AND COI[Gene]', retmax = max_jaguar_COI_hits, use_history = TRUE)
# tiger_COI_seqs <- entrez_search(db="nuccore", term='"Panthera tigris"[Organism] AND COI[Gene]', retmax = max_leopard_COI_hits, use_history = TRUE)
# 
# # # This exploratory function provides a visual confirmation of number of sequences found in returned lists of sequences
# sequence_count <- function(...) {
#   items <- list(...)
#   for (i in items) {
#     cat("Total sequences found:", i$count, "\n")
#     cat("Sample sequence IDs:", head(i$ids, 3), "\n")
#     cat("\n")
#   }
# }
# sequence_count(leopard_COI_seqs, snow_COI_seqs, lion_COI_seqs, jaguar_COI_seqs, tiger_COI_seqs)
# 
# 
# # #Obtain sequence data in FASTA format
# batch_size <- 500 #Batch size is chosen to be 500 to efficiently retrieve data in manageable chunks without overloading the server
# # 
# # #Create function which Fetches sequence data by looping through batches, retrieving sequence info, and accumulation in FASTA format
# fetch_data <- function(entrez_search_results, batch_size) {
#   search_results_total_count <- entrez_search_results$count
#   complete_sequence_data <- c()
# 
#   for (start in seq(0, search_results_total_count, by = batch_size)) {
#   batch <- entrez_fetch(db = "nuccore",
#                         web_history = entrez_search_results$web_history,
#                         rettype = "fasta",
#                         retmax = batch_size,
#                         retstart = start)
#   complete_sequence_data <- c(complete_sequence_data, batch)
# }
# 
# return(complete_sequence_data)
# }
# 
# #Pass in panthera species search results to fetch_data function in order to output a fasta of all sequences to later be converted to biostirngs
# leopard_fasta_seqs <- fetch_data(leopard_COI_seqs, batch_size)
# write(leopard_fasta_seqs, "../data/leopard.fasta", sep = "\n")
# 
# snow_fasta_seqs <- fetch_data(snow_COI_seqs, batch_size)
# write(snow_fasta_seqs, "../data/snow.fasta", sep = "\n")
# 
# lion_fasta_seqs <- fetch_data(lion_COI_seqs, batch_size)
# write(lion_fasta_seqs, "../data/lion.fasta", sep = "\n")
# 
# jaguar_fasta_seqs <- fetch_data(jaguar_COI_seqs, batch_size)
# write(jaguar_fasta_seqs, "../data/jaguar.fasta", sep = "\n")
# 
# tiger_fasta_seqs <- fetch_data(tiger_COI_seqs, batch_size)
# write(tiger_fasta_seqs, "../data/tiger.fasta", sep = "\n")


# Loop through each FASTA file, create data frames with standardized naming conventions for efficient filtering, manipulation, and analysis of COI sequences 
for (file in list.files("./data", full.names=TRUE)) {
  extension <- file_ext(file)
  if(extension != "fasta") { #Since we are only dealing with FASTA files, all other files can be ignored
    next
  } else {
    #create stringset variables with correct naming conventions to then be converted to dataframes
    base_name <- file_path_sans_ext(basename(file))
    variable_name_stringset <- paste0(base_name, "_coi_stringset") 
    stringset_data <- readDNAStringSet(file)
    assign(variable_name_stringset, stringset_data)
    
    variable_name_df <- paste0("df_", base_name)
    df_data <- data.frame(COI_title = names(stringset_data), COI_sequence = paste(stringset_data))
    assign(variable_name_df, df_data)
    
    rm(list = variable_name_stringset, stringset_data, df_data) # Remove unneeded stringset variables from the environment for clarity
    
  }
}



# 2. DATA EXPLORATION, FILTERING, AND QUALITY CONTROL ----

# FILTER and EXPLORE sequence data, geographic data, and sequence data from BOLD and NCBI


# Create a data frame containing only species names and  corresponding countries, as these columns are relevant for downstream analysis
panthera_only_countries_df <- BOLD_data %>% 
  select(species_name, country) 

# The 'skim' function offers a comprehensive data frame summary (i.e. NA counts, completeness, unique values) - It will be used for exploratory data analysis frequently
skim(panthera_only_countries_df)

# Filter rows where animals were recorded in countries outside their wild range, as these are likely non-wild individuals. This loop identifies unique countries for each species, allowing subsequent exclusion of non-native locations in the data sets.
for (species in unique(panthera_only_countries_df$species_name)) {
    cat("Species:", species, "\n") 
    cat("Countries:", unique(panthera_only_countries_df$country[panthera_only_countries_df$species_name == species]), "\n\n")
}

# After researching species' ranges, update data sets to exclude the following countries that are outside of species' wild range: 
pardus_countries_to_exclude <- c("Germany", "Brazil")
tigris_countries_to_exclude <- c("South Africa", "Canada", "Brazil")
leo_countries_to_exclude <- c("Thailand", "Russia", "Canada")


#Creating Data frames for MAIN visualization #1

# update 'panthera_only_countries_df' to exclude non-native countries (and locations that are not countries at all + NA values)
panthera_only_countries_df <- panthera_only_countries_df %>% 
  filter(!(species_name == "Panthera pardus" & country %in% pardus_countries_to_exclude)) %>% 
  filter(!(species_name == "Panthera tigris" & country %in% tigris_countries_to_exclude)) %>% 
  filter(!(species_name == "Panthera leo" & country %in% leo_countries_to_exclude)) %>% 
  filter(country != "Unrecoverable", country != "Exception - Zoological Park") %>% 
  na.omit()

skim(panthera_only_countries_df)


# This density and subsequent world_map data frames will be be used downstream for creation of chloropleth map
map_density_df <- panthera_only_countries_df %>% 
  group_by(country) %>% 
  summarise(unique_species_count = n_distinct(species_name)) %>% 
  arrange(country) # Sort the dataframe by species in alphabetical order for convenient manual checking of dataframe

world_map <- map_data("world")
merged_map_density_df <- world_map %>% # This dataframe will ultimately be used to map density data 
  left_join(map_density_df, by = c("region" = "country")) %>% 
  select(-subregion)

# This latitude/longitude data frames will be used downstream for creation of point map
panthera_LatLon_df <- BOLD_data %>% 
  select(species_name, lat, lon, country) %>%
  filter(!(species_name == "Panthera pardus" & country %in% pardus_countries_to_exclude)) %>% 
  filter(!(species_name == "Panthera tigris" & country %in% tigris_countries_to_exclude)) %>% 
  filter(!(species_name == "Panthera leo" & country %in% leo_countries_to_exclude)) %>% 
  filter(country != "Unrecoverable", country != "Exception - Zoological Park") %>% 
  select(species_name, lat, lon, country) %>% 
  na.omit()
  
skim(panthera_LatLon_df)


world_map <- map_data("world")
merged_map_coordinates_df <- world_map %>% #This dataframe will ultimately be used to map coordinates
  left_join(panthera_LatLon_df, by = c("region" = "country")) %>% 
  select(-subregion)

# Creating Data frames for MAIN visualization #3

# For downstream similarity analysis, a similarity matrix is created from species and country data
# Create a list of unique countries for each species
species_countries <- panthera_only_countries_df %>%
  group_by(species_name) %>%
  summarise(countries = list(unique(country)))

# Initialize an empty matrix to store the similarities
species_names <- unique(panthera_only_countries_df$species_name)
similarity_matrix <- matrix(0, nrow = length(species_names), ncol = length(species_names))
rownames(similarity_matrix) <- species_names
colnames(similarity_matrix) <- species_names

# Calculate the similarity (shared countries) between each pair of species
for (i in 1:length(species_names)) {
  for (j in i:length(species_names)) {
    # Get the countries for species i and species j
    countries_i <- species_countries$countries[species_countries$species_name == species_names[i]][[1]]
    countries_j <- species_countries$countries[species_countries$species_name == species_names[j]][[1]]
    
    # Calculate the intersection (shared countries) between the two species
    shared_countries <- length(intersect(countries_i, countries_j))
    
    # Find the maximum number of countries for the species pair to normalize
    max_countries <- max(length(countries_i), length(countries_j))
    
    # Fill the matrix with the normalized shared country count as a decimal to be used in a heatmap
    similarity_matrix[i, j] <- shared_countries / max_countries
    similarity_matrix[j, i] <- shared_countries / max_countries  # The matrix is symmetric
  }
}

# Common names for the species - used for naming conventions
common_names <- c("Snow Leopard", "Tiger", "Leopard", "Jaguar", "Lion")

# Create new row and column names by combining scientific names and common names
new_names <- paste(rownames(similarity_matrix), " (", common_names, ")", sep = "")

rownames(similarity_matrix) <- new_names
colnames(similarity_matrix) <- new_names



# EXPLORATORY FIG 1: Create a figure to visually display the number of unique countries each big cat species inhabits within its native range
summary_country_df <- panthera_only_countries_df %>%
  group_by(species_name) %>%
  summarise(number_of_unique_countries = n_distinct(country))

ggplot(summary_country_df, aes(x = species_name, y = number_of_unique_countries, fill = species_name)) + 
  geom_bar(stat = "identity", fill = "#FF6F00", color = "black") +
  labs(
    title = expression(bold("Number of Unique Countries Inhabited by") ~ bold(italic("Panthera")) ~ bold("Species")),
    x = "Species",
    y = "Number of Unique Countries"
  ) +
  scale_x_discrete(labels = c(
    "Panthera leo" = expression(italic("P. leo") ~ "(Lion) 🦁"), 
    "Panthera onca" = expression(italic("P. onca") ~ "(Jaguar) 🐾"),
    "Panthera pardus" = expression(italic("P. pardus") ~ "(Leopard) 🐆"),
    "Panthera tigris" = expression(italic("P. tigris") ~ "(Tiger) 🐯"),
    "Panthera uncia" = expression(italic("P. uncia") ~ "(Snow Leopard) ❄️")
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),           
    axis.title.y = element_text(size = 25, face = "bold"),            
    axis.text.x = element_text(size = 21),                          
    axis.text.y = element_text(size = 21),
    legend.position = "none"
  )



# Filtering Sequence Data

# This function explores all data frame sequences, generating an info table, in order to explore the sequences before and after filtering
present_sequence_info <- function(list_of_dataframes) {
  # Create list of data frame names for generating sequence summaries and establishing naming conventions
  names <- c("Panthera leo (Lion)", "Panthera onca (Jaguar)", "Panthera pardus (Leopard)", "Panthera tigris (Tiger)", "Panthera Uncia (Snow Leopard)")
  
  # Initialize an empty summary table
  summary_table <- data.frame(Species = character(),
                              Num_Sequences = numeric(),
                              Min_Length = numeric(),
                              Max_Length = numeric(),
                              Mean_Length = numeric(),
                              Median_Length = numeric(),
                              stringsAsFactors = FALSE)
  
  # Loop through all data frames to fill in the summary table for the sequences
  for (i in seq_along(list_of_dataframes)) {
    df <- list_of_dataframes[[i]]
    sequence_name <- names[i]
    
    # Create DNAStringSet object from COI_sequence column so that dnastringset summary functions can be used
    sequences <- DNAStringSet(df$COI_sequence)
    
    # Create a DNAStringSet object from the COI_sequence column to enable the use of DNAStringSet summary functions
    sequence_lengths <- width(sequences)
    num_sequences <- length(sequences)
    min_length <- min(sequence_lengths)
    max_length <- max(sequence_lengths)
    mean_length <- mean(sequence_lengths)
    median_length <- median(sequence_lengths)
    
    # Append the results to the summary_table to visualize summary information in table format
    summary_table <- rbind(summary_table, data.frame(
      Species = sequence_name,
      Num_Sequences = num_sequences,
      Min_Length = min_length,
      Max_Length = max_length,
      Mean_Length = mean_length,
      Median_Length = median_length
    ))
  }
  
  print(kable(summary_table))
}


# Define constants for sequence filtering
missing_data <- 0.01 # A missing_data value of 0.01 was chosen because a 1% threshold is commonly used in sequence data processing
length_var <- 50 # A length_var of 50 was chosen as a reasonable range around the median to account for natural variation in sequence lengths

# Assign names to the dataframes in the following list so we can access each data frame by name directly in upcoming loop
named_dataframes <- list(df_lion = df_lion, df_jaguar = df_jaguar, df_leopard = df_leopard, 
                   df_tiger = df_tiger, df_snow = df_snow)

# Iterate over the list and apply sequence filtering to each one (created 5 new filtered data frames)
for (df_name in names(named_dataframes)) {
  df <- named_dataframes[[df_name]]  # Get the dataframe from the list
  
  # Perform filtering operations
  filtered_df <- df %>% 
    filter(str_length(COI_sequence) <= 1000) %>% #Filter sequences over 1000 sequences to ensure no genomes are kept in filtered data set
    mutate(COI_sequence = str_remove_all(COI_sequence, "^N+|N+$|-")) %>% #Remove all 
    filter(str_count(COI_sequence, "N") <= (missing_data * str_count(COI_sequence))) %>% # To ensure high quality data, remove any sequences with over 1% of ambiguous nucleotides
    filter(str_count(COI_sequence) >= median(str_count(COI_sequence)) - length_var & 
             str_count(COI_sequence) <= median(str_count(COI_sequence)) + length_var) %>% # Remove sequences that are not close to the median to homogenise data
    select(COI_title, COI_sequence)
  
  # Dynamically name the new dataframe and assign it the filtered data
  assign(paste0("filtered_", df_name), filtered_df)
  
  rm(df, filtered_df)
}


# Compare sequence information before and after filtering for all five panthera dataframes
dataframes <- list(df_lion, df_jaguar, df_leopard, df_tiger, df_snow)
present_sequence_info(dataframes)

filtered_dataframes <- list(filtered_df_lion, filtered_df_jaguar, filtered_df_leopard, filtered_df_tiger, filtered_df_snow)
present_sequence_info(filtered_dataframes)



# EXPLORATORY FIG 2: Create a 'vioplot' exploratory figure to showcase distribution of sequence lengths before and after filtering

#This function counts the nucleotides in each sequence for subsequent violin plots
count_nucleotides <- function(df) {
  df %>%
    filter(str_length(COI_sequence) <= 1000) %>% # Exclude whole genome sequences (This is relevant for jaguar dataframe)
    mutate(Sequence_length = width(DNAStringSet(as.character(COI_sequence))))
}

# Apply the function to the original datasets (before filtering)
df_lion_w_counts <- count_nucleotides(df_lion)
df_jaguar_w_counts <- count_nucleotides(df_jaguar)
df_leopard_w_counts <- count_nucleotides(df_leopard)
df_tiger_w_counts <- count_nucleotides(df_tiger)
df_snow_w_counts <- count_nucleotides(df_snow)

sequence_lengths_original <- list(
  Lion = df_lion_w_counts$Sequence_length,
  Jaguar = df_jaguar_w_counts$Sequence_length,
  Leopard = df_leopard_w_counts$Sequence_length,
  Tiger = df_tiger_w_counts$Sequence_length,
  Snow_Leopard = df_snow_w_counts$Sequence_length
)

# Apply the function to the filtered datasets
filtered_df_lion_w_counts <- count_nucleotides(filtered_df_lion)
filtered_df_jaguar_w_counts <- count_nucleotides(filtered_df_jaguar)
filtered_df_leopard_w_counts <- count_nucleotides(filtered_df_leopard)
filtered_df_tiger_w_counts <- count_nucleotides(filtered_df_tiger)
filtered_df_snow_w_counts <- count_nucleotides(filtered_df_snow)

sequence_lengths_filtered <- list(
  Lion = filtered_df_lion_w_counts$Sequence_length,
  Jaguar = filtered_df_jaguar_w_counts$Sequence_length,
  Leopard = filtered_df_leopard_w_counts$Sequence_length,
  Tiger = filtered_df_tiger_w_counts$Sequence_length,
  Snow_Leopard = filtered_df_snow_w_counts$Sequence_length
)

# Plot the combined violin plots
par(mfrow = c(1, 1), cex.axis = 1.2, cex.main = 2) #Adjust the plot and axes titles for optimal size 

vioplot(
  sequence_lengths_original$Lion, sequence_lengths_original$Jaguar,
  sequence_lengths_original$Leopard, sequence_lengths_original$Tiger,
  sequence_lengths_original$Snow_Leopard,
  names = c("Lion", "Jaguar", "Leopard", "Tiger", "Snow Leopard"),
  col = "#56B4E9", border = "black", ylim = c(140, 840),
  main = "Combined Violin Plot of Sequence Lengths Before and After Filtering"
)

# Add the second violin plot to the same axis to showcase comparison
vioplot(
  sequence_lengths_filtered$Lion, sequence_lengths_filtered$Jaguar,
  sequence_lengths_filtered$Leopard, sequence_lengths_filtered$Tiger,
  sequence_lengths_filtered$Snow_Leopard,
  names = c("Lion", "Jaguar", "Leopard", "Tiger", "Snow Leopard"),
  col = "#E41A1C", border = "black", add = TRUE
)

text(2, 680, "*No Change", cex = 1, col = "black", font = 2) 
text(3, 680, "*No Change", cex = 1, col = "black", font = 2) 
text(5, 680, "*No Change", cex = 1, col = "black", font = 2) 


# Manually adjust the size and position of the x and y axes to improve clarity and interpretability of the plot
mtext("Species", side = 1, line = 3, cex = 1.5, font = 2)
mtext("Sequence Length (base pairs)", side = 2, line = 2.5, cex = 1.5, las = 0, font = 2) 

# Add a legend to distinguish between the violin visualizations before filtering vs after (font size and positioning are adjusted as required) 
legend(
  x = 4.72, y = 915, 
  legend = c("Before Filtering", "After Filtering"),
  fill = c("#56B4E9", "#E41A1C"), 
  bty = "n",             # No border around the legend
  cex = 1.37,             # Increase font size
  x.intersp = 0.2,       # Reduce space between text and symbols
  y.intersp = 0.6,
  pt.cex = 1.5,          # Controls the size of the symbols, impacts the box size
  text.font = 2,
  xpd = TRUE
)

# Remove temporary dataframes used for exploratory plots to clean up the environment
rm(df_jaguar_w_counts, df_leopard_w_counts, df_lion_w_counts, df_snow_w_counts, df_tiger_w_counts, filtered_df_jaguar_w_counts, filtered_df_leopard_w_counts, filtered_df_lion_w_counts, filtered_df_snow_w_counts, filtered_df_tiger_w_counts, sequence_lengths_filtered, sequence_lengths_original, summary_country_df)



# 3. MAIN ANALYSIS ----
# MAIN VISUALIZATION 1: Create a map to visualize the geographic distribution of Panthera species

# Create a choropleth map showing species density by country and plot coordinates to visualize range overlap
map <- ggplot() +
  geom_polygon(data = merged_map_density_df, aes(x = long, y = lat, group = group, fill = unique_species_count), color = "grey") +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen", na.value = "grey80") +
  
  # Add species points using 'lat.y' and 'lon' (geom_point) using merged_map_coordinates_df
  geom_point(data = merged_map_coordinates_df, aes(x = lon, y = lat.y, color = species_name), size = 5, alpha = 0.7) + 
  
  scale_color_manual(
    values = c(
      "Panthera onca" = "#E41A1C",  
      "Panthera pardus" = "#984EA3",
      "Panthera tigris" = "#FF7F00",
      "Panthera leo" = "#FFFF33",   
      "Panthera uncia" = "#377EB8"    
    ),
    labels = c(
      "Panthera onca" = "Jaguar",
      "Panthera pardus" = "Leopard",
      "Panthera tigris" = "Tiger",
      "Panthera leo" = "Lion",
      "Panthera uncia" = "Snow Leopard"
    )
  ) +
  
  # Add the legend for the colored species points
  guides(color = guide_legend(title = "Species", title.theme = element_text(face = "bold")), 
         fill = guide_legend(title = "Unique Species Count per Country", title.theme = element_text(face = "bold"))) +  
  
  labs(
    title = expression("Geographic Distribution and Unique Species Count by Country for" ~ italic("Panthera") ~ "Genus"),
    fill = "Unique Species Count per Country",
    color = "Species"
  ) +
  
  theme_void() +
  theme( #The following adjustments are done to improve aesthetics and clarity of visualization
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", family = "Helvetica"),
    legend.position = c(0.17, 0.4),  
    legend.title = element_text(size = 14, face = "bold", hjust = 0, margin = margin(r = 10)),  
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.key.width = unit(0.6, "cm"),  
    legend.key.height = unit(0.6, "cm"), 
    legend.box.margin = margin(10, 10, 10, 10),  
    legend.box.spacing = unit(1.5, "cm") 
  )

print(map) #NOTE: This might around a minute to properly load


# MAIN VISUALIZATION 2: Create a phylogeny to showcase evolutionary relatedness between species in the Panthera genus


#First, create a phylogeny from sequence data

# Initialize empty vector for sequences
sequence_data <- character()

# Add one sequence from each filtered data frame to create the phylogeny.
# Since all sequences have been filtered, any sequence from each data frame can be chosen.
set.seed(119)
for (filtered_df in filtered_dataframes) {
  random_row <- filtered_df[sample(nrow(filtered_df), 1), ]
  random_row_COI <- random_row$COI_sequence
  sequence_data <- c(sequence_data, random_row_COI)
}

sequence_names <- c("Panthera leo (Lion)", 
                    "Panthera onca (Jaguar)", 
                    "Panthera pardus (Leopard)", 
                    "Panthera tigris (Tiger)", 
                    "Panthera uncia (Snow Leopard)")

# Perform multiple sequence alignment using ClustalW method, as this method is commonly used for multi-sequence alignment
alignment <- msa::msa(sequence_data, method = "ClustalW", type = "dna")

aligned_sequences <- msa::msaConvert(alignment, type = "seq") # Convert the alignment to a character matrix

aligned_sequences <- as.DNAbin(aligned_sequences)

# Create the phylogenetic tree using distance matrix (K80 model) as this model is very commonly used in models for DNA sequence evolution
dist_matrix <- dist.dna(aligned_sequences, model = "K80")

# Create and plot the phylogenetic tree using the Neighbor Joining method (this method was chosen because of it's efficiency)
phylo_tree <- ape::nj(dist_matrix)

phylo_tree$tip.label <- sequence_names

par(mar = c(2, 2, 4, 2))  # Adjust margins (bottom, left, top, right)

# Plot the phylogeny
plot(
  phylo_tree,
  main = "Phylogenetic Tree of Panthera Species (COI)", 
  tip.label = gsub("_", " ", phylo_tree$tip.label), 
  cex = 1.5, 
  edge.width = 2,
  label.offset = 0.0005, 
  no.margin = FALSE 
)



# MAIN VISUALIZATION 3: Create a heatmap using the similarity matrix to showcase how similar countries are that species live in

heatmap_legend_param <- list(
  title = "Similarity",           
  title_gp = gpar(fontsize = 12),  
  labels_gp = gpar(fontsize = 10),
  grid_width = unit(20, "mm"),    
  grid_height = unit(8, "mm")    
)

# customize the heatmap for aesthetics and clarity
heatmap <- Heatmap(similarity_matrix, 
                   name = "Similarity",
                   col = colorRamp2(c(0, 1), c("white", "blue")),  
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   width = unit(7, "in"),                          
                   height = unit(7, "in"),                        
                   rect_gp = gpar(col = "black", lwd = 0.5),       
                   heatmap_legend_param = heatmap_legend_param,   
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     # Add text annotations for each cell
                     grid.text(sprintf("%.3f", similarity_matrix[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
                   })

draw(heatmap, heatmap_legend_side = "right")

grid.text("Country Overlap Similarity Among Panthera Species", 
          x = unit(0.5, "npc"),    # Center the title horizontally
          y = unit(0.95, "npc"),   # Position the title slightly below the top of the plotting area
          gp = gpar(fontsize = 16, fontface = "bold"))  # Title font size and style








# BIBLIO
# https://wildcatconservation.org/wild-cats/africa/leopard/
# https://www.worldwildlife.org/species/snow-leopard
# https://www.researchgate.net/figure/Range-of-the-tiger-Panthera-tigris-Historical-range-of-tiger-distribution-is-shown-in_fig1_23785740
# https://www.fws.gov/species/jaguar-panthera-onca
# https://www.worldwildlife.org/species/lion--19

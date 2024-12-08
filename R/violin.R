library(ggplot2)

# Update x-axis labels to include "P. panthera"
sequence_data$ind <- paste("P. panthera", sequence_data$ind)
filtered_sequence_data$ind <- paste("P. panthera", filtered_sequence_data$ind)

# Add a filter status column
sequence_data$Filter_Status <- "Before Filtering"
filtered_sequence_data$Filter_Status <- "After Filtering"

# Combine both datasets into one data frame
combined_data <- rbind(sequence_data, filtered_sequence_data)

# Reorder Filter_Status for correct facet positioning
combined_data$Filter_Status <- factor(combined_data$Filter_Status, levels = c("Before Filtering", "After Filtering"))

# Create the ggplot violin plot with shared x-axis
ggplot(combined_data, aes(x = ind, y = values, fill = Filter_Status)) +
  geom_violin(trim = TRUE, alpha = 0.8, color = "darkblue") +
  facet_grid(Filter_Status ~ ., scales = "free_y") + # Use facet_grid for shared x-axis
  labs(
    title = "Sequence Length Distribution Before and After Filtering",
    x = "Species",
    y = "Sequence Length"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(size = 14, face = "bold"), # Facet labels on the right
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_fill_manual(values = c("Before Filtering" = "lightblue", "After Filtering" = "lightcoral")) +
  scale_y_continuous(limits = c(140, 810)) # Set shared y-axis scale

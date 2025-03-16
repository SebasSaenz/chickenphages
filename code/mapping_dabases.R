# Load library
library(tidyverse)


# Function to extract words from a file
extract_words <- function(file_path, pattern_net, pattern_pattern2, db) {
  # Extract only the filename from the file path
  file_name <- basename(file_path)
  
  # Read the text file
  file_contents <- readLines(file_path)
  
  # Extract the words based on patterns
  extracted_net <- regmatches(file_contents, gregexpr(pattern_net, file_contents, perl = TRUE))
  extracted_pattern2 <- regmatches(file_contents, gregexpr(pattern_pattern2, file_contents, perl = TRUE))
  
  # Create a data frame with extracted words and filename
  df <- data.frame(
    file = file_name,
    data_base = db,
    mapped = as.numeric(unlist(extracted_net)),
    total_reads = as.numeric(unlist(extracted_pattern2))
  )
  
  return(df)
}

# List of file paths
filenames <- list.files(path = "data/map_db/refseq_database/",
                        pattern = "*sh*",
                        full.names = TRUE)

filenames_cm <- list.files(path = "data/map_db/cmkmc_stats/",
                        pattern = "*sh*",
                        full.names = TRUE)

filenames_cpc <- list.files(path = "data/map_db/cpc_database/",
                           pattern = "*sh*",
                           full.names = TRUE)

# Define patterns
pattern_net <- "(?<=found\\s)\\w+"
pattern_pattern2 <- "(?<=out of\\s)\\w+"
db <- "refseq"
# Apply the function to each file and combine the results
results_df <- map_dfr(filenames, ~ extract_words(.x, pattern_net, pattern_pattern2, db))

db <- "cmkmc"
results_df_cm <- map_dfr(filenames_cm, ~ extract_words(.x, pattern_net, pattern_pattern2, db))

db <- "cpc"
results_df_cpc <- map_dfr(filenames_cpc, ~ extract_words(.x, pattern_net, pattern_pattern2, db))

# Display the results

rbind(results_df, results_df_cpc) |> 
  mutate(relative = 100*(mapped/total_reads)) |>  
  ggplot(aes(x = data_base,
             y = relative)) +
  geom_boxplot(width = 0.4, outlier.colour = "white") +
  geom_jitter(width = 0.07, alpha = 0.7, color = "grey", size = 1.5) +
  scale_y_continuous(
    limits = c(0, 14),
                     breaks = seq(0, 14, 2)) +
  labs(y = "Reads mapped (%)",
       x = "") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major.y = element_line(linetype = 2, linewidth = 0.5,),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(),
        text = element_text(size = 16))

ggsave(filename = "plots/comparision_db.png", width = 5, height = 5, dpi = 400)



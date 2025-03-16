# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(here)

# Coverage log dataframe -------------------------------------------------------
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
filenames <- list.files(path = "data/coverage_samples_log/",
                        pattern = "*sh*",
                        full.names = TRUE)

# Define patterns
pattern_net <- "(?<=found\\s)\\w+"
pattern_pattern2 <- "(?<=out of\\s)\\w+"
db <- "cpc"
# Apply the function to each file and combine the results
results_df <- map_dfr(filenames, ~ extract_words(.x, pattern_net, pattern_pattern2, db))

write_tsv(results_df, file = "output_clean/coverage_log.txt")


# vOTUs coverage dataframe------------------------------------------------------

data_join <- list.files(path = "data/coverage_samples/", # Identify all tsv files
                        pattern = "*.txt", full.names = TRUE) %>%
  lapply(read_tsv) %>%                              # iterate read
  reduce(full_join, by = "Contig")

# Save the joined data to a new CSV file
write_tsv(data_join, file = "output_clean/coverm_abundance_votu.tsv")






representatives <- read_tsv("output/representatives_id.txt")

tax_df_family <- read_csv("data/final_prediction_phagcn.csv") #PhageGCN2

filter_phages <- read_tsv("data/filter_phage_class.txt") |> 
  filter(host_bacteria == FALSE) |> 
  pull(class)


process_genomad_files <- function(data_path = "data/genomad_files/", file_pattern = "*.tsv") {
  # List paths
  all_paths <- list.files(
    path = data_path,
    pattern = file_pattern,
    full.names = TRUE
  )
  
  # Join files
  all_content <- lapply(all_paths, read.table, header = TRUE, sep = "\t", encoding = "UTF-8")
  
  # Manipulate file paths
  all_filenames <- basename(all_paths)
  
  # Combine content and filenames
  all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
  
  # Combine into a data table
  all_result <- data.table::rbindlist(all_lists, fill = TRUE)
  
  return(all_result)
}

result <- process_genomad_files(data_path = "data/genomad_files/", file_pattern = "*.tsv")


names(result)[12] <- "sample" # change column name


class_df <- representatives |> 
  mutate(fix = if_else(str_detect(cluster_representative, ".*_1$"), TRUE, FALSE),
         cluster_representative = str_replace(cluster_representative, "_1$", "")) |>
  inner_join(result, by =c("cluster_representative"="seq_name")) |> 
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |> 
  select(cluster_representative, class, fix) |> 
  filter(!class %in% filter_phages) |> 
  pull(cluster_representative)
 
# Fix id for seqkit

fix_id_bacteriophages <- representatives |> 
mutate(fix = if_else(str_detect(cluster_representative, ".*_1$"), TRUE, FALSE),
       cluster_representative = str_replace(cluster_representative, "_1$", "")) |>
  inner_join(result, by =c("cluster_representative"="seq_name")) |> 
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |> 
  select(cluster_representative, class, fix) |> 
  filter(!class %in% filter_phages) |> 
  mutate(cluster_representative = if_else(fix, str_replace(cluster_representative, "$", "_1"), cluster_representative)) |> 
  select(cluster_representative)

write_tsv(fix_id_bacteriophages, file = "output/bacteriophage_id.txt")

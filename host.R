# Load libraries
library(tidyverse)
library(data.table)

# Load dataframes
representatives <- read_tsv("output/representatives_id.txt")

# Load genomad files for taxonomy and number of genes

all_pathpanel.grid.major.x = all_paths <- list.files(
  path = "data/iphop/",
  pattern = "*.csv", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
         header = TRUE,
         sep = ",",
         encoding = "UTF-8"
  )

all_filenames <- all_paths %>% # Manipulate file paths
  basename() %>%
  as.list()

all_lists <- mapply(c,
                    all_content,
                    all_filenames,
                    SIMPLIFY = FALSE
)

all_result <- rbindlist(all_lists, fill = T)

all_result |>
  group_by(Virus) |> 
  slice(1) |>
  ungroup() |> 
  select(Virus, Host.genus) |> 
  separate(Host.genus,
           into = c("domain", "phylum", "class", "order", "family", "genus"),
           sep = ";") |> 
  select(Virus, phylum) |> 
  right_join(representatives, by=c("Virus"="cluster_representative")) |> 
  mutate(phylum =if_else(is.na(phylum), "Unknown", phylum),
         phylum = str_remove(phylum, "p__"),
         phylum = str_remove(phylum, "_[A-Z]")) |> 
  count(phylum) |> 
  ggplot(aes(x = n,
             y = fct_reorder(phylum, n))) +
  geom_col() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 9000),
                     breaks = seq(0, 9000, 2000)) +
  theme_bw()

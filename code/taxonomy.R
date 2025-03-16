# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(data.table)

base_color <- c('#01665e', '#762a83', '#2166ac','#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf')

# Load data --------------------------------------------------------------------
representatives <- read_tsv("output/representatives_id.txt")

tax_df_family <- read_csv("data/final_prediction_phagcn.csv") #PhageGCN2

all_pathpanel.grid.major.x = all_paths <- list.files(
  path = "data/genomad_files/",
  pattern = "*.tsv", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
         header = TRUE,
         sep = "\t",
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
names(all_result)[12] <- "sample" # change column name

# Pie plot Class ---------------------------------------------------------------
class_df <- representatives |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(all_result, by =c("cluster_representative"="seq_name")) |> 
  select(taxonomy) |>  
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |> 
  count(class) |>
  mutate(rel_abund = 100*(n/sum(n)),
         class = if_else(rel_abund < 92, "Other", class)) |> 
  group_by(class) |> 
  summarise(sum_rel = sum(rel_abund))
  
class_df |> 
mutate(ymax = cumsum(sum_rel),
       ymin = c(0, head(ymax, n=-1))) |> 
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=class)) +
  geom_rect(color = "white") +
  coord_polar(theta="y") + 
  scale_fill_manual(values = base_color) +
  xlim(c(2, 4)) +
  theme_void()

ggsave(filename = "plots/figure2/viral_class.pdf", width = 5, height = 5, dpi = 400)

# Bar plot Family --------------------------------------------------------------
x <- tax_df |> 
  mutate(prediction = str_replace(prediction, ".*_like", "Other")) |> 
  count(prediction) |> 
  mutate(precentage = 100*(n/sum(n)))




class_df <- representatives |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(all_result, by =c("cluster_representative"="seq_name")) |> 
  select(cluster_representative, taxonomy) |>  
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |> 
  count(order)

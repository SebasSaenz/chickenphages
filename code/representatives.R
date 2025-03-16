# Load libraries
library(tidyverse)
library(data.table)
library(ggridges)
library(here)

# Load files
clusters <- read_tsv("data/my_clusters.tsv", col_names = c("cluster_representative", "cluster_members"))
checkv <- read_tsv("data/quality_summary-2.tsv")

  
# Save a file with representatives IDs
representatives <- clusters |> 
  select(cluster_representative) #|> 
  write_tsv(file = "output/representatives_id.txt")

# Create a df with the representatives and checkv info
cluster_df <- clusters |> 
  select(cluster_representative) |>
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(checkv, by = c("cluster_representative"="contig_id")) |> 
  mutate(checkv_quality = factor(checkv_quality,
                                 levels = c("Low-quality", "Medium-quality","High-quality", "Complete")))

# Number of vOTUs by quality
cluster_df |>  
  count(checkv_quality) |> 
  ggplot(aes(y = checkv_quality,
         x = n)) +
  geom_col() +
  scale_x_continuous(limits = c(0, 9500),
                     breaks = seq(0, 9500, 1000)) +
  scale_y_discrete(labels = c("High-quality"="High\nquality",
                              "Medium-quality" = "Medium\nquality",
                              "Low-quality" = "Low\nquality")) +
  labs(x = "Number of vOTUs",
       y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 16))


ggsave(file = "plots/figure2/quality.pdf", width = 6, height = 4, dpi = 400)

# Contig size by quality
cluster_df |>
  select(checkv_quality, contig_length) |> 
  mutate(contig_length_mb = contig_length/1000) |> 
  ggplot(aes(y = checkv_quality,
             x = contig_length_mb)) +
  geom_boxplot(outlier.colour = "grey",
               outlier.shape = 1,staplewidth = 0.3,
               width = 0.5) +
  scale_x_continuous(limits = c(0 ,450),
                     breaks = seq(0 ,450, 50)) +
  scale_y_discrete(labels = c("High-quality"="High\nquality",
                              "Medium-quality" = "Medium\nquality",
                              "Low-quality" = "Low\nquality")) +
  labs(y = NULL,
       x = "Lenght (kb)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 16))

ggsave(file = "plots/figure2/size.pdf", width = 6, height = 4, dpi = 400)

# Number of Provirus
cluster_df |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  select(cluster_representative) |> 
  inner_join(checkv, by = c("cluster_representative"="contig_id")) |>
  select(provirus) |> 
  count(provirus) |> 
  mutate(provirus = factor(provirus,
                           levels = c("Yes", "No"))) |> 
  ggplot(aes(x = provirus,
             y = n)) +
  geom_col(width = 0.5) +
  labs(y = "vOTUs",
       x = "Provirus") +
  scale_y_continuous(limits = c(0, 20000),
                     breaks = seq(0, 20000, 2500)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size =14))

ggsave(file = "plots/provirus.png", width = 4, height = 4, dpi = 400)


# Number of viral genes
cluster_df|> 
  select(checkv_quality, viral_genes) |> 
  ggplot(aes(x = viral_genes,
             y = checkv_quality)) +
  geom_boxplot(outlier.colour = "grey",
               outlier.shape = 1,
               staplewidth = 0.3,
               width = 0.5)  +
  labs(y = NULL,
       x = "Viral genes") +
  scale_x_continuous(limits = c(0, 275),
                     breaks = seq(0, 275, 25)) +
  scale_y_discrete(labels = c("High-quality"="High\nquality",
                              "Medium-quality" = "Medium\nquality",
                              "Low-quality" = "Low\nquality")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid = element_line(linetype = 2),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 16))
  
ggsave(file = "plots/figure2/viral_genes.pdf", width = 6, height = 4, dpi = 400)


# Load genomad files for taxonomy and number of genes

all_paths <- list.files(
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

# COunt viral contigs per samples -------
clusters |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(all_result, by =c("cluster_representative"="seq_name")) |> 
  select(taxonomy) |>  
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |>  
  count(class) |> 
  mutate(class = if_else(is.na(order), "Unclassified", class),
         percent = n/sum(n)*100) |> 
  ggplot(aes(y = fct_reorder(class, n),
             x = n)) +
  geom_col() +
  scale_x_log10() +
  labs(y="Viral classes",
       x = "vOTUs") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size =14)) + 
  annotation_logticks(side = "b") 

ggsave(file = "plots/suplementary/viral_classes.png", width = 6, height = 4, dpi = 400)


clusters |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(all_result, by =c("cluster_representative"="seq_name")) |> 
  select(taxonomy) |>  
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |>  
  count(order) |>
  mutate(order = if_else(is.na(order), "Unclassified", order)) |> 
  filter(grepl("les", order) | order == "Unclassified") |> 
  ggplot(aes(y = fct_reorder(order, n),
             x = n)) +
  geom_col() +
  scale_x_log10() +
  labs(y="Viral order",
       x = "vOTUs") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size =14)) + 
  annotation_logticks(side = "b") 

ggsave(file = "plots/viral_order.png", width = 6, height = 4, dpi = 400)

clusters |> 
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) |> 
  inner_join(all_result, by =c("cluster_representative"="seq_name")) |> 
  select(taxonomy) |>  
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") |>  
  count(order, family) |>
  mutate(family = if_else(is.na(family), order, family),
         family = if_else(is.na(family), "Unclassified", family)) |> 
  filter(grepl("ae", family) | family == "Unclassified") |> 
  ggplot(aes(y = fct_reorder(family, n),
             x = n)) +
  geom_col() +
  scale_x_log10() +
  labs(y="Viral family",
       x = "vOTUs") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size =14)) + 
  annotation_logticks(side = "b") 

ggsave(file = "plots/viral_family.png", width = 6, height = 5, dpi = 400)


#Prophages ----------------

x <- clusters %>%  
  mutate(cluster_representative = str_replace(cluster_representative, "_1$", "")) %>%  
  inner_join(all_result, by =c("cluster_representative"="seq_name")) %>% 
  count(topology)

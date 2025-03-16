library(tidyverse)

# Select the ID of dereplciate vOTU refseq

votu_id_refseq <- read_tsv("data/votu_match/my_clusters_refseq_dereplicate.tsv",
                           col_names = c("parent", "cluster")) |> 
  select(parent)

write_tsv(votu_id_refseq, file = "output/votu_id_refseq.txt")

df <- read_tsv("data/votu_match/my_clusters_rvd.tsv", col_names = c("parent", "cluster")) |> 
  mutate(parent = str_replace_all(parent, "\\|" ,"_"),
         cluster = str_replace_all(cluster, "\\|" ,"_"))

cpc <- read_tsv("output/representatives_id.txt") |> 
  mutate(cluster_representative = str_replace(cluster_representative, "\\|" ,"_")) |> 
  pull(cluster_representative)

refseq <- read_tsv("output/votu_id_refseq.txt") |> 
  pull(parent)

id_mgv <- read_tsv("data/votu_match/rvd_id.txt", col_names = "id") |> 
  separate(id,
           into = "id",
           sep = " ") |> 
  mutate(id = str_remove(id, ">"),
         id = str_replace_all(id, "\\|" ,"_")) |> 
  pull(id)

x <- df |> 
  filter(!parent == cluster) |>
  mutate(refseq_parent = if_else(parent %in% id_mgv, TRUE, FALSE),
         cpc_cluster = if_else(str_detect(cluster, paste(cpc, collapse = "|")), TRUE, FALSE),
         cluster_1 = if_else(refseq_parent == TRUE & cpc_cluster == TRUE , TRUE, FALSE),
         cpc_parent = if_else(parent %in% cpc, TRUE, FALSE),
         refseq_cluster = if_else(str_detect(cluster, paste(id_mgv, collapse = "|")), TRUE, FALSE),
         cluster_2 = if_else(cpc_parent == TRUE & refseq_cluster == TRUE, TRUE, FALSE),
  )


x |> 
  count(cluster_2)

x <- df |> 
  #filter(parent == cluster) |> 
  inner_join(cpc, by = c("parent"="cluster_representative")) |> 
  mutate(refseq_cluster = if_else(str_detect(cluster, paste(id_refseq, collapse = "|")), TRUE, FALSE)) |> 
  filter(refseq_cluster == TRUE)

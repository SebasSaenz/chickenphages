# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Load checkv data frame -------------------------------------------------------
df <- read_tsv("data/quality_summary-2.tsv")

# Select best contigs ----------------------------------------------------------

medium_q <- df %>% 
  filter(completeness >= 50,
         viral_genes > 0,
         viral_genes > 3 * host_genes,
         is.na(warnings))

         
low_q <- df %>% 
  filter(checkv_quality == "Low-quality",
         contig_length >= 20000,
         completeness >= 40 ,
         viral_genes > 0,
         viral_genes > 3 * host_genes,
         is.na(warnings))

# Make an ID list
viral_contigs <- rbind(medium_q, low_q) %>%
  select(contig_id)

write_tsv(viral_contigs, file = "output/viral_contigs_hq.txt")

# there is a problem with the IDS need to add _1 to missing ones
id_check <- read_tsv("output/id_check.txt", col_names = "contig_id") %>% 
  mutate(contig_id = str_remove(contig_id, ">"))

commonID<-intersect(viral_contigs$contig_id,id_check$contig_id)

x <- viral_contigs[!viral_contigs$contig_id %in% commonID,]

contigs_missing <- x %>% 
  mutate(contig_id = str_replace(contig_id, "$", "_1"))

write_tsv(contigs_missing, "output/id_missing.txt")


# total viral contigs by quality -----------------------


viral_contigs %>%
  select(checkv_quality) %>% 
  distinct()
  filter(checkv_quality == "Complete") %>% 
  select(completeness_method) %>% 
  count(completeness_method)

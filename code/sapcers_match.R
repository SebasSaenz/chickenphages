# Load library -----------------------------------------------------------------
library(tidyverse)

# Load data --------------------------------------------------------------------

abundance_df <- read_tsv("output_clean/coverm_abundance_votu.tsv")

metadata <- read_csv("data/MAGs_dataset_complete.csv")|> 
  select(Run, isolation_source, geo_loc_name_country)

ileum_meta <- read_tsv("data/PFOWL2-3.tsv")

mag_tax <- read_tsv("data/MAGs_taxonomy_tpm.tsv") %>%  
  select(mag_id=`Bin ID`, Tax, starts_with("TPM"))

gtdb <- read_tsv("data/illeum_gtdb.tsv")

# Wrangle data -----------------------------------------------------------------


gtdb_df <- gtdb %>% 
  select(user_genome, classification) %>% 
  separate(classification,
           into = c("kindom","phylum", "class", "order", "family", "genus", "species"),
           sep = ";")
  


# Pivot MAG abundance 
 mag_abundance <- mag_tax %>% 
  pivot_longer(-c(mag_id, Tax), names_to = "sample", values_to = "TPM") %>% 
  separate(Tax,
           into = c("kindom","phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>%
  select(sample, genus, TPM) %>% 
  mutate(genus = str_remove(genus, "g_"),
         genus = if_else(genus =="", "Unclassified", genus),
         genus = if_else(is.na(genus), "Unclassified", genus),
         sample = str_remove(sample, "TPM ")) %>% 
  group_by(sample, genus) %>% 
  summarise(sum_tpm = sum(TPM), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(mag_rel = 100 * sum_tpm/sum(sum_tpm))
  
mag_abundance %>% 
  filter(genus == "Ligilactobacillus" | genus == "Lactobacillus"| genus == "Limosilactobacillus") %>%
  inner_join(ileum_meta, by =c("sample"="SampleID")) %>% 
  ggplot(aes(x = strain,
             y = mag_rel,
             color = genus)) + 
  geom_boxplot()
  

otu_id_ileum <- abundance_df %>% 
  pivot_longer(cols = !Contig, 
               names_to = "sample", 
               values_to = "abundance") %>% 
  mutate(presence = if_else(abundance > 0, 1, 0)) %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  filter(grepl("UHO", sample)) %>% 
  group_by(Contig) %>% 
  summarise(sum = sum(presence), .groups = "drop") %>% 
  filter(sum > 0) %>% 
  pull(Contig)

ordination_df <- abundance_df %>% 
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  left_join(x, by = c("Contig"="sseqid")) %>% 
  group_by(sample, genus) %>% 
  summarise(sum_rpkm = sum(rpkm)) %>% 
  mutate(rel_abun = 100 * (sum_rpkm / sum(sum_rpkm))) %>% 
  filter(grepl("UHO", sample))



ordination_df %>% 
filter(genus == "Ligilactobacillus" | genus == "Lactobacillus"| genus == "Limosilactobacillus") %>% 
  inner_join(ileum_meta, by =c("sample"="SampleID")) %>% 
  ggplot(aes(x = strain,
             y = rel_abun,
             color = genus)) + 
  geom_boxplot()
  
  









spacers_match <- read_tsv("data/spacers_match.txt") %>% 
  select(qseqid, sseqid) %>% 
  mutate(qseqid = str_remove(qseqid, ".fa_.*"))


x <- spacers_match %>% 
  inner_join(gtdb_df, by = c("qseqid"="user_genome")) %>% 
  select(qseqid, sseqid, genus) %>% 
  mutate(genus = str_remove(genus, "g__"))
  count(genus)


x <- spacers_match %>% 
  count(sseqid)

x <- sapcers_match %>% 
  inner_join(mag_tax, by = c("qseqid"="mag_id")) %>% 
  separate(Tax,
           into = c("kindom","phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>% 
  select(sseqid, genus) %>% 
  mutate(genus = if_else(is.na(genus), "Unknown", genus),
         genus = if_else(genus == "", "Unknown", genus)) 

  left_join(prevalent_votu_ile, by = c("sseqid"="Contig")) %>% 
  select(qseqid, sseqid, genus, prevalence)


  
  
x %>% 
  count(sseqid, genus) %>% 
  left_join(ordination_df, by =c("sseqid"="Contig")) %>% 
  drop_na() %>% 
  inner_join(ileum_meta, by =c("sample"="SampleID")) %>% 
  filter(genus == "g_Ligilactobacillus" | genus == "g_Lactobacillus"| genus == "g_Limosilactobacillus") %>% 
  ggplot(aes(x = strain,
             y = rel_abun,
             color = genus)) + 
  geom_boxplot()
  
  
  
  left_join(mag_abundance, by =c("qseqid"="mag_id", "sample")) %>% 
 

write_tsv(x, file = "clean_data/network.txt") 
  

z <- x %>% 
  count(genus)

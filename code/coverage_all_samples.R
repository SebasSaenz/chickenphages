# Load libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(archive)
library(patchwork)
library(vegan)
library(ggalign)

# Load data frames  ------------------------------------------------------------

coverage_log <- read_tsv("output_clean/coverage_log.txt")

metadata <- read_csv("data/MAGs_dataset_complete.csv")|> 
  select(Run, isolation_source, geo_loc_name_country)

filter_id_virome <-  read_tsv("data/filter_viromen_id.txt") |> 
  pull(run)

tax_df <- read_csv("data/final_prediction_phagcn.csv")

abundance_df <- read_tsv("output_clean/coverm_abundance_votu.tsv")


# Absence presence plot ----------------------------------------------------------------------


presence_plot <- abundance_df %>% 
  pivot_longer(cols = !Contig, 
               names_to = "sample", 
               values_to = "abundance") %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM"),
         presence = if_else(abundance > 0, 1, 0)) %>% 
  #filter(!sample %in% filter_id_virome) %>% 
  group_by(Contig) %>% 
  summarise(sum_presence = sum(presence)) %>% 
  ggplot(aes(x = sum_presence)) +
  geom_histogram(binwidth = 5) +
  scale_x_continuous(breaks = seq(0, 600, 50),
                     expand = c(0, 0)) +
  scale_y_log10(breaks = c(1, 10,100, 1000, 10000),
                guide = guide_axis_logticks(),
                expand =c(0,0)) +
  labs(y = "vOTUs",
       x = "Number of samples") +
  theme_bw() + 
  theme( panel.grid = element_blank(),
         panel.grid.major.y = element_line(linetype = 2, linewidth = 0.2, colour = "grey"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 14),
      )

ggsave(filename = "plots/votus_presence.png",
       width = 7,
       height = 4,
       dpi = 400)  



#Number of samples
number_saples <- abundance_df %>% 
  pivot_longer(cols = !Contig, 
               names_to = "sample", 
               values_to = "abundance") %>% 
        mutate(sample = str_remove(sample, "representative_viral.fasta/"),
               sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  select(sample) %>% 
  unique() %>% 
  inner_join(metadata, by = c("sample"="Run")) %>% 
  count(isolation_source)



# Number of vOTUs per sample -------------------------------
x <- abundance_df %>% 
        pivot_longer(cols = !Contig, 
                     names_to = "sample", 
                     values_to = "abundance") %>% 
        mutate(sample = str_remove(sample, "representative_viral.fasta/"),
               sample = str_remove(sample, "_1.fq.gz RPKM"),
               presence = if_else(abundance > 0, 1, 0)) %>% 
  inner_join(metadata, by = c("sample"="Run")) %>% 
  select(Contig, presence, isolation_source) %>% 
  group_by(Contig, isolation_source) %>% 
  summarise(sum_presence = sum(presence), .groups = "drop") %>% 
  inner_join(number_saples, by="isolation_source") %>% 
  mutate(percent = 100 * (sum_presence/n)) %>% 
  mutate(more_25 = if_else(percent >= 25, 1, 0),
        more_50 = if_else(percent >= 50, 1, 0),
        more_75 = if_else(percent >= 75, 1, 0)) %>% 
  select(isolation_source, more_25, more_50, more_75) %>% 
  pivot_longer(-isolation_source, names_to = "percent", values_to = "value") %>% 
  group_by(isolation_source,percent ) %>% 
  summarise(sum_value = sum(value))


plot_samplae_votu <- x %>% 
  mutate(isolation_source = str_to_sentence(isolation_source),
        isolation_source = factor(isolation_source,
        levels = c("Crop", "Duodenum", "Jejunum", "Ileum", "Caeca", "Colorectum", "Faeces"))) %>% 
  ggplot(aes(x = isolation_source,
             y = sum_value,
            fill = percent)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(limits = c(0 ,225),
                    breaks = seq(0, 225, 25)) +
  scale_fill_manual(breaks = c("more_25", "more_50", "more_75"),
labels = c(">25%", ">50%", ">75%"),
values = c('#66c2a5','#fc8d62','#8da0cb')) +
  labs(x = NULL,
  y = "vOTUs") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.2, colour = "grey"),
    legend.title = element_blank(),
    legend.position = c(0.8, 0.7))




fig_5 <- presence_plot / plot_samplae_votu




ggsave(fig_5, filename = "plots/fig5/votus.png", width = 6, height = 6, dpi = 300)



# ------------------------------------------------------------
core_phages <- x %>% 
  filter(sum_presence >= 300) %>% 
  pull(Contig)



core <- abundance_df %>% 
  pivot_longer(-Contig, names_to = "Run", values_to = "RPKM") %>% 
  mutate(Run = str_remove(Run, "representative_viral.fasta/"),
         Run = str_remove(Run, "_1.fq.gz RPKM")) %>% 
  filter(!Run %in% filter_id_virome) %>% 
  inner_join(metadata, by ="Run") %>%
  group_by(Contig) %>% 
  summarise(mean_rpkm = mean(RPKM), .groups = "drop")
  
core_set <- core %>% 
  filter(mean_rpkm > 1) %>% 
  pull(Contig)


  
  core_df <- abundance_df %>%
  filter(Contig %in% core_set) %>% 
  pivot_longer(-Contig, names_to = "Run", values_to = "RPKM") %>% 
  mutate(Run = str_remove(Run, "representative_viral.fasta/"),
         Run = str_remove(Run, "_1.fq.gz RPKM")) %>% 
  inner_join(metadata, by ="Run") %>%
  #filter(isolation_source != "faeces") %>% 
  mutate(log10_RPKM = log(RPKM+1),
         presence = if_else(RPKM >1, 1, 0),
         phage_source = str_remove(Contig, "_k141.*"),
         isolation_source = factor(isolation_source,
                                   levels = c("crop", "duodenum", "jejunum",
                                               "ileum", "caeca", "colorectum", "faeces"))) %>% 
inner_join(metadata, by=c("phage_source"="Run")) %>% 
  mutate(color_code = case_when(isolation_source.y =="ileum" ~"red"))
  


  core_df %>%
    filter(
      #isolation_source.x == "crop",
      #isolation_source.y == "ileum" | isolation_source.y == "crop"
           ) %>%
    ggplot(aes(x = Run,
             y = Contig,
             fill = presence)) +
  geom_tile() +
  facet_grid(isolation_source.y~geo_loc_name_country.x, scales = "free", space = "free", switch = "x") +
  labs(x = NULL,
       y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none", 
        strip.text.x = element_text(angle = 90, size = 5),
        strip.text.y = element_text(angle = 0, size = 5))

ggsave("plots/heatmap.png", width = 7, height = 7, dpi = 300, units = "in")
  
  
core_df %>% 
  filter(Contig == "SRR6323258_k141_47341") %>% 
  ggplot(aes(y = isolation_source.x,
             x = log10_RPKM)) +
  geom_boxplot()
  
colors_barplot <- rev(c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4'))

df_bar <- core %>% 
  mutate(isolation_country = str_remove(Contig, "_k141.*")) %>% 
  inner_join(metadata, by = c("isolation_country"="Run")) %>% 
  count(isolation_source, geo_loc_name_country)

vector_order <- df_bar %>% 
  group_by(geo_loc_name_country) %>% 
  summarise(total = sum(n)) %>% 
  arrange(total) %>% 
  pull(geo_loc_name_country)


  df_bar %>% 
    mutate(geo_loc_name_country = factor(geo_loc_name_country,
                                         levels = vector_order)) %>% 
    ggplot(aes(y = geo_loc_name_country,
             x = n,
             fill = isolation_source)) +
  geom_col() +
  scale_x_continuous(limits = c(0, 8000),
                     breaks = seq(0, 8000, 1000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colors_barplot) +
  labs(x = "vOTUs",
       y = NULL) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 2),
        panel.grid.minor.x = element_blank())


ggsave("plots/origin_votus.png", width = 6, height = 4, dpi = 300, units = "in")
  
metadata %>% 
  count(geo_loc_name_country, isolation_source) %>% 
  ggplot(aes(x = geo_loc_name_country,
             y = n,
             fill = isolation_source)) +
  geom_col()
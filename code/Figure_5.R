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


color_base <- c('#d73027','#7a0177', '#fc8d59','#fee090','#01665e','#0c2c84','#91bfdb','#4575b4')

# Absence presence plot -------------------------------------------------------------------

read_map <- coverage_log %>% 
  mutate(file = str_remove(file, "_coverM.sh.e.*")) %>% 
  inner_join(metadata, by = c("file" = "Run")) %>% 
  filter(!file %in% filter_id_virome) %>% 
  mutate(isolation_source = factor(isolation_source,
                                   levels = c("crop",
                                              "duodenum",
                                              "jejunum",
                                              "ileum",
                                              "caeca",
                                              "colorectum",
                                              "faeces"),
                                   labels = c("Crop",
                                              "Duodenum",
                                              "Jejunum",
                                              "Ileum",
                                              "Caeca",
                                              "Colorectum",
                                              "Faeces")),
         rel_abun = 100*(mapped/total_reads)) %>% 
  select(file, isolation_source, mapped, total_reads, rel_abun, data_base) %>% 
  #filter(total_reads > 0) %>% 
  filter(isolation_source != "Faeces") %>% 
  ggplot(aes(x = isolation_source,
             y = rel_abun)) +
  geom_violin(aes(fill = data_base, fill = after_scale(colorspace::lighten(fill, .5))),
              size = 1.2, color = NA,
              bw = 0.5
  ) +
  # geom_jitter(width  = 0.1,
  #             color = "grey",
  #             alpha = 0.7,
  #             size = 1,
  #             shape =21) +
  geom_boxplot(width = 0.4,
               alpha = 0.5,
               outlier.colour = "black",
               outlier.size = 1) +
  scale_y_continuous(limits = c(0, 22.5),
                     breaks = seq(0, 22.5, 2.5)) +
  scale_fill_manual(values = c("grey"))+
  labs(x = NULL,
       y = "Reads mapped (%)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.5,),
    panel.grid.major  = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    text = element_text(size = 12),
    axis.text.x = element_text(size = 8)
)

ggsave(filename = "plots/section_viral.png", width = 5, height = 3, dpi = 400)  

 
# Diversity --------------------------------------------------------------------------------------

diversity_df <- abundance_df %>% 
  pivot_longer(cols = !Contig, 
               names_to = "sample", 
               values_to = "abundance") %>% 
  mutate(presence = if_else(abundance > 0, 1, 0)) %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  group_by(sample) %>% 
  summarise(observed = sum(presence), .groups = "drop") %>% 
  left_join(metadata, by = c("sample"="Run")) 


diversity_plot <- diversity_df %>% 
  filter(isolation_source != "faeces") %>%
  mutate(
    isolation_source =str_to_sentence(isolation_source),
    isolation_source = factor(isolation_source,
                                   levels =c("Crop", "Duodenum",
                                             "Jejunum", "Ileum",
                                             "Caeca", "Colorectum"))) %>% 
  ggplot(aes(y = isolation_source,
             x = observed)) +
  geom_violin(aes(fill = isolation_source,
                  fill = after_scale(colorspace::lighten(fill, .2))),
              size = 1.2, color = NA,
              bw = 5) +
  # geom_jitter(width = 0.3,
  #             alpha = 0.8,
  #             aes(color = isolation_source)) +
  geom_boxplot(width = 0.25,           
               alpha = 0,
              outlier.alpha = 1) +
  scale_x_continuous(limits = c(0, 500),
                     breaks = seq(0, 500, 100)) +
  scale_fill_manual(values=color_base) +
  labs(x = NULL,
       y = "Observed vOTUs") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "none",
    legend.box = "horizontal",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  guides(colour = guide_legend(nrow = 1))

ggsave(filename = "plots/observed_votu.png",
       width = 6,
       height = 4,
       dpi = 400)


# Ordination NMDS ------------------------------------------------------------------------------------------------- 

ordination_df <- abundance_df %>% 
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  filter(!sample %in% filter_id_virome) %>% 
   group_by(sample) %>% 
   mutate(
     rel_abun = 100 * (rpkm / sum(rpkm))) %>%
  drop_na() %>% 
   pivot_wider(id_cols = Contig,
               names_from = sample, values_from = rel_abun)

matrix_ordination <- ordination_df %>% 
  select(!Contig) %>% 
  t()


dist <- vegdist(matrix_ordination, method = "bray")

pcoa <- cmdscale(dist, k=2, eig = TRUE, add = TRUE) 
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")

100* pcoa$eig / sum(pcoa$eig)


nmds_plot <- positions %>% 
    as_tibble(rownames="samples") %>% 
    left_join(metadata, by = c("samples"="Run")) %>% 
    filter(isolation_source != "faeces") %>% 
    mutate(
      isolation_source =str_to_sentence(isolation_source),
      isolation_source = factor(isolation_source,
                                levels =c("Crop", "Duodenum",
                                          "Jejunum", "Ileum",
                                          "Colorectum", "Caeca"))) %>% 
    ggplot(aes(x=pcoa1,
               y=pcoa2,
               color = isolation_source)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = color_base) +
    labs(x = "PCo 1 (9.88%)",
         y = "PCo 2 (3.76%)") +
    theme(
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
      legend.box = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
    

# Area plot by section -------------------------------------------------------

bar_color <- c('#ffff99','#b15928','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#a6cee3','#1f78b4')


pivot_abundance <- abundance_df %>% 
  #mutate(Contig = str_replace(Contig, "_1$", "")) %>% 
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>% 
  left_join(host_df, by = c("Contig"="Virus")) %>% 
  select(sample, rpkm, family) %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
                             sample = str_remove(sample, "_1.fq.gz RPKM"),
                             family = if_else(is.na(family), "Unclassified", family),
                             family = str_remove(family, "f__")) %>% 
                filter(!sample %in% filter_id_virome,
                       #prediction != "Unclassified"
                       ) %>% 
                inner_join(metadata, by =c("sample"="Run")) %>% 
                group_by(sample, family, isolation_source) %>% 
                summarise(sum_rpkm = sum(rpkm), .groups = "drop") %>% 
                group_by(sample) %>% 
                mutate(rel_abund = 100*(sum_rpkm/sum(sum_rpkm))) %>% 
                ungroup() %>% 
                drop_na()
              
              
pool_df <- pivot_abundance |>  
  mutate(rel_abund = if_else(is.na(rel_abund), 0, rel_abund)) |> 
  group_by(family) |> 
  summarise(max_abundance = max(rel_abund),
            mean = mean(rel_abund),
            .groups = "drop") |> 
  mutate(pool = if_else(mean < 1, TRUE, FALSE)) |> 
  select(family, pool, mean)
              

sample_order <- pivot_abundance  %>% 
                filter(isolation_source == "colorectum") %>% 
                inner_join(pool_df, by = "family") %>% 
                filter(family == "Unclassified") %>% 
                arrange(desc(rel_abund)) %>%
                mutate(order_samples = 1:nrow(.)) %>% 
                select(sample, order_samples)
              
              
              
colo_plot <- pivot_abundance |> 
                filter(isolation_source == "colorectum") |> 
                inner_join(pool_df, by = "family") |> 
                mutate(family = if_else(pool, "Other", family)) |> 
                group_by(sample, isolation_source, family) |> 
                summarise(sum_abundance = sum(rel_abund),
                          mean = sum(mean),
                          .groups = "drop") |> 
                mutate(family = factor(family),
              family = fct_reorder(family, mean)) |> 
                inner_join(sample_order, by = "sample") |> 
                mutate(sample = factor(sample),
                       sample = fct_reorder(sample, order_samples)) |> 
  #filter(order == "Unclassified") |>
                ggplot(aes(x = sample,
                           y = sum_abundance,
                           fill = family)) +
                geom_col(width  = 1) +
                scale_y_continuous(expand = c(0, 0),
              breaks = seq(0, 100, 20)) +
                scale_fill_manual(values = bar_color) +
                labs(y = "Relative abundance (%)",
                     x = NULL,
                    title = "Colorectum") +
                theme_bw() + 
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      legend.title = element_blank(),
                      legend.text = element_text(size = 8),
                    title = element_text(size = 8, margin = margin(0,0,0,0, "cm")))


areaplot <- crop_plot + duo_plot + jeju_plot + ile_plot + cae_plot + colo_plot  + 
  plot_layout(guides = 'collect', axes = 'collect') &
  theme(legend.key.size = unit(0.4, "cm"))


ggsave(filename = "plots/figure6/areplot.png", width = 8, height = 4)





# Compose plot ---------------------------------------------------------------------------------------------

figure_6 <- (diversity_plot | nmds_plot) /  areaplot +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
  theme(axis.title.y = element_text(size = 12))

ggsave(figure_6, filename = "plots/figure6/figure_6.png", width = 8, height = 6, dpi = 300)

library(tidyverse)

coverage_log <- read_tsv("output_clean/coverage_log.txt")

metadata <- read_csv("data/MAGs_dataset_complete.csv")|> 
  select(Run, isolation_source, geo_loc_name_country)

filter_id_virome <-  read_tsv("data/filter_viromen_id.txt") |> 
  pull(run)


abundance_df <- read_tsv("output_clean/coverm_abundance_votu.tsv")


filter_meta <- metadata %>% 
  filter(Run %in% filter_id_virome)


cover <- coverage_log %>% 
  mutate(file = str_remove(file, "_coverM.sh.e.*")) %>%
  filter(file %in% filter_id_virome) %>% 
  mutate(percent = 100 * (mapped/total_reads))


x <- abundance_df %>% 
  pivot_longer(cols = !Contig, 
               names_to = "sample", 
               values_to = "abundance") %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM"),
         presence = if_else(abundance > 0, 1, 0)) %>% 
  filter(sample %in% filter_id_virome) %>% 
  group_by(Contig) %>% 
  summarise(sum_presence = sum(presence)) %>% 
  ggplot(aes(x = sum_presence)) +
  geom_histogram(binwidth = 5) +
  scale_x_continuous(breaks = seq(0, 650, 25),
                     expand = c(0, 0)) +
  scale_y_log10(breaks = c(1, 10,100, 1000, 10000),
                guide = guide_axis_logticks(),
                expand =c(0,0)) +
  labs(y = "vOTUs",
       x = "Number of samples") +
  theme_bw() + 
  theme(panel.grid = element_line(linetype = 2),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

ggsave(filename = "plots/votus_presence.png",
       width = 7,
       height = 4,
       dpi = 400)  





ordination_df <- abundance_df %>% 
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>% 
  mutate(sample = str_remove(sample, "representative_viral.fasta/"),
         sample = str_remove(sample, "_1.fq.gz RPKM")) %>% 
  filter(sample %in% filter_id_virome) %>% 
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

color_base <- c('#d73027','#7a0177', '#fc8d59','#fee090','#01665e','#0c2c84','#91bfdb','#4575b4')

positions %>% 
  as_tibble(rownames="samples") %>% 
  left_join(filter_meta, by = c("samples"="Run")) %>% 
  ggplot(aes(x=pcoa1,
             y=pcoa2,
             color = isolation_source)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = color_base) +
  labs(x = "PCo 1 (9.88%)",
       y = "PCo 2 (3.76%)") +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

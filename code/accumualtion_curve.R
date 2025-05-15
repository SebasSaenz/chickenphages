# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(vegan)
library(patchwork)

# Load data --------------------------------------------------------------------
abundance_df <- read_tsv("output_clean/coverm_abundance_votu.tsv")

filter_id_virome <- read_tsv("data/filter_viromen_id.txt") |>
  pull(run)

genus_cluster <- read.table(
  "data/genus_clusters.txt",
  header = FALSE,
  sep = "\t",
  fill = TRUE
)

family_cluster <- read.table(
  "data/family_clusters.txt",
  header = FALSE,
  sep = "\t",
  fill = TRUE
)


# Wrangle data -----------------------------------------------------------------

# Saturation curve vOTUs
presence <- abundance_df %>%
  pivot_longer(cols = !Contig, names_to = "sample", values_to = "abundance") %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM"),
    presence = if_else(abundance > 0, 1, 0)
  )

df <- presence %>%
  pivot_wider(id_cols = Contig, names_from = sample, values_from = presence) %>%
  column_to_rownames(var = "Contig") %>%
  t()

accum_curve <- specaccum(df, method = "random", permutations = 1000)

accum_species_plot <- data.frame(
  samples = accum_curve$sites,
  richness = accum_curve$richness,
  sd = accum_curve$sd
) %>%
  # Calculate confidence intervals
  mutate(lower = richness - sd, upper = richness + sd) %>%
  ggplot(aes(x = samples, y = richness)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "skyblue") +
  scale_x_continuous(limits = c(0, 1500)) +
  labs(
    x = NULL,
    y = "Cumulative species-level vOTUs",
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
# Saturation curve genus

genus_df <- genus_cluster %>%
  mutate(cluster = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "position",
    values_to = "Contig"
  ) %>%
  select(cluster, Contig) %>%
  mutate(Contig = str_remove(Contig, "_CDS")) %>%
  filter(!is.na(Contig) & Contig != "") %>%
  left_join(presence, by = "Contig") %>%
  group_by(sample, cluster) %>%
  summarise(sum_presence = sum(presence), .groups = "drop") %>%
  mutate(presence = if_else(sum_presence > 0, 1, 0)) %>%
  pivot_wider(
    id_cols = cluster,
    names_from = sample,
    values_from = presence
  ) %>%
  column_to_rownames(var = "cluster") %>%
  t()

accum_curve_genus <- specaccum(genus_df, method = "random", permutations = 1000)


accum_genus_plot <- data.frame(
  samples = accum_curve_genus$sites,
  richness = accum_curve_genus$richness,
  sd = accum_curve_genus$sd
) %>%
  # Calculate confidence intervals
  mutate(lower = richness - sd, upper = richness + sd) %>%
  ggplot(aes(x = samples, y = richness)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "skyblue") +
  scale_x_continuous(limits = c(0, 1500)) +
  labs(
    x = NULL,
    y = "Cumulative genus-level vOTUs",
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )


# Saturation curve family-level
family_df <- family_cluster %>%
  mutate(cluster = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "position",
    values_to = "Contig"
  ) %>%
  select(cluster, Contig) %>%
  mutate(Contig = str_remove(Contig, "_CDS")) %>%
  filter(!is.na(Contig) & Contig != "") %>%
  left_join(presence, by = "Contig") %>%
  group_by(sample, cluster) %>%
  summarise(sum_presence = sum(presence), .groups = "drop") %>%
  mutate(presence = if_else(sum_presence > 0, 1, 0)) %>%
  pivot_wider(
    id_cols = cluster,
    names_from = sample,
    values_from = presence
  ) %>%
  column_to_rownames(var = "cluster") %>%
  t()

accum_curve_family <- specaccum(
  family_df,
  method = "random",
  permutations = 1000
)


accum_family_plot <- data.frame(
  samples = accum_curve_family$sites,
  richness = accum_curve_family$richness,
  sd = accum_curve_family$sd
) %>%
  # Calculate confidence intervals
  mutate(lower = richness - sd, upper = richness + sd) %>%
  ggplot(aes(x = samples, y = richness)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "skyblue") +
  scale_x_continuous(limits = c(0, 1500)) +
  labs(
    x = "Number of Samples",
    y = "Cumulative family-level vOTUs",
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# Compose plot -----------------------------------------------------------------

accum_species_plot /
  accum_genus_plot /
  accum_family_plot +
  plot_annotation(tag_levels = 'A')

ggsave(
  filename = "plots/accumulation_curves.png",
  dpi = 300,
  width = 5,
  height = 7
)

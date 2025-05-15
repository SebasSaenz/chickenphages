# Load libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(archive)
library(patchwork)
library(vegan)
library(ggalign)
library(ggtext)

# Load data frames  ------------------------------------------------------------

coverage_log <- read_tsv("output_clean/coverage_log.txt")

ileum_meta <- read_tsv("data/PFOWL2-3.tsv")

metadata <- read_csv("data/MAGs_dataset_complete.csv") |>
  select(Run, isolation_source, geo_loc_name_country)

filter_id_virome <- read_tsv("data/filter_viromen_id.txt") |>
  pull(run)

tax_df <- read_csv("data/final_prediction_phagcn.csv")

abundance_df <- read_tsv("output_clean/coverm_abundance_votu.tsv")


color_base <- c(
  '#d73027',
  '#7a0177',
  '#fc8d59',
  '#fee090',
  '#01665e',
  '#0c2c84',
  '#91bfdb',
  '#4575b4'
)

# Absence presence plot -------------------------------------------------------------------

read_map <- coverage_log %>%
  mutate(file = str_remove(file, "_coverM.sh.e.*")) %>%
  inner_join(metadata, by = c("file" = "Run")) %>%
  filter(!file %in% filter_id_virome) %>%
  mutate(
    isolation_source = factor(
      isolation_source,
      levels = c(
        "crop",
        "duodenum",
        "jejunum",
        "ileum",
        "caeca",
        "colorectum",
        "faeces"
      ),
      labels = c(
        "Crop",
        "Duodenum",
        "Jejunum",
        "Ileum",
        "Caeca",
        "Colorectum",
        "Faeces"
      )
    ),
    rel_abun = 100 * (mapped / total_reads)
  ) %>%
  select(file, isolation_source, mapped, total_reads, rel_abun, data_base) %>%
  #filter(total_reads > 0) %>%
  filter(isolation_source != "Faeces") %>%
  ggplot(aes(x = isolation_source, y = rel_abun)) +
  geom_violin(
    aes(fill = data_base, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2,
    color = NA,
    bw = 0.5
  ) +
  # geom_jitter(width  = 0.1,
  #             color = "grey",
  #             alpha = 0.7,
  #             size = 1,
  #             shape =21) +
  geom_boxplot(
    width = 0.4,
    alpha = 0.5,
    outlier.colour = "black",
    outlier.size = 1
  ) +
  scale_y_continuous(limits = c(0, 22.5), breaks = seq(0, 22.5, 2.5)) +
  scale_fill_manual(values = c("grey")) +
  labs(x = NULL, y = "Reads mapped (%)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.5, ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    text = element_text(size = 12),
    axis.text.x = element_text(size = 8)
  )

ggsave(filename = "plots/section_viral.png", width = 5, height = 3, dpi = 400)


# Diversity --------------------------------------------------------------------------------------

diversity_df <- abundance_df %>%
  pivot_longer(cols = !Contig, names_to = "sample", values_to = "abundance") %>%
  mutate(presence = if_else(abundance > 0, 1, 0)) %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM")
  ) %>%
  group_by(sample) %>%
  summarise(observed = sum(presence), .groups = "drop") %>%
  left_join(metadata, by = c("sample" = "Run"))


diversity_plot <- diversity_df %>%
  filter(isolation_source != "faeces") %>%
  mutate(
    isolation_source = str_to_sentence(isolation_source),
    isolation_source = factor(
      isolation_source,
      levels = rev(c(
        "Crop",
        "Duodenum",
        "Jejunum",
        "Ileum",
        "Caeca",
        "Colorectum"
      ))
    )
  ) %>%
  ggplot(aes(y = isolation_source, x = observed)) +
  geom_violin(
    aes(
      fill = isolation_source,
      fill = after_scale(colorspace::lighten(fill, .2))
    ),
    size = 1.2,
    color = NA,
    bw = 5
  ) +
  # geom_jitter(width = 0.3,
  #             alpha = 0.8,
  #             aes(color = isolation_source)) +
  geom_boxplot(width = 0.25, alpha = 0, outlier.alpha = 1) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  scale_fill_manual(values = color_base) +
  labs(y = NULL, x = "Observed vOTUs") +
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

ggsave(filename = "plots/observed_votu.png", width = 6, height = 4, dpi = 400)


observed_diet_plot <- diversity_df %>%
  filter(grepl("UHO", sample)) %>%
  inner_join(ileum_meta, by = c("sample" = "SampleID")) %>%
  ggplot(aes(x = diet, y = observed)) +
  # geom_jitter(size = 2, width = 0.3, color = "grey") +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 280, 40)) +
  geom_line(
    data = tibble(x = c(1, 2), y = c(260, 260)),
    aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tibble(x = 1.5, y = 275),
    aes(x = x, y = y, label = "*"),
    inherit.aes = FALSE
  ) +
  labs(y = " Number of observed vOTUs", x = NULL) +
  theme_bw()

kruskal.test(observed ~ strain, data = observed_strain_plot)

observed_strain_plot <- diversity_df %>%
  filter(grepl("UHO", sample)) %>%
  inner_join(ileum_meta, by = c("sample" = "SampleID")) %>%
  ggplot(aes(x = strain, y = observed)) +
  # geom_jitter(size = 2, width = 0.3, color = "grey") +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 280, 40)) +
  geom_line(
    data = tibble(x = c(1, 2), y = c(260, 260)),
    aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tibble(x = 1.5, y = 275),
    aes(x = x, y = y, label = "***"),
    inherit.aes = FALSE
  ) +
  labs(y = " Number of observed vOTUs", x = NULL) +
  theme_bw()


library(MASS)
glm_interaction_nb <- glm.nb(
  observed ~ strain * diet,
  data = observed_strain_plot
)

summary(glm_interaction_nb)

plot(glm_interaction$residuals)


anova(glm_interaction, glm_model, test = "Chisq")

diversity_df %>%
  filter(grepl("UHO", sample)) %>%
  inner_join(metadata, by = c("sample" = "Run")) %>%
  ggplot(aes(x = Host_age, y = observed)) +
  geom_jitter(aes(color = Host_age), size = 2, width = 0.3) +
  geom_boxplot(alpha = 0, width = 0.5) +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 280, 30)) +
  facet_grid(~`feed_addictive_/growth_promoter`) +
  theme_classic()


ggsave(filename = "plots/test_ileum.png", width = 3, height = 5, dpi = 300)


# Diversity within samples -----------------------------------------------------

otu_id_ileum <- abundance_df %>%
  pivot_longer(cols = !Contig, names_to = "sample", values_to = "abundance") %>%
  mutate(presence = if_else(abundance > 0, 1, 0)) %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM")
  ) %>%
  filter(grepl("UHO", sample)) %>%
  group_by(Contig) %>%
  summarise(sum = sum(presence), .groups = "drop") %>%
  filter(sum > 0) %>%
  pull(Contig)

diversity_jaccard <- abundance_df %>%
  pivot_longer(cols = !Contig, names_to = "sample", values_to = "abundance") %>%
  mutate(presence = if_else(abundance > 0, 1, 0)) %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM")
  ) %>%
  filter(grepl("UHO", sample), Contig %in% otu_id_ileum) %>%
  pivot_wider(id_cols = Contig, names_from = sample, values_from = presence)

# Function to calculate shared OTUs between sample pairs and return a table
shared_otus_table <- function(df) {
  # Remove first column (OTU IDs) to get only sample data
  sample_data <- df[, -1]

  # Get sample names
  sample_names <- colnames(sample_data)

  # Generate sample pairs
  sample_pairs <- combn(sample_names, 2, simplify = FALSE)

  # Create empty result list
  results <- list()

  for (i in seq_along(sample_pairs)) {
    pair <- sample_pairs[[i]]
    shared_count <- sum(sample_data[[pair[1]]] > 0 & sample_data[[pair[2]]] > 0)

    results[[i]] <- data.frame(
      Sample1 = pair[1],
      Sample2 = pair[2],
      Shared_OTUs = shared_count
    )
  }

  # Combine results into a single dataframe
  shared_df <- do.call(rbind, results)

  return(shared_df)
}


# Generate shared OTUs table
shared_otu_table <- shared_otus_table(diversity_jaccard)


strain_type <- ileum_meta %>%
  dplyr::select(SampleID, strain)

shared_diet_plot <- shared_otu_table %>%
  inner_join(strain_type, by = c("Sample1" = "SampleID")) %>%
  inner_join(strain_type, by = c("Sample2" = "SampleID")) %>%
  mutate(
    result = if_else(diet.x == diet.y, "Same<br>diet", "Different<br>diet")
  ) %>%
  ggplot(aes(x = result, y = Shared_OTUs)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(0, 220)) +
  geom_line(
    data = tibble(x = c(1, 2), y = c(200, 200)),
    aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tibble(x = 1.5, y = 210),
    aes(x = x, y = y, label = "***"),
    inherit.aes = FALSE
  ) +
  labs(y = "Number of shared vOTUs", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_markdown())

kruskal.test(Shared_OTUs ~ result, data = shared_diet_plot)

shared_strain_plot <- shared_otu_table %>%
  inner_join(strain_type, by = c("Sample1" = "SampleID")) %>%
  inner_join(strain_type, by = c("Sample2" = "SampleID")) %>%
  mutate(
    result = if_else(
      strain.x == strain.y,
      "Same<br>strain",
      "Different<br>strain"
    )
  ) %>%
  ggplot(aes(x = result, y = Shared_OTUs)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(0, 220)) +
  geom_line(
    data = tibble(x = c(1, 2), y = c(200, 200)),
    aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tibble(x = 1.5, y = 210),
    aes(x = x, y = y, label = "***"),
    inherit.aes = FALSE
  ) +
  labs(y = "Number of shared vOTUs", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_markdown())

kruskal.test(Shared_OTUs ~ result, data = shared_strain_plot)

library(rstatix)

kruskal_result <- kruskal.test(Shared_OTUs ~ result, data = shared_strain_plot)
kruskal_effsize(Shared_OTUs ~ result, data = shared_strain_plot)


glm_model <- glm(
  Shared_OTUs ~ result,
  data = shared_strain_plot,
  family = Gamma(link = "log")
)
summary(glm_model)


illeum_votu <- (observed_strain_plot +
  observed_diet_plot +
  shared_strain_plot) +
  shared_diet_plot +
  plot_layout(
    ncol = 2,
    axis_titles = "collect_y",
    guides = 'collect',
    axes = 'collect'
  )


# vOTUs present in breed

prevalent_votu_ile <- diversity_jaccard %>%
  inner_join(ileum_meta, by = c("sample" = "SampleID")) %>%
  select(Contig, sample, strain, presence) %>%
  group_by(strain, Contig) %>%
  summarise(sum_presence = sum(presence), .groups = "drop") %>%
  pivot_wider(names_from = strain, values_from = sum_presence) %>%
  mutate(
    prevalence_LSL = LB < 108 & LSL > 108,
    prevalence_LB = LB > 108 & LSL < 108,
    prevalence = case_when(
      prevalence_LSL == TRUE & prevalence_LB == FALSE ~ "Prevalent in LSL",
      prevalence_LSL == FALSE & prevalence_LB == TRUE ~ "Prevalent in LB"
    )
  )

ileum_meta %>%
  count(strain)


# Taxonomy vOTUs illeum

illeum_abundance <- abundance_df %>%
  filter(Contig %in% otu_id_ileum) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM")
  ) %>%
  filter(!sample %in% filter_id_virome) %>%
  filter(grepl("UHO", sample)) %>%
  left_join(ileum_meta, by = c("sample" = "SampleID")) %>%
  left_join(host_df, by = c("Contig" = "Virus")) %>%
  mutate(
    genus = if_else(is.na(genus), "Unclassified", genus),
    genus = str_remove(genus, "g__")
  ) %>%
  group_by(sample, genus, strain) %>%
  summarise(sum_rpkm = sum(rpkm), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(rel_abund = 100 * (sum_rpkm / sum(sum_rpkm))) %>%
  ungroup() %>%
  drop_na()


ileum_abunda_votu <- illeum_abundance %>%
  filter(
    genus == "Ligilactobacillus" |
      genus == "Limosilactobacillus" |
      genus == "Lactobacillus"
  ) %>%
  ggplot(aes(x = strain, y = rel_abund)) +
  geom_jitter(
    aes(color = genus),
    position = position_jitterdodge(jitter.width = 0.2),
    size = 0.8,
    shape = 21
  ) +
  geom_boxplot(
    aes(group = interaction(strain, genus)),
    alpha = 0,
    outlier.alpha = 0,
    width = 0.7
  ) +
  coord_cartesian(expand = FALSE) +
  scale_y_continuous(limits = c(0, 95), breaks = seq(0, 90, 15)) +
  scale_color_manual(
    values = c('#a6611a', '#dfc27d', '#80cdc1', 'grey'),
    labels = c(
      expression(italic("Lactobacillus")),
      expression(italic("Ligilactobacillus")),
      expression(italic("Limosilactobacillus"))
    )
  ) +
  labs(x = NULL, y = "Relative abundance %") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = 2),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.key.size = unit(1, "pt"),
    legend.justification = c(1, 1)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) # Makes legend points bigger


ggsave(
  filename = "plots/strain_host_abundance.png",
  dpi = 300,
  width = 3,
  height = 3
)


pool_df <- illeum_abundance |>
  mutate(rel_abund = if_else(is.na(rel_abund), 0, rel_abund)) |>
  group_by(genus) |>
  summarise(
    max_abundance = max(rel_abund),
    mean = mean(rel_abund),
    .groups = "drop"
  ) |>
  mutate(pool = if_else(mean < 0.4, TRUE, FALSE)) |>
  select(genus, pool, mean)


sample_order <- illeum_abundance %>%
  filter(strain == "LSL") %>%
  inner_join(pool_df, by = "genus") %>%
  filter(genus == "Unclassified") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order_samples = 1:nrow(.)) %>%
  select(sample, order_samples)


illeum_abundance |>
  filter(strain == "LSL") |>
  inner_join(pool_df, by = "genus") |>
  mutate(genus = if_else(pool, "Other", genus)) |>
  group_by(sample, strain, genus) |>
  summarise(
    sum_abundance = sum(rel_abund),
    mean = sum(mean),
    .groups = "drop"
  ) |>
  mutate(genus = factor(genus), genus = fct_reorder(genus, mean)) |>
  inner_join(sample_order, by = "sample") |>
  mutate(
    sample = factor(sample),
    sample = fct_reorder(sample, order_samples)
  ) |>
  #filter(order == "Unclassified") |>
  ggplot(aes(x = sample, y = sum_abundance, fill = genus)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = bar_color) +
  labs(y = "Relative abundance (%)", x = NULL, title = "LSL") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    title = element_text(size = 8, margin = margin(0, 0, 0, 0, "cm"))
  )


ggsave(filename = "plots/illeum_tax_lsl.png", dpi = 300, height = 3, width = 6)
# Ordination NMDS --------------------------------------------------------------

ordination_df <- abundance_df %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM")
  ) %>%
  filter(!sample %in% filter_id_virome) %>%
  group_by(sample) %>%
  mutate(
    rel_abun = 100 * (rpkm / sum(rpkm))
  ) %>%
  drop_na() %>%
  pivot_wider(id_cols = Contig, names_from = sample, values_from = rel_abun)

matrix_ordination <- ordination_df %>%
  select(!Contig) %>%
  t()


dist <- vegdist(matrix_ordination, method = "bray")

pcoa <- cmdscale(dist, k = 2, eig = TRUE, add = TRUE)
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")

100 * pcoa$eig / sum(pcoa$eig)


nmds_plot <- positions %>%
  as_tibble(rownames = "samples") %>%
  left_join(metadata, by = c("samples" = "Run")) %>%
  filter(isolation_source != "faeces") %>%
  mutate(
    isolation_source = str_to_sentence(isolation_source),
    isolation_source = factor(
      isolation_source,
      levels = c("Crop", "Duodenum", "Jejunum", "Ileum", "Colorectum", "Caeca")
    )
  ) %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = isolation_source)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_color_manual(values = color_base) +
  labs(x = "PCo 1 (9.88%)", y = "PCo 2 (3.76%)") +
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

bar_color <- c(
  '#ffff99',
  '#b15928',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#a6cee3',
  '#1f78b4'
)


pivot_abundance <- abundance_df %>%
  #mutate(Contig = str_replace(Contig, "_1$", "")) %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "rpkm") %>%
  left_join(host_df, by = c("Contig" = "Virus")) %>%
  select(sample, rpkm, family) %>%
  mutate(
    sample = str_remove(sample, "representative_viral.fasta/"),
    sample = str_remove(sample, "_1.fq.gz RPKM"),
    family = if_else(is.na(family), "Unclassified", family),
    family = str_remove(family, "f__")
  ) %>%
  filter(
    !sample %in% filter_id_virome,
    #prediction != "Unclassified"
  ) %>%
  inner_join(metadata, by = c("sample" = "Run")) %>%
  group_by(sample, family, isolation_source) %>%
  summarise(sum_rpkm = sum(rpkm), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(rel_abund = 100 * (sum_rpkm / sum(sum_rpkm))) %>%
  ungroup() %>%
  drop_na()


pool_df <- pivot_abundance |>
  mutate(rel_abund = if_else(is.na(rel_abund), 0, rel_abund)) |>
  group_by(family) |>
  summarise(
    max_abundance = max(rel_abund),
    mean = mean(rel_abund),
    .groups = "drop"
  ) |>
  mutate(pool = if_else(mean < 1, TRUE, FALSE)) |>
  select(family, pool, mean)


sample_order <- pivot_abundance %>%
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
  summarise(
    sum_abundance = sum(rel_abund),
    mean = sum(mean),
    .groups = "drop"
  ) |>
  mutate(family = factor(family), family = fct_reorder(family, mean)) |>
  inner_join(sample_order, by = "sample") |>
  mutate(
    sample = factor(sample),
    sample = fct_reorder(sample, order_samples)
  ) |>
  #filter(order == "Unclassified") |>
  ggplot(aes(x = sample, y = sum_abundance, fill = family)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = bar_color) +
  labs(y = "Relative abundance (%)", x = NULL, title = "Colorectum") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    title = element_text(size = 8, margin = margin(0, 0, 0, 0, "cm"))
  )


areaplot <- crop_plot +
  duo_plot +
  jeju_plot +
  ile_plot +
  cae_plot +
  colo_plot +
  plot_layout(guides = 'collect', axes = 'collect') &
  theme(legend.key.size = unit(0.4, "cm"))


ggsave(filename = "plots/figure6/areplot.png", width = 8, height = 4)


# Compose plot ---------------------------------------------------------------------------------------------
# figure_6 <- (diversity_plot | nmds_plot) /  areaplot +
#   plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
#   theme(axis.title.y = element_text(size = 12))

figure_6 <- (diversity_plot + plot_spacer()) /
  areaplot +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(axis.title.y = element_text(size = 8))

ggsave(
  figure_6,
  filename = "plots/figure6/figure_6.png",
  width = 8,
  height = 6,
  dpi = 300
)

wrap_ileum_abunda_votu <- wrap_elements(ileum_abunda_votu)
wrap_illeum_votu <- wrap_elements(illeum_votu)
figure_7 <- (wrap_illeum_votu | wrap_ileum_abunda_votu) +
  plot_annotation(tag_levels = "A")

ggsave(
  figure_7,
  filename = "plots/figure_7.png",
  width = 6,
  height = 5,
  dpi = 300
)

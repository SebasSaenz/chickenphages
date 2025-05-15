# Load libraries and data frames
library(tidyverse)
library(ggrepel)
library(patchwork)

style_color <- c("#01665e", "#762a83", "#2166ac")

tax_color <- c('#7570b3', '#e7298a', '#66a61e', '#1b9e77', '#d95f02')

votu <- read_tsv("output/representatives_id.txt") |>
  rename(Virus = cluster_representative)

# Merge files
data_join <- list.files(
  path = "data/iphop/", # Identify all tsv files
  pattern = "*.csv",
  full.names = TRUE
) %>%
  lapply(read_csv) %>% # iterate read
  reduce(rbind) |>
  rename(host_tax = `Host genus`)

# Clean data and add unknowm host
host_df <- data_join |>
  group_by(Virus) |>
  slice(1) |>
  ungroup() |>
  right_join(votu, by = c("Virus")) |>
  select(Virus, host_tax) |>
  separate(
    host_tax,
    into = c("domain", "phylum", "class", "order", "family", "genus"),
    sep = ";"
  )

host_df |>
  select(domain) |>
  mutate(
    domain = str_remove(domain, "d__"),
    domain = if_else(is.na(domain), "Unknown", domain)
  ) |>
  count(domain) |>
  mutate(
    rel_abund = 100 * (n / sum(n)),
    ymax = cumsum(rel_abund),
    ymin = c(0, head(ymax, n = -1))
  ) |>
  ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = domain)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  #scale_fill_manual(values = base_color) +
  xlim(c(2, 4)) +
  theme_void()

ggsave(filename = "plots/figure2/host.pdf", width = 5, height = 5, dpi = 400)


host_df |>
  select(phylum) |>
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = if_else(is.na(phylum), "Unknown", phylum),
    phylum = str_remove(phylum, "_[A-Z]")
  ) |>
  count(phylum) |>
  mutate(phylum = if_else(n < 100, "Others", phylum)) |>
  group_by(phylum) |>
  summarise(sum_n = sum(n)) |>
  mutate(
    rel_abund = 100 * (sum_n / sum(sum_n)),
    csum = rev(cumsum(rev(rel_abund))),
    pos = rel_abund / 2 + lead(csum, 1),
    pos = if_else(is.na(pos), rel_abund / 2, pos)
  ) |>
  ggplot(aes(x = 2, y = rel_abund, fill = phylum)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  geom_text_repel(
    aes(
      label = paste0(round(rel_abund, 2), "%"),
      y = pos
    ),
    nudge_y = 2,
    nudge_x = 0.5
  ) +
  theme_void() +
  theme()

##### Figure 3

genome_size <- read_tsv("data/quality_summary-2.tsv")


tax_votu <- read_csv("data/final_prediction_phagcn.csv")
life_style <- read_tsv("data/representative_viral.fasta_bacphlip.txt") |>
  rename(Virus = `...1`)

host_counts <- host_df |>
  select(phylum) |>
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = str_remove(phylum, "_[A-Z]"),
    phylum = if_else(is.na(phylum), "Unknown", phylum)
  ) |>
  count(phylum)


host_pool_df <- host_counts |>
  mutate(host_pool = if_else(n < 6, TRUE, FALSE)) |>
  select(phylum, host_pool)

host_plot <- host_counts |>
  inner_join(host_pool_df, by = "phylum") |>
  mutate(phylum = if_else(host_pool, "Other", phylum)) |>
  group_by(phylum) |>
  summarise(sum_n = sum(n)) |>
  ungroup() |>
  ggplot(aes(
    y = fct_reorder(phylum, sum_n),
    x = sum_n
  )) +
  geom_col() +
  geom_text(aes(y = phylum, x = sum_n / 2, label = sum_n), color = "black") +
  scale_x_log10(guide = guide_axis_logticks(), expand = c(0, 0.3)) +
  labs(
    y = NULL,
    x = "vOTUs"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_line(linetype = 2),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  )

host_plot
# Genome size by host

phylum_factor <- host_counts |>
  inner_join(host_pool_df, by = "phylum") |>
  mutate(phylum = if_else(host_pool, "Other", phylum)) |>
  group_by(phylum) |>
  summarise(sum_n = sum(n))

phylum_factor <- phylum_factor[order(phylum_factor$sum_n), ] |>
  pull(phylum)

size_plot <- host_df |>
  mutate(Virus = str_replace(Virus, "_1$", "")) |>
  inner_join(genome_size, by = c("Virus" = "contig_id")) |>
  select(contig_length, phylum) |>
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = str_remove(phylum, "_[A-Z]"),
    phylum = if_else(is.na(phylum), "Unknown", phylum),
    contig_length_mb = contig_length / 1000
  ) |>
  inner_join(host_pool_df, by = "phylum") |>
  mutate(
    phylum = if_else(host_pool, "Other", phylum),
    phylum = factor(phylum, levels = phylum_factor)
  ) |>
  ggplot(aes(
    y = phylum,
    x = contig_length_mb
  )) +
  geom_boxplot(
    outlier.colour = "grey",
    outlier.shape = 1,
    staplewidth = 0.3,
    width = 0.5
  ) +
  scale_x_log10(guide = guide_axis_logticks()) +
  # scale_x_continuous(
  #   limits = c(0, 450),
  #   breaks = seq(0, 450, 100)
  # ) +
  labs(
    y = NULL,
    x = "Lenght (kb)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_line(linetype = 2),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  )
size_plot

# taxonomy
taxonomy <- read_tsv("data/representative_viral_taxonomy_2025.tsv") %>%
  separate(
    lineage,
    into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
    sep = ";"
  )

tax_df <- host_df |>
  left_join(taxonomy, by = c("Virus" = "seq_name")) |>
  select(phylum.x, class.y) |>
  mutate(
    phylum.x = str_remove(phylum.x, "p__"),
    phylum.x = str_remove(phylum.x, "_[A-Z]"),
    phylum.x = if_else(is.na(phylum.x), "Unknown", phylum.x),
    class.y = if_else(class.y == "", "Unclassified", class.y)
  ) |>
  count(phylum.x, class.y) |>
  inner_join(host_pool_df, by = c("phylum.x" = "phylum")) |>
  mutate(
    phylum.x = if_else(host_pool, "Other", phylum.x),
    phylum.x = factor(phylum.x, levels = phylum_factor)
  ) |>
  group_by(phylum.x) |>
  mutate(rel_abund = 100 * (n / sum(n))) |>
  ungroup()

pool <- tax_df |>
  select(class.y, rel_abund) |>
  group_by(class.y) |>
  summarise(max_prediction = max(rel_abund)) |>
  mutate(pool = if_else(max_prediction > 1, TRUE, FALSE)) |>
  select(class.y, pool)

tax_plot <- tax_df |>
  inner_join(pool, by = "class.y") |>
  mutate(class.y = if_else(pool, class.y, "Other")) |>
  group_by(phylum.x, class.y) |>
  summarise(sum_rel = sum(rel_abund), .groups = "drop") |>
  ggplot(aes(y = phylum.x, x = sum_rel, fill = class.y)) +
  geom_col() +
  scale_fill_manual(values = tax_color) +
  labs(y = NULL, x = "vOTUs (%)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_line(linetype = 2),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5, 5, 1, 5),
    legend.box.margin = margin(0, 0, 0, -9)
  )

tax_plot
# Life style
style_plot <- host_df |>
  full_join(life_style, by = "Virus") |>
  select(phylum, Virulent, Temperate) |>
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = str_remove(phylum, "_[A-Z]"),
    phylum = if_else(is.na(phylum), "Unknown", phylum),
    style = case_when(
      Virulent >= 0.95 ~ "Virulent",
      Temperate >= 0.95 ~ "Temperate",
      Temperate < 0.95 & Temperate < 0.95 ~ "Uncertain"
    ),
    style = factor(style, levels = c("Temperate", "Virulent", "Uncertain"))
  ) |>
  count(phylum, style) |>
  inner_join(host_pool_df, by = "phylum") |>
  mutate(
    phylum = if_else(host_pool, "Other", phylum),
    phylum = factor(phylum, levels = phylum_factor)
  ) |>
  group_by(phylum) |>
  mutate(rel_abund = 100 * (n / sum(n))) |>
  ungroup() |>
  ggplot(aes(y = phylum, x = rel_abund, fill = style)) +
  geom_col() +
  scale_fill_manual(values = style_color) +
  labs(y = NULL, x = "vOTUs (%)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_line(linetype = 2),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(10, "pt"),
    legend.spacing.x = unit(0.1, "pt")
  ) +
  guides(fill = guide_legend(ncol = 2))


# Compose plot

(host_plot | size_plot | style_plot | tax_plot) +
  plot_annotation(tag_levels = 'A')

ggsave(
  filename = "plots/figure3/Figure_3.png",
  width = 12,
  height = 6,
  dpi = 400
)


# Host by Genus
host_counts_genus <- host_df |>
  select(genus) |>
  mutate(
    genus = str_remove(genus, "g__"),
    genus = str_remove(genus, "_[A-Z]"),
    genus = if_else(is.na(genus), "Unknown", genus)
  ) |>
  count(genus)


host_pool_df_genus <- host_counts_genus |>
  mutate(host_pool = if_else(n < 100, TRUE, FALSE)) |>
  select(genus, host_pool)

host_plot_genus <- host_counts_genus |>
  inner_join(host_pool_df_genus, by = "genus") |>
  mutate(genus = if_else(host_pool, "Other", genus)) |>
  group_by(genus) |>
  summarise(sum_n = sum(n)) |>
  ungroup() |>
  ggplot(aes(
    y = fct_reorder(genus, sum_n),
    x = sum_n
  )) +
  geom_col() +
  geom_text(aes(y = genus, x = sum_n / 2, label = sum_n), color = "black") +
  scale_x_log10() +
  labs(
    y = NULL,
    x = "vOTUs"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_line(linetype = 2),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  )


ggsave(
  host_plot_genus,
  file = "plots/host_genus.png",
  width = 5,
  height = 8,
  units = "in",
  dpi = 300
)


# Polyvalent phages
polyvalent <- data_join |>
  right_join(votu, by = c("Virus")) |>
  filter(grepl("CRISPR", `List of methods`)) |>
  select(Virus, host_tax) |>
  separate(
    host_tax,
    into = c("domain", "phylum", "class", "order", "family", "genus"),
    sep = ";"
  ) %>%
  select(Virus, genus) %>%
  unique() %>%
  count(Virus)

polyvalent <- data_join |>
  right_join(votu, by = c("Virus")) |>
  filter(grepl("CRISPR", `List of methods`)) |>
  select(Virus, host_tax) |>
  separate(
    host_tax,
    into = c("domain", "phylum", "class", "order", "family", "genus"),
    sep = ";"
  ) %>%
  select(Virus, order) %>%
  unique() %>%
  count(Virus)


polyvalent_plot <- data.frame(
  tax_level = c("Genus", "Family", "Order"),
  value = c(5.601, 0.353, 0.055)
) %>%
  mutate(
    tax_level = factor(tax_level, levels = c("Genus", "Family", "Order"))
  ) %>%
  ggplot(aes(x = tax_level, y = value)) +
  geom_col(width = 0.6) +
  geom_text(
    aes(x = tax_level, y = value + 0.2, label = round(value, 2)),
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, 0.5),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Number of phages (%)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

ggsave(
  polyvalent_plot,
  file = "plots/polyvalent.png",
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 300
)

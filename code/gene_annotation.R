# Load libraries
library(tidyverse)
library(patchwork)

# Load data and set colours
annotation_df <- read_tsv("data/pharokka_cds_final_merged_output.tsv")

base_color <- c("black", "grey")


# Pie plot total genes ---------------------------------------------------------------------
pie_df <- annotation_df %>%
  count(category) %>%
  mutate(
    new_cat = if_else(
      category == "unknown function",
      "Unknown function",
      "Known function"
    )
  ) %>%
  group_by(new_cat) %>%
  summarise(sum_n = sum(n), .groups = "drop") %>%
  mutate(percent = 100 * (sum_n / sum(sum_n))) %>%
  mutate(ymax = cumsum(percent), ymin = c(0, head(ymax, n = -1)))

# Position of labels
pie_df$labelPosition <- (pie_df$ymax + pie_df$ymin) / 2
pie_df$label <- paste0(pie_df$new_cat, "\n(", round(pie_df$percent, 1), "%)")

pie_plot <- pie_df %>%
  ggplot(aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 3.5,
    xmin = 2.5,
    fill = new_cat
  )) +
  geom_rect(color = "white") +
  geom_text(
    x = 4.3,
    aes(y = labelPosition, label = label),
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    x = 1.5,
    y = 0.2,
    label = "Genes\n999.462",
    size = 6,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = base_color) +
  xlim(c(1.5, 4)) +
  theme_void() +
  theme(legend.position = "none")

pie_plot

ggsave(pie_plot, filename = "plots/figure4/pie_genes.pdf", dpi = 300)

# Pie plot category ------------------------------------------------------------------------------------------
cat_colors <- c(
  'black',
  '#ffffd9',
  '#c7e9b4',
  '#525252',
  '#7fcdbb',
  '#d9d9d9',
  '#1d91c0',
  '#225ea8',
  '#0c2c84'
)

pie_cat <- annotation_df %>%
  count(category) %>%
  filter(!category == "unknown function") %>%
  group_by(category) %>%
  summarise(sum_n = sum(n), .groups = "drop") %>%
  mutate(percent = 100 * (sum_n / sum(sum_n))) %>%
  mutate(
    ymax = cumsum(percent),
    ymin = c(0, head(ymax, n = -1)),
    category = case_when(
      category == "connector" ~ "Connector",
      category == "DNA, RNA and nucleotide metabolism" ~
        "Nucleotide metabolism",
      category == "head and packaging" ~ "Head/Packaging",
      category == "integration and excision" ~ "Integration/Excision",
      category == "lysis" ~ "Lysis",
      category == "moron, auxiliary metabolic gene and host takeover" ~
        "Auxiliary metabolic",
      category == "other" ~ "Other",
      category == "tail" ~ "Tail",
      category == "transcription regulation" ~ "Transcription"
    )
  )

pie_cat$labelPosition <- (pie_cat$ymax + pie_cat$ymin) / 2
pie_cat$label <- paste0(
  pie_cat$category,
  "\n(",
  round(pie_cat$percent, 1),
  "%)"
)

pie_cat_plot <- pie_cat %>%
  ggplot(aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 3.5,
    xmin = 2.5,
    fill = category
  )) +
  geom_rect(color = "white") +
  geom_text(
    x = 3.8,
    aes(y = labelPosition, label = label),
    size = 3,
    fontface = "bold"
  ) +
  geom_text(
    x = 1.5,
    y = 0.2,
    label = "Genes\n187.886",
    size = 6,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cat_colors) +
  xlim(c(1.5, 4)) +
  theme_void() +
  theme(legend.position = "none")


ggsave(pie_cat_plot, filename = "plots/figure4/pie_cat.pdf", dpi = 300)

# Bar plot functions ------------------------------------------------------------
cat_colors_bars <- c(
  '#ffffd9',
  '#d9d9d9',
  '#c7e9b4',
  '#525252',
  '#7fcdbb',
  '#1d91c0',
  '#225ea8',
  '#0c2c84'
)

barplot <- annotation_df %>%
  count(annot, category) %>%
  mutate(
    percent = 100 * (n / sum(n)),
    annot = str_to_sentence(annot),
    annot = str_replace(annot, "Dna", "DNA"),
    annot = str_replace(annot, "Rna", "RNA"),
    annot = str_replace(annot, "rna", "RNA"),
    category = str_to_sentence(category)
  ) %>%
  #filter(category != "Unknown function") %>%
  filter(percent > 0.1, annot != "Hypothetical protein") %>%
  ggplot(aes(x = percent, y = fct_reorder(annot, n), fill = category)) +
  geom_col() +
  geom_text(
    aes(x = percent + 0.02, y = annot, label = annot),
    hjust = 0,
    size = 3.7
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 1.6),
    breaks = seq(0, 1.6, 0.2)
  ) +
  scale_fill_manual(
    values = cat_colors_bars,
    labels = c(
      "Dna, rna and nucleotide metabolism" = "Nucleotide metabolism",
      "Head and packaging" = "Head/Packaging",
      "Integration and excision" = "Integration/Excision"
    ),
  ) +
  labs(y = NULL, x = "Number of genes (%)") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    #legend.box = "horizontal",
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

barplot

ggsave(
  barplot,
  filename = "plots/figure4/barplot_genes.pdf",
  width = 3,
  height = 8,
  dpi = 300
)
# Defense plot -----------------------------------------------------------------

defense_genus <- defense_df %>%
  select(sys_beg, activity) %>%
  mutate(sys_beg = str_remove(sys_beg, "_[^_]+$")) %>%
  #count(sys_beg)
  left_join(host_df, by = c("sys_beg" = "Virus")) %>%
  count(genus)


defense_df <- read_tsv("data/representative_viral_defense_finder_systems.tsv")

defense_trans <- defense_df %>%
  select(sys_beg, activity) %>%
  mutate(sys_beg = str_remove(sys_beg, "_[0-9]+$")) %>%
  count(sys_beg, activity) %>%
  count(activity) %>%
  mutate(percent = n / 19778)

color_code <- c(
  '#9e0142',
  '#d53e4f',
  '#f46d43',
  '#fdae61',
  '#fee08b',
  'black',
  '#e6f598',
  '#abdda4',
  '#66c2a5',
  '#3288bd',
  '#5e4fa2',
  '#8dd3c7',
  '#ffffb3',
  '#bebada',
  '#fb8072',
  '#80b1d3',
  '#fdb462',
  '#b3de69',
  '#fccde5',
  '#d9d9d9',
  '#bc80bd',
  '#ccebc5',
  '#ffed6f'
)


defense_plot <- defense_df %>%
  select(activity, type) %>%
  mutate(
    type = str_replace(type, "Abi.*", "Abi"),
    type = str_replace(type, "PD.*", "PD"),
    type = str_remove(type, "Anti_")
  ) %>%
  count(activity, type) %>%
  group_by(activity) %>%
  mutate(
    percent = 100 * (n / sum(n)),
    group = if_else(percent < 1, "Other", type)
  ) %>%
  ungroup() %>%
  group_by(group, activity) %>%
  summarise(sum_percent = sum(percent), .groups = "drop") %>%
  ggplot(aes(x = activity, y = sum_percent, fill = group)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
  scale_fill_manual(values = color_code) +
  labs(x = NULL, y = "Number of system types (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.2),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

ggsave(
  defense_plot,
  filename = "plots/figure4/defense.pdf",
  width = 6,
  height = 5,
  dpi = 300
)

# Compose plot -----------------------------------------------------------------

fig_4 <- (((pie_plot | pie_cat_plot) / defense_plot) | barplot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18))


ggsave(
  fig_4,
  filename = "plots/figure4/Figure_4.png",
  width = 10,
  height = 10,
  dpi = 300
)


pies <- wrap_elements(pie_plot + pie_cat_plot)

wrap_plots((pies / defense_plot) | barplot, ncol = 1)

################

qannotation_df %>%
  count(category) %>%
  mutate(percent = 100 * (n / sum(n)), category = str_to_sentence(category)) %>%
  #filter(category != "Unknown function") %>%
  ggplot(aes(x = percent, y = fct_reorder(category, n))) +
  geom_col() +
  geom_text(
    aes(x = percent + 4, y = category, label = paste0(round(percent, 1), "%")),
    size = 2.5
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 90),
    breaks = seq(0, 90, 10)
  ) +
  scale_y_discrete(
    labels = c(
      "Dna, rna and nucleotide metabolism" = "DNA, RNA and\n nucleotide metabolism",
      "Moron, auxiliary metabolic gene and host takeover" = "Auxiliary metabolic gene\nand host takeover"
    )
  ) +
  labs(y = NULL, x = "Number of genes (%)") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  )


ggsave(filename = "plots/gene_categories.png", width = 5, height = 5, dpi = 400)

test <- annotation_df %>%
  filter(
    category == "moron, auxiliary metabolic gene and host takeover",
    annot == "phosphoadenosine phosphosulfate reductase"
  ) %>%
  count(contig)

test2 <- test %>%
  left_join(host_df, by = c("contig" = "Virus")) %>%
  count(genus)


ggsave(
  file = "plots/gene_annotation.png",
  width = 4,
  height = 7,
  units = "in",
  dpi = 300
)


annotation_df %>%
  filter(category == "moron, auxiliary metabolic gene and host takeover") %>%
  count(annot, category) %>%
  mutate(percent = 100 * (n / sum(n)), annot = str_to_sentence(annot)) %>%
  filter(percent > 1) %>%
  ggplot(aes(x = percent, y = fct_reorder(annot, n))) +
  geom_col() +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)
  ) +
  labs(y = NULL, x = "Number of genes (%)") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    #legend.box = "horizontal",
    panel.grid.major.x = element_line(linetype = 2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

ggsave(
  file = "plots/gene_annotation_amg.png",
  width = 4,
  height = 6,
  units = "in",
  dpi = 300
)

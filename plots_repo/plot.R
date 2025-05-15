# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)

# Load metadata ----------------------------------------------------------------
df <- read_csv("data_availability/metagenomic_samples_metadata.csv")
color_continent <- c('#a6611a', '#dfc27d', '#80cdc1', '#018571')
# Wrangle metadata -------------------------------------------------------------

country <- df %>%
  filter(geo_loc_name_country != "United Kingdom") %>% # We did not use UK samples for collection
  select(geo_loc_name_country, geo_loc_name_country_continent) %>%
  count(geo_loc_name_country, geo_loc_name_country_continent) %>%
  ggplot(aes(
    x = n,
    y = fct_reorder(geo_loc_name_country, n),
    fill = geo_loc_name_country_continent
  )) +
  geom_col() +
  geom_text(
    aes(x = n + 5, y = geo_loc_name_country, label = n),
    hjust = 0,
    size = 3
  ) +
  scale_fill_manual(values = color_continent) +
  scale_x_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  labs(y = NULL, x = "Number of samples") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(linetype = 2, linewidth = 0.2),
    axis.ticks.y = element_blank(),
    legend.position = c(0.7, 0.3)
  )

gut_region <- df %>%
  filter(geo_loc_name_country != "United Kingdom") %>% # We did not use UK samples for collection
  select(isolation_source) %>%
  count(isolation_source) %>%
  mutate(isolation_source = str_to_sentence(isolation_source)) %>%
  ggplot(aes(x = n, y = fct_reorder(isolation_source, n))) +
  geom_col() +
  geom_text(
    aes(x = n + 5, y = isolation_source, label = n),
    hjust = 0,
    size = 3
  ) +
  scale_x_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  labs(y = NULL, x = "Number of samples") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(linetype = 2, linewidth = 0.2),
    axis.ticks.y = element_blank()
  )


country + gut_region

ggsave(
  filename = "plots_repo/summary_samples.png",
  dpi = 300,
  width = 7,
  height = 3
)

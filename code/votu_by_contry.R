# Determine the origin of the viral genomes

library(tidyverse)

# load files
df <- read_tsv("output/representatives_id.txt")

metadata <- read_csv("data/MAGs_dataset_complete.csv") %>% 
  select(Run, geo_loc_name_country, isolation_source)

# Wrangling data --------------------------------------------------------------

virus_by_contry <- df %>% 
  mutate(Run = str_remove(cluster_representative, "_k141_.*|_k119_.*|_k99_.*")
         ) %>% 
  left_join(metadata, by= "Run") %>% 
  count(geo_loc_name_country) %>% 
  mutate(percent = n/sum(n) *100)

virus_by_contry %>% 
  ggplot(aes(y = fct_reorder(geo_loc_name_country, n),
             x = n)) +
  geom_col() +
  geom_text(aes(x = n + 300,
            y = geo_loc_name_country,
            label = n),
            size = 3) +
  scale_x_continuous(breaks = seq(0, 9000, 1000),
                     ) +
  labs(y = NULL,
       x = "vOTUs") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 2),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "white",
                                        fill = "white"),
        plot.background = element_rect(colour = "white",
                                        fill = "white"))


ggsave(filename = "plots/suplementary/virus_country.png",
       width = 6,
       height = 4,
       dpi = 300)

# vOTUs by region --------------------------------------------------------------

virus_by_region <- df %>% 
  mutate(Run = str_remove(cluster_representative, "_k141_.*|_k119_.*|_k99_.*")
  ) %>% 
  left_join(metadata, by= "Run") %>% 
  count(isolation_source) %>% 
  mutate(percent = n/sum(n) *100,
         isolation_source = str_to_sentence(isolation_source))

virus_by_region %>% 
  ggplot(aes(y = fct_reorder(isolation_source, n),
             x = n)) +
  geom_col() +
  geom_text(aes(x = n + 500,
                y = isolation_source,
                label = n),
            size = 3) +
  scale_x_continuous(breaks = seq(0, 16000, 2000),
  ) +
  labs(y = NULL,
       x = "vOTUs") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 2),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "white",
                                        fill = "white"),
        plot.background = element_rect(colour = "white",
                                       fill = "white"))

ggsave(filename = "plots/suplementary/virus_region.png",
       width = 6,
       height = 4,
       dpi = 300)

library(tidyverse)

metadata <- read_csv("data/MAGs_dataset_complete.csv")

# Find number of samples per continent
metadata %>% 
  select(geo_loc_name_country_continent) %>% 
  count(geo_loc_name_country_continent)


metadata %>% 
  select(Run, geo_loc_name_country, Bytes, Bases) %>%
  filter(geo_loc_name_country != "United Kingdom") %>% 
  drop_na() %>% 
  summarise(sum_bytes = sum(Bytes),
            sum_bases = sum(Bases))


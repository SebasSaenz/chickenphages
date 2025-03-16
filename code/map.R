library(tidyverse)
library(patchwork)
library(here)
library(dplyr)
library(stringr)
library(ggplot2)
library(maps)
library(ggtext)
library(ggrepel)

options(scipen = 999) ## To disable scientific notation

world <- map_data("world")
head(world)

metadata <- read_csv("data/MAGs_dataset_complete.csv")

df <- metadata |> 
  rename(region = "geographic_location_(country_and/or_sea)") |>  
  select(region) |> 
  count(region)

labels <- inner_join(df, world, by = "region") %>% 
  group_by(region) %>% 
  summarise(avr_long = mean(long),
         avr_lat = mean(lat))


color_map <- c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494','lightgrey')

world %>%
  filter(region != "Antarctica") %>% 
  left_join(., df, by = "region") %>%
  mutate(n = case_when(n <20 ~ "1-19",
                       n >= 20 & n < 40 ~ "20-39",
                       n >= 40 & n < 50 ~ "40-49",
                       n >= 50 & n <= 60 ~ "50-60",
                       n >= 500 & n < 600 ~ "500-599",
                       n >= 600 & n < 700 ~ "600-700",
                       is.na(n) ~ "No data")) %>% 
  ggplot(aes(x=long, y = lat, group = group)) +
  #coord_fixed(1.2) +
  geom_polygon(aes(fill = n), color = "white") +
  geom_text_repel(data = labels, aes(x = avr_long,
                               y = avr_lat,
                               label = region),
                  min.segment.length = 1,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  size = 3,
                  inherit.aes = FALSE) +
  scale_fill_manual(values = color_map) +
  theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.4),
    legend.position = c(0.5,0.5),
    legend.direction= "horizontal",
    legend.title = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")) +
  guides(fill = guide_legend(nrow = 1))

worldplot

ggsave(filename = "plots/map.png", width = 12, height = 6, dpi = 400)


git <- metadata |> 
  mutate(study = case_when(grepl("UHO", Run)~"This study", 
                           grepl("LH1", Run)~"This study",
                           .default = "Other study")) |>
  count(isolation_source, study) |> 
  mutate(isolation_source = str_to_title(isolation_source)) |> 
  ggplot(aes(y = factor(isolation_source, 
                        levels = c("Crop", "Colorectum", "Duodenum", "Jejunum",
                                   "Caeca","Faeces","Ileum")),
             x = n,
             fill = study)) +
  geom_col()+
  scale_x_continuous(limits = c(0, 600),
                     breaks = seq(0 ,600, 100))+
  scale_fill_manual(values = c("black", "grey")) +
  labs(y = "Sample source",
         x = "Number of samples") +
  theme_bw() +
  theme(legend.title= element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size = 12, , family = "Helvetica"))

ggsave(file= "plots/samnple_source.png", width = 5, height = 4, dpi = 450)


sex <- metadata |>
  mutate(study = case_when(grepl("UHO", Run)~"This study", 
                           grepl("LH1", Run)~"This study",
                           .default = "Other study")) |> 
  count(host_gender, study) |> 
  mutate(host_gender = str_to_title(host_gender),
         host_gender = if_else(is.na(host_gender), "Unknown", host_gender)) |> 
  ggplot(aes(y = fct_reorder(host_gender, n),
             x = n,
             fill = study)) +
  geom_col(width = 0.5) +
  scale_x_continuous(limits = c(0, 700),
                     breaks = seq(0 ,700, 100)) +
  scale_fill_manual(values = c("black", "grey")) +
  labs(x = "Number of samples",
       y = "Sex")+
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(linetype = 2),
        text = element_text(size = 12, family = "Helvetica"))

git + sex +
  plot_layout(guides = 'collect') & theme(legend.position = "bottom")

ggsave(filename = "plots/sex_samples.png", width = 6.5, height = 3, dpi = 300)





x <- metadata |> 
  filter(geo_loc_name_country == "Germany",
         grepl("UHO", Run)) |> 
  count(Host_age)
  
metadata |> 
  count(geo_loc_name_country)

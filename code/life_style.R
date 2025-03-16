library(tidyverse)
library(ggrepel)

base_color <- c('#01665e', '#762a83', '#2166ac','#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf')

life_style <- read_tsv("data/representative_viral.fasta_bacphlip.txt")

style_df <- life_style |>
  rename(virus =...1) |> 
  mutate(style = case_when(Virulent >= 0.95 ~ "Virulent",
                           Temperate >= 0.95 ~ "Temperate",
                           Temperate < 0.95 & Temperate < 0.95  ~ "Uncertain"))


style_df |> 
  count(style) |> 
  mutate(percentage = 100 * (n/sum(n)),
         ymax = cumsum(percentage),
         ymin = c(0, head(ymax, n=-1))) |> 
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=style)) +
  geom_rect(color = "white") +
  coord_polar(theta="y") + 
  scale_fill_manual(values = base_color) +
  xlim(c(2, 4)) +
  theme_void()

ggsave(filename = "plots/figure2/lifestyle.pdf", width = 5, height = 5, dpi = 400)

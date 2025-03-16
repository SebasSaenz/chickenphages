# Calculate the number of defence systems in the vOTUs
# Figure number 4

library(tidyverse)

df <- read_tsv("data/representative_viral_defense_finder_systems.tsv") 

x <- df %>% select(sys_beg, activity) %>% 
  mutate(sys_beg = str_remove(sys_beg, "_[0-9]+$")) %>% 
  count(sys_beg, activity) %>% 
  count(activity) %>% 
  mutate(percent = n /19778)

color_code <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','black' ,'#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f' )


df %>% select(activity, type) %>% 
  mutate(type = str_replace(type, "Abi.*", "Abi"),
         type = str_replace(type, "PD.*", "PD"),
         type = str_remove(type, "Anti_")) %>% 
  count(activity, type) %>%
  mutate(group = if_else(n < 5, "Other", type)) %>% 
  ggplot(aes(x = n,
             y = group,
             fill = activity)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 800, 100),
                     limits = c(0, 820),
                     expand = c(0, 0)) +
    scale_fill_manual(values = c("black", "grey")) +
  labs(y = NULL,
       x = "Number of systems") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = 2, colour = "grey", linewidth = 0.2),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.9),
        text = element_text(size=14))



ggsave(filename = "plots/figure4/defense.pdf", width = 4, height = 5, dpi = 350)


df %>% select(activity, type) %>% 
  mutate(type = str_replace(type, "Abi.*", "Abi"),
         type = str_replace(type, "PD.*", "PD"),
         type = str_remove(type, "Anti_")) %>% 
  count(activity, type) %>%
  group_by(activity) %>% 
  mutate(percent = 100 * (n/sum(n)),
         group = if_else(percent < 1, "Other", type)) %>%
  ungroup() %>% 
  ggplot(aes(x = activity,
         y = percent,
         fill = group)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0, 100, 10),
                     expand = c(0, 0)) +
  scale_fill_manual(values = color_code) +
  labs(x = NULL,
       y = "Number of system types (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(linetype = 2, linewidth  = 0.2),
        legend.title = element_blank(),
        axis.ticks.x = element_blank())
  

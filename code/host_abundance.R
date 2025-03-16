library(tidyverse)


filenames <- list.files(path = "data/coverage_host/logs/",
                        pattern = "*sh*",
                        full.names = TRUE)


pattern_net <- "(?<=found\\s)\\w+"
pattern_pattern2 <- "(?<=out of\\s)\\w+"
db <- "chicken"
# Apply the function to each file and combine the results

results_df <- map_dfr(filenames, ~ extract_words(.x, pattern_net, pattern_pattern2, db)) %>% 
  mutate(Perc_mapped = 100*(mapped/total_reads)) %>% 
  rename(sample=file, Host=data_base, Mapped=mapped, Total=total_reads)

host <- read_tsv("data/coverage_host/CoverM_mapped_host.tsv") %>% 
  rename("sample"="...1")



rbind(host, results_df) %>%
  mutate(Host = factor(Host,
                       levels = rev(c("turkey", "chicken", "duck", "human", "pig", "cow")),
                       labels = rev(c("Turkey (n=69)", "Chicken (n=50)", "Duck (n=41)",
                                      "Human (n=44)", "Pig (n=48)", "Cow (n=37)")))) %>% 
  ggplot(aes(y = Host,
             x = Perc_mapped)) +
  geom_jitter(shape = 21, color = "black", fill = "grey", size = 0.5, height = 0.2) +
  geom_boxplot(outlier.size = 1, width = 0.5, alpha = 0) +
  scale_x_continuous(limit = c(0 ,7),
                     breaks = seq(0, 7 ,1)) +
  labs(x = "Reads mapped (%)",
       y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = 2, linewidth = 0.2, color = "grey"),
        axis.text.y = element_text(size = 10))



ggsave(filename = "plots/fig5/host.png", width = 4, height = 3, dpi = 300)



df <- read_tsv("data/coverage_host/CoverM_RPKM_host.tsv")

metadata <- read_tsv("data/coverage_host/CoverM_mapped_host.tsv") %>% 
  rename("sample"="...1") %>% 
  select(sample, Host)
  



df %>% 
  pivot_longer(-Contig, names_to = "sample", values_to = "value") %>% 
  inner_join(metadata, by ="sample") %>% 
  mutate(presence = if_else(value > 1, 1, 0)) %>% 
  ggplot(aes(x = sample,
             y = Contig,
             fill = presence)) +
  geom_tile() +
  scale_fill_gradient(high="black", low="white") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())


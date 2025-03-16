pacman::p_load(
  tidyverse,
  data.table,
  ggridges)


all_paths <- list.files(
  path = "data/genomad_files/",
  pattern = "*.tsv", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8"
  )

all_filenames <- all_paths %>% # Manipulate file paths
  basename() %>%
  as.list()

all_lists <- mapply(c,
                    all_content,
                    all_filenames,
                    SIMPLIFY = FALSE
)

all_result <- rbindlist(all_lists, fill = T)
names(all_result)[12] <- "sample" # change column name

# COunt viral contigs per samples -------
 x <- all_result %>% 
  select(length, taxonomy) %>% 
  separate(taxonomy,
           into = c("root", "realm", "kingdom", "phylum", "class", "order", "family"),
           sep = ";") %>% 
  count(class, order) %>% 
  filter(grepl("ae", order))

x %>% 
  ggplot(aes(x = n)) + 
  geom_density_ridges() +
  scale_x_continuous(breaks = seq(0, 1300, 100))



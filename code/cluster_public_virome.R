# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(VennDiagram)
library(grid)         # for grid.newpage()

# 1. Read & parse the TSV
df <- read_tsv("data/my_clusters_CGVD.tsv",
               col_names = c("rep","members"),
               col_types = cols(
                 rep     = col_character(),
                 members = col_character()
               )) %>%
  mutate(cluster = row_number(),
         # strip quotes & split the comma-list into a character vector
         members = str_remove_all(members, '"') %>% str_split(","))

# 2. Unnest reps + members into one long table
df_long <- df %>%
  # rep as one row per cluster…
  transmute(cluster, genome = rep) %>%
  bind_rows(
    # … plus every member as its own row
    df %>% select(cluster, genome = members) %>% unnest(genome)
  ) %>%
  # 3. flag which are “Germany” vs. “Other”
  mutate(is_germany = str_detect(genome, regex("germany", ignore_case = TRUE)))

# 4. Compute sets
all_germany     <- df_long %>% filter(is_germany)     %>% distinct(genome) %>% pull(genome)
all_non_germany <- df_long %>% filter(!is_germany)    %>% distinct(genome) %>% pull(genome)

shared_clusters <- df_long %>%
  group_by(cluster) %>%
  summarize(has_germany = any(is_germany),
            has_other   = any(!is_germany)) %>%
  filter(has_germany & has_other) %>%
  pull(cluster)

shared_germany     <- df_long %>% filter(cluster %in% shared_clusters, is_germany)  %>% distinct(genome) %>% pull(genome)
shared_non_germany <- df_long %>% filter(cluster %in% shared_clusters, !is_germany) %>% distinct(genome) %>% pull(genome)

unique_germany     <- setdiff(all_germany,     shared_germany)
unique_non_germany <- setdiff(all_non_germany, shared_non_germany)

# 5. Report counts
tibble(
  category = c("Germany‐only", "Other‐only", "Shared"),
  count    = c(length(unique_germany),
               length(unique_non_germany),
               length(shared_germany))
) %>% mutate(percent = 100 * count/sum(count)) %>%  print()

# 6. Venn diagram (2-set)
grid.newpage()
draw.pairwise.venn(
  area1      = length(all_germany),
  area2      = length(all_non_germany),
  cross.area = length(shared_germany),
  category   = c("Germany Genomes", "Other Genomes"),
  cat.pos    = c(-20, 20),
  euler.d    = FALSE,
  scaled     = FALSE
)

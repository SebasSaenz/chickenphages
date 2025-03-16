library(tidyverse)

df<- read_csv("final_assignments.csv")

x <- df |> 
  filter(is.na(RefSeqID)) |> 
  select(`family (prediction)`) |> 
  mutate(`family (prediction)` = if_else(grepl("novel", `family (prediction)`), "Novel", `family (prediction)`)) |> 
  count(`family (prediction)`) 

         
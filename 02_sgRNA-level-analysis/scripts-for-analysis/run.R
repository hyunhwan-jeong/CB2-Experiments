library(CB2)
library(readxl)
library(tidyverse)

read_xlsx("Dolcetto.xlsx", sheet = 1, skip = 2) %>% print

df_crispri_a375_count <- read_xlsx("Dolcetto.xlsx", sheet = 1, skip = 2) %>%
  select(-3,-4,-5) %>% select(sgRNA=1, pDNA=2, RepA=3, RepB=4, RepC=5)
df_crispri_a375_info <- read_xlsx("Dolcetto.xlsx", sheet = 3) %>%
  select(sgRNA=1, Gene=2)

print(df_crispri_a375_count)
print(df_crispri_a375_info)

left_join(df_crispri_a375_count, df_crispri_a375_info, by = "sgRNA")  %>% 
  select(sgRNA, Gene, everything()) %>% 
  write_delim("input.txt", delim="\t")


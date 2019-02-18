library(here)
library(tidyverse)

msg <- read_tsv("RRA_tmp/sample1.sgrna_summary.txt")
msg$control_mean %>% summary
msg$treat_mean %>% summary

mg <- read_tsv("RRA_tmp/sample1.gene_summary.txt")
mg$`neg|goodsgrna` %>% sum
nrow(msg)

library(CB2)

read_tsv("RRA_tmp/input.tsv") %>% 
  unite(id, c("gene", "sgRNA")) %>% column_to_rownames("id") -> df_count

df_design <- tribble(
  ~sample_name, ~group,
  "B1", "before",
  "B2", "before",
  "B3", "before",
  "A1", "after",
  "A2", "after",
  "A3", "after"
)

sgrna <- run_estimation(df_count, df_design, "before", "after")

sgrna %>% 
  mutate(pool = "list") %>%
  mutate(prob = 1) %>% 
  mutate(chosen = 1 ) %>% 
  mutate(plow = t_value %>% unlist) %>% 
  select(sgrna = sgRNA,
         symbol = gene,
         pool, 
         plow,
         prob, chosen) %>% arrange(plow[,1]) %>% write_tsv("RRA_tmp/cb_rra.txt")

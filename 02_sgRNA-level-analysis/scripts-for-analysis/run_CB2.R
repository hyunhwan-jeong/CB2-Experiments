#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(here)

df_count <- read_tsv("input.txt")

df_count %>% unite("id", c("Gene", "sgRNA")) %>% column_to_rownames("id") -> df_count

ctl <- c("pDNA")
trt <- c("RepA", "RepB", "RepC")

df_design <-
    rbind(
        data.frame(group = "ctl", sample_name = ctl, stringsAsFactors = FALSE),
        data.frame(group = "trt", sample_name = trt, stringsAsFactors = FALSE)
    )

print(df_design)
sgRNA_stat <- CB2::run_estimation(df_count, df_design, "ctl", "trt")
print(sgRNA_stat %>% select(sgRNA, t_value, p_pa, p_pb))
write_csv(sgRNA_stat, "output_cb2.csv")

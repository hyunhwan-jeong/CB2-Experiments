#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(here)

source(here("utils", "parse_args.R"))

df_count %>% unite("id", c("gene", "sgRNA")) %>% column_to_rownames("id") -> df_count

df_design <-
    rbind(
        data.frame(group = "trt", sample_name = trt),
        data.frame(group = "ctl", sample_name = ctl)
    )

sgRNA_stat <- CB2::run_estimation(df_count, df_design, "ctl", "trt")
gene_stat <- CB2::measure_gene_stats(sgRNA_stat)

gene_stat %>% select(gene = gene, fdr = fdr_pb, stat = p_pb) %>% write_csv(outfile)

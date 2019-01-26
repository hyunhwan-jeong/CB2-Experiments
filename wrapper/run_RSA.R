#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(here)
library(glue)

source(here("utils", "parse_args.R"))


df_count %>% unite("id", c("gene", "sgRNA")) %>% column_to_rownames("id") -> df_count

df_design <-
    rbind(
        data.frame(group = "trt", sample_name = trt),
        data.frame(group = "ctl", sample_name = ctl)
    )

sgrna_stat <- CB2::run_estimation(df_count, df_design, "ctl", "trt")


PY3_PATH <- "/Users/hyunhwan/miniconda3/bin/python3"
RSA_PATH <- "/Users/hyunhwan/Projects/InProgress/RSA"


TMP_INPUT <- tempfile()
TMP_OUTPUT <- tempfile()


TMP_INPUT <- tempfile()
TMP_OUTPUT <- tempfile()

sgrna_stat %>%
  mutate(FC = 2^logFC) %>%
  select(gene, sgRNA, FC) -> df_RSA

df_RSA %>% group_by(gene) %>% summarise(n = n()) %>% filter(n > 1) %>% ungroup() %>% pull(gene) -> keeping_genes
df_RSA %>%
  #filter(gene %in% keeping_genes) %>%
  write_csv(TMP_INPUT)

"{PY3_PATH} {RSA_PATH}/RSA.py {TMP_INPUT} --lb 0 --ub 1e8 -o {TMP_OUTPUT} -g gene -s FC -b" %>% glue %>% system


read_csv(TMP_OUTPUT) %>%
  group_by(gene) %>%
  summarise(logFDR = mean(Log_q),
            Hit = mean(RSA_Hit),
            logP = mean(LogP)) %>%
  mutate(fdr = 10^logFDR, stat = 10^logP) -> ret_RSA

print(ret_RSA)

ret_RSA %>% select(gene, fdr, stat) %>% write_csv(outfile)

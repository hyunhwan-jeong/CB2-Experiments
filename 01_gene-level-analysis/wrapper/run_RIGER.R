#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(here)
library(glue)

source(here("utils", "parse_args.R"))


df_count %>% unite("id", c("gene", "sgRNA")) %>% column_to_rownames("id") -> df_count

df_design <-
    rbind(
        data.frame(group = "ctl", sample_name = ctl),
        data.frame(group = "trt", sample_name = trt)
    )

sgrna_stat <- CB2::run_estimation(df_count, df_design, "ctl", "trt")

RIGER_PATH <- "/Users/hyunhwan/Projects/InProgress/rigerj"
RIGER_VERSION <- "2.0.2"


TMP_INPUT <- tempfile()
TMP_OUTPUT <- tempfile()

sgrna_stat %>%
  select(Construct = sgRNA, GeneSymbol = gene, NormalizedScore = logFC) %>%
  mutate(NormalizedScore = NormalizedScore) %>%
  mutate("Construct Rank" = rank(NormalizedScore, ties.method = "first")) -> df_RIGER

df_RIGER %>% group_by(GeneSymbol) %>% summarise(n = n()) %>% filter(n > 1) %>% ungroup() %>% pull(GeneSymbol) -> keeping_genes

df_RIGER %>% filter(GeneSymbol %in% keeping_genes) %>%
  write_delim(TMP_INPUT, delim="\t")

"java -jar {RIGER_PATH}/target/rigerj-{RIGER_VERSION}-assembly.jar -inputFile {TMP_INPUT} -outputFile {TMP_OUTPUT} -numRandomScoresPerGeneSetSize 100000" %>% glue %>% system

read_delim(TMP_OUTPUT, delim="\t") %>%
  dplyr::rename(gene = "Gene Name") %>%
  dplyr::select(-1) %>%
  mutate(`p-value` = ifelse(`p-value` < 0, 0, `p-value`)) %>%
  mutate(fdr = p.adjust(`p-value`, method="fdr")) -> ret_RIGER

print(ret_RIGER)
ret_RIGER %>% select(gene, fdr, stat = `p-value`) %>% write_csv(outfile)

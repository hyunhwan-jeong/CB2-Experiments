#!/usr/bin/env Rscript
library(tidyverse)
library(glue)
library(here)
source(here("utils", "parse_args.R"))


trt <- paste0(trt, collapse = ",")
ctl <- paste0(ctl, collapse = ",")

tmp_input <- tempfile()
tmp_output <- tempfile()

df_count %>% select(sgRNA = sgRNA, gene=gene, everything()) %>% write_delim(tmp_input, delim="\t")
print(params)


"mageck test -k {tmp_input} -t {trt} -c {ctl} -n {tmp_output} {params}" %>% glue -> cmd 
system(cmd)
read_delim("{tmp_output}.gene_summary.txt" %>% glue, delim="\t") %>% 
    select(gene = id, fdr = `neg|fdr`, stat = `neg|score`) %>% write_csv(outfile)


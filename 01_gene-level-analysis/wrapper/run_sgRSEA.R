#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(PBNPA)
library(here)

source(here("utils", "parse_args.R"))
df_count <- as.data.frame(df_count)

multiplier <- 100
eval(parse(text=params))

library(tidyverse)
library(sgRSEA)

if(length(ctl)==1) {
    ctl <- rep(ctl, length(trt))
}

dat <- df_count %>% select(sgRNA = sgRNA, gene =gene, everything()) %>% as.data.frame %>% 
    UQnormalize(trt=trt, ctrl=ctl, print=T)

ret_sgRSEA <- sgRSEA(dat, r.seed = 123, multiplier = multiplier)$gene.neg %>% as.data.frame

print(ret_sgRSEA %>% head)

ret_sgRSEA %>% rownames_to_column("gene") %>% 
    select(gene = gene, fdr = FDR.neg, stat = NScore) %>%
    write_csv(outfile)


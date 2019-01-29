#!/usr/bin/env Rscript
library(tidyverse)
library(ScreenBEAM)
library(glue)
library(here)

source(here("utils", "parse_args.R"))

tmp_input <- tempfile()
tmp_output <- tempfile()

df_count %>% select(sgRNA = sgRNA, gene=gene, everything()) %>% write_delim(tmp_input, delim="\t")


burnin <- 5000
nitt <- 150000

eval(parse(text=params))


###NGS data
r<-ScreenBEAM(
  input.file=tmp_input,
  control.samples=ctl,
  case.samples=trt,
  control.groupname='ctl',
  case.groupname='trt',
  data.type='NGS',
  do.normalization=TRUE,
  filterLowCount=TRUE,
  filterBy = 'control',
  count.cutoff=4,
  ###Bayesian computing
  nitt=nitt,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=burnin#number of burnin in MCMC sampling, 5000 is default
)

print(r %>% head)

r %>% select(gene, fdr = 7, beta = 4, stat = 6) %>% 
    mutate(fdr = ifelse(beta>0, 1, fdr)) %>%
    mutate(stat = ifelse(beta>0, 1, stat)) %>%
    select(gene, fdr, stat) %>%
    write_csv(outfile)

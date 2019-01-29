#!/usr/bin/env Rscript
library(tidyverse)
library(CB2)
library(PBNPA)
library(here)
source(here("utils", "parse_args.R"))

df_count <- as.data.frame(df_count)

if(length(ctl)==1) {
    ctl <- rep(ctl, length(trt))
}
datlist <- list()
for(i in 1:length(ctl)) {
    datlist[[i]] <- data.frame(sgRNA = df_count$sgRNA,
                               Gene = df_count$gene,
                               initial.count = df_count[,ctl[i]],
                               final.count = df_count[,trt[i]])
}

sim.no <- 10

eval(parse(text=params))
print(sim.no)

result <- PBNPA(datlist, sim.no = sim.no)$final.result
print(result %>% head)

result %>% 
    select(gene = Gene, fdr = neg.fdr, stat = neg.pvalue) %>%
    write_csv(outfile)

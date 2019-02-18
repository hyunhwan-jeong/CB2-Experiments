#!/usr/bin/env Rscript
library(tidyverse)
library(here)
library(CRISPhieRmix)
source(here("utils", "parse_args.R"))
eval(parse(text=params))

library(DESeq2)
df_design <-
    rbind(
        data.frame(group = "ctl", sample_name = ctl, stringsAsFactors = FALSE),
        data.frame(group = "trt", sample_name = trt, stringsAsFactors = FALSE)
    )


df_count <- df_count %>% unite("id", c("gene", "sgRNA")) %>% column_to_rownames("id")

dds <- DESeqDataSetFromMatrix(
    df_count,
    df_design,
    ~group)
dds <- DESeq(dds)
df_ret <- results(dds) %>% as.data.frame %>% rownames_to_column("id")
print(df_ret %>% head)

if(!exists("neg_prefix")) {
    neg_prefix <- "NEGCTRL"
}

df_ret$log2FoldChange[is.na(df_ret$log2FoldChange)] <- 0
neg_ctrls <- df_ret %>% pull(id) %>% startsWith(neg_prefix, .)

neg_log2FC <- df_ret$log2FoldChange[neg_ctrls]
df_ret <- df_ret %>% filter(!neg_ctrls) %>% separate(id, c("gene", "seq"), sep = "_") 
print(df_ret %>% head)

genes <- factor(df_ret$gene, levels = unique(df_ret$gene))
log2FC <- df_ret$log2FoldChange

if(length(neg_log2FC) == 0) {
    neg_log2FC <- NULL
}
ret <- CRISPhieRmix(x=log2FC,  geneIds = genes, negCtrl = neg_log2FC, PLOT=F, VERBOSE = T)

#print(ret$genes)
#print(ret$FDR)

tibble(gene = ret$genes, fdr = ret$FDR, stat = ret$locfdr) %>% write_csv(outfile)

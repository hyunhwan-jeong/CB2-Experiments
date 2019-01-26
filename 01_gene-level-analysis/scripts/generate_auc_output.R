#!/usr/bin/env Rscript
library(argparse)
library(tidyverse)
library(CB2)

parser <- ArgumentParser()
parser$add_argument("infile", nargs=1, help="file has been analyzed.")
parser$add_argument("-e", "--essential", help="File Path of the essential gene list.",  required=T)
parser$add_argument("-n", "--nonessential", help="File Path of the non-essential gene list.",  required=T)
parser$add_argument("-o", "--outfile", help="File Path to the output",  required=T)

args <- parser$parse_args()
infile <- args$infile
outfile <- args$outfile

e <- scan(args$essential, what = "character")
n <- scan(args$nonessential, what = "character")

df_fdr <- read_csv(infile)
gene <- df_fdr$gene
fdr <- df_fdr$fdr

df_auc <- rbind(
  data.frame(gene = gene[gene %in% e], fdr = fdr[gene %in% e], category = 1),
  data.frame(gene = gene[gene %in% n], fdr = fdr[gene %in% n], category = 0))

#if(sum(!(e %in% gene))) {
#  df_auc <- rbind(df_auc, data.frame(gene=e[!(e %in% gene)], fdr=NA, category = 1))
#}

#if(sum(!(n %in% gene))) {
#  df_auc <- rbind(df_auc, data.frame(gene=n[!(n %in% gene)], fdr=NA, category = 0))
#}

write_csv(df_auc, outfile)

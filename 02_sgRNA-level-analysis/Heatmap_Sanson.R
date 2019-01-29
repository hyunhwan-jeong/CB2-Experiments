library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)

CRISPRn <- list()
read_delim("dat/CRISPRn_readcount.txt", delim = "\t") %>%
  rename(sgRNA = sgRNA,
         gene = Gene,
         "plasmid" = pDNA,
         "Rep A" = RepA,
         "Rep B" = RepB,
         "Rep C" = RepC) %>% unite("sgRNA", c("gene", "sgRNA")) -> CRISPRn$readcount

CRISPRi <- list()
read_delim("dat/CRISPRi_readcount.txt", delim = "\t") %>%
  rename(sgRNA = sgRNA,
         gene = Gene,
         "plasmid" = pDNA,
         "Rep A" = RepA,
         "Rep B" = RepB,
         "Rep C" = RepC) %>% unite("sgRNA", c("gene", "sgRNA")) -> CRISPRi$readcount


read_and_merge <- function(obj, screen) {
  sprintf("dat/%s_cb2_sgrna.csv", screen) %>% read_csv() -> obj$cb2_sg
  sprintf("dat/%s_mageck_sgrna.txt", screen) %>% read.table(sep="\t", stringsAsFactors = F, header = T) %>% unite("id", c("Gene", "sgrna")) -> obj$mageck_sg
  obj$cb2_sg %>% left_join(obj$mageck_sg, by = c("sgRNA"="id")) -> obj$merged_sg
  obj
}

CRISPRn %>% read_and_merge("CRISPRn") -> CRISPRn
CRISPRi %>% read_and_merge("CRISPRi") -> CRISPRi

normalize <- function(obj) {
  for(i in 3:ncol(obj$readcount)) {
    s <- sum(obj$readcount[,i])
    obj$readcount[,i] <- obj$readcount[,i] / s * 10^6
  }
  obj
}

CRISPRn %>% normalize() -> CRISPRn
CRISPRi %>% normalize() -> CRISPRi

plot_heatmap <- function(obj, df_sg, main_title) {
  dplyr::select(df_sg, sgRNA) %>% 
    left_join(obj$readcount, by = "sgRNA") %>% 
    as.data.frame %>%
    column_to_rownames("sgRNA") -> df_rc 
  print(df_rc)
  df_rc %>% 
    `+`(1) %>% log2 %>%
    pheatmap(
      #cluster_cols = F,
      show_rownames = F,
      border_color = NA,
      legend = T,
      treeheight_row = 0,
      main = main_title,
      silent = T) 
}

generate_figure <- function(obj, cutoff = 0.01, labels = c("A", "B")) {
  df_sg <- obj$merged_sg  %>% 
    filter(p_ts < cutoff, p.twosided > cutoff) 
  obj$merged_sg  %>% 
    filter(p_ts < cutoff, p.twosided > cutoff) %>% 
    plot_heatmap(obj, ., "CB2") -> hm.cc
  obj$merged_sg %>% 
    filter(p_ts > cutoff, p.twosided < cutoff) %>% 
    plot_heatmap(obj, ., "MAGeCK") -> hm.mg
  
  obj$merged_sg %>% 
    select(CB2 = p_ts, MAGeCK = p.twosided) %>%
    mutate(CB2 = CB2 < cutoff,
           MAGeCK = MAGeCK < cutoff) %>% 
    euler %>% 
    plot(quantities=T,
         fill = c("#f8e4e7", "#fbfbe5"),
         edges = c("#b83f40", "#bfc160")) -> vd
  plot_grid(plot_grid(hm.cc$gtable, hm.mg$gtable, scale = 0.8, nrow = 1),
            vd,   
            nrow=1,
            rel_widths = c(2.5,1),
            labels = labels)
}


plot_grid(
  CRISPRn %>% generate_figure(labels = c("A","B")),
  CRISPRi %>% generate_figure(labels = c("C","D")),
  nrow = 2,
  labels = "AUTO") %>% save_plot("figures/heatmap_Sanson.png", ., base_width = 6, base_height = 8)


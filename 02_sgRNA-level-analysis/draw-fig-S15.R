library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)
library(here)
CRISPR.RT112 <- list()
CRISPR.UMUC3 <- list()
CRISPRi.RT112 <- list()
read_csv("dat/CRISPR.RT112.csv") %>% dplyr::rename(gene=Gene) -> CRISPR.RT112$readcount
read_csv("dat/CRISPR.UMUC3.csv") -> CRISPR.UMUC3$readcount
read_csv("dat/CRISPRi.RT112.csv") -> CRISPRi.RT112$readcount

rename_cols <- function(obj) {
  obj$readcount %>%
    dplyr::rename("T0 1" = A1,
                  "T0 2" = A2,
                  "T0 3" = A3,
                  "T1 1" = B1,
                  "T1 2" = B2,
                  "T1 3" = B3) -> obj$readcount
  obj
}

CRISPR.RT112 %>% rename_cols() -> CRISPR.RT112
CRISPR.UMUC3 %>% rename_cols() -> CRISPR.UMUC3
CRISPRi.RT112 %>% rename_cols() -> CRISPRi.RT112

read_and_merge <- function(obj, screen) {
  sprintf("dat/evers_sgoutput/%s/CB2_sgRNA.csv", screen) %>% read_csv() -> obj$cb2_sg
  sprintf("dat/evers_sgoutput/%s/MAGeCK_sgRNA.csv", screen) %>%
    read_csv() -> obj$mageck_sg
  obj$cb2_sg %>% left_join(obj$mageck_sg,
                           by = c("sgRNA"="sgrna")) -> obj$merged_sg
  obj
}


CRISPR.RT112 %>% read_and_merge("CRISPR.RT112") -> CRISPR.RT112
CRISPR.UMUC3 %>% read_and_merge("CRISPR.UMUC3") -> CRISPR.UMUC3
CRISPRi.RT112 %>% read_and_merge("CRISPRi.RT112") -> CRISPRi.RT112

normalize <- function(obj) {
  for(i in 3:ncol(obj$readcount)) {
    s <- sum(obj$readcount[,i])
    obj$readcount[,i] <- obj$readcount[,i] / s * 10^6
  }
  obj
}

CRISPR.RT112 %>% normalize() -> CRISPR.RT112
CRISPR.UMUC3 %>% normalize() -> CRISPR.UMUC3
CRISPRi.RT112 %>% normalize() -> CRISPR.RT112i

plot_heatmap <- function(obj, df_sg, main_title) {
  if(nrow(df_sg) == 0) {
    return(NULL)
  }
  dplyr::select(df_sg, sgRNA) %>%
    left_join(obj$readcount, by = "sgRNA") %>%
    as.data.frame %>%
    column_to_rownames("sgRNA") %>%
    dplyr::select(-gene) -> df_rc
  df_rc %>%
    `+`(1) %>% log2 %>%
    pheatmap(
      scale = "row",
      #cluster_cols = F,
      #cluster_rows = F,
      show_rownames = F,
      border_color = NA,
      legend = F,
      #annotation_row = df_an,
      clustering_method = "average",
      #clustering_distance_cols = "correlation",
      treeheight_row = 0,
      #treeheight_col = 0,
      main = main_title,
      silent = T) #%>% .$gtable
}

generate_figure <- function(obj, cutoff = 0.01, labels = c("A", "B")) {
  #df_sg <- obj$merged_sg  %>%
  #  filter(p_value_neg < cutoff, p.low > cutoff)
  obj$merged_sg  %>%
    filter(p_value_neg < cutoff, p.low < cutoff) %>%
    plot_heatmap(obj, ., "CC2 & MAGeCK") -> hm.cc
  obj$merged_sg %>%
    filter(p_value_neg > cutoff, p.low < cutoff) %>%
    plot_heatmap(obj, ., "MAGeCK") -> hm.mg
  
  obj$merged_sg %>%
    select(CC2 = p_value_neg, MAGeCK = p.low) %>%
    mutate(CC2 = CC2 < cutoff,
           MAGeCK = MAGeCK < cutoff) %>%
    euler %>%
    plot(quantities=T,
         fill = c("#f8e4e7", "#fbfbe5"),
         edges = c("#b83f40", "#bfc160")) -> vd
  plot_grid(plot_grid(hm.cc$gtable, hm.mg$gtable, scale = 0.8, nrow = 1),
            plot_grid(vd, scale=0.8),
            nrow=1,
            rel_widths = c(2,1),
            labels = labels)
}

CRISPRi.RT112 %>% generate_figure()

plot_grid(
  CRISPR.RT112 %>% generate_figure(labels = c("A", "B")),
  CRISPR.UMUC3 %>% generate_figure(labels = c("C", "D")),
  CRISPRi.RT112 %>% generate_figure(labels = c("E", "F")),
  ncol=1) %>%
  save_plot("figures/fig-S15.pdf",., base_width = 8, base_height = 12)


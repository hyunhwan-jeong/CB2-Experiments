library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)

Gsk3 <- list()
read_delim("dat/Gsk3_readcount.txt", delim = "\t") %>%
  rename(sgRNA = gRNA,
         gene = Gene,
         "Presort 1" = Before_1,
         "Presort 2" = Before_2,
         "Presort 3" = Before_3,
         "Presort 4" = Before_4,
         "GFP+ 1" = Pos_1,
         "GFP+ 2" = Pos_2,
         "GFP+ 3" = Pos_3,
         "GFP+ 4" = Pos_4) -> Gsk3$readcount

Tnf <- list()
read_delim("dat/Tnf_readcount.tsv", delim = "\t")  %>% 
  select(-starts_with("Pre_")) %>%
  rename("TNF-Low 1" = Low_TNF_EX1,
         "TNF-Low 2" = Low_TNF_EX2,
         "TNF-Low 3" = Low_TNF_EX3,
         "TNF-High 1" = High_TNF_EX1,
         "TNF-High 2" = High_TNF_EX2,
         "TNF-High 3" = High_TNF_EX3) -> Tnf$readcount

summary(Gsk3$readcount)
summary(Tnf$readcount)

read_and_merge <- function(obj, screen) {
  sprintf("dat/%s_cb2_sgrna.csv", screen) %>% read_csv() -> obj$cb2_sg
  sprintf("dat/%s_mageck_sgrna.txt", screen) %>% read.table(sep="\t", stringsAsFactors = F, header = T) -> obj$mageck_sg
  obj$cb2_sg %>% left_join(obj$mageck_sg, by = c("sgRNA"="sgrna")) -> obj$merged_sg
  obj
}

Gsk3 %>% read_and_merge("Gsk3") -> Gsk3
  Tnf %>% read_and_merge("Tnf") -> Tnf

normalize <- function(obj) {
  for(i in 3:ncol(obj$readcount)) {
    s <- sum(obj$readcount[,i])
    obj$readcount[,i] <- obj$readcount[,i] / s * 10^6
  }
  obj
}

Gsk3 %>% normalize() -> Gsk3
Tnf %>% normalize() -> Tnf

plot_heatmap <- function(obj, df_sg, main_title) {
  dplyr::select(df_sg, sgRNA) %>% 
    left_join(obj$readcount, by = "sgRNA") %>% 
    as.data.frame %>%
    column_to_rownames("sgRNA") %>%
    dplyr::select(-gene) -> df_rc 
  df_rc %>% 
    `+`(1) %>% log2 %>%
    pheatmap(
      #scale = "row",
      cluster_cols = F,
      #cluster_rows = F,
      show_rownames = F,
      border_color = NA,
      legend = F,
      #annotation_row = df_an,
      #clustering_method = "average",
      #clustering_distance_cols = "correlation",
      treeheight_row = 0,
      #treeheight_col = 0,
      main = main_title,
      silent = T) #%>% .$gtable
}

generate_figure <- function(obj, cutoff = 0.01, labels = c("A", "B")) {
  df_sg <- obj$merged_sg  %>% 
    filter(p_value_twosided < cutoff, p.twosided > cutoff) 
  obj$merged_sg  %>% 
    filter(p_value_twosided < cutoff, p.twosided > cutoff) %>% 
    plot_heatmap(obj, ., "cb2") -> hm.cc
  obj$merged_sg %>% 
    filter(p_value_twosided > cutoff, p.twosided < cutoff) %>% 
    plot_heatmap(obj, ., "MAGeCK") -> hm.mg
  
  obj$merged_sg %>% 
    select(cb2 = p_value_twosided, MAGeCK = p.twosided) %>%
    mutate(cb2 = cb2 < cutoff,
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


#Gsk3 %>% generate_figure() %>% save_plot("Gsk3.png", ., base_width = 6, base_height = 4)
#Tnf %>% generate_figure() %>% save_plot("Tnf.png", ., base_width = 6, base_height = 4)

plot_grid(
  Gsk3 %>% generate_figure(labels = c("A","B")), 
  Tnf %>% generate_figure(labels = c("C", "D")),
  nrow = 2,
  labels = "AUTO") %>% save_plot("figures/heatmap.png", ., base_width = 6, base_height = 8)

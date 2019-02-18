setwd("01_gene-level-analysis/")
lineplot_perf_cb2 <- function(all_df, col_ord) {
  require(tidyverse)
  require(cowplot)
  ct <- c(0.1, 0.05, 0.01, 0.005, 0.001)
  df.prof <- tibble()
  dset <- "CRISPRn-RT112"
  #mat <- "CB2"
  for(mat in unique(all_df$method)) {
    tmp <- all_df %>% filter(method == mat)
    for(fdr in ct) {
      #fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$essential == 1))
      FP <- sum((tmp$fdr <= fdr) & (tmp$essential == 0))
      FN <- sum((tmp$fdr > fdr) & (tmp$essential == 1))
      
      precision <- TP / max(1,(TP+FP))
      recall <- TP / max(1,(TP+FN))
      fmeasure <- 2*(precision*recall)/(precision+recall)
      if(is.na(fmeasure)) fmeasure <- 0
      df.prof <- bind_rows(df.prof, tibble(dataset=dset, method=mat, FDR=fdr, value=precision, measure="precision", TP=TP, FP=FP, FN=FN))
      df.prof <- bind_rows(df.prof, tibble(dataset=dset, method=mat, FDR=fdr, value=recall, measure="recall", TP=TP, FP=FP, FN=FN))
      df.prof <- bind_rows(df.prof, tibble(dataset=dset,method=mat, FDR=fdr, value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
  print(df.prof)
  df.prof <- df.prof %>% as.data.frame
  df.prof$method <- df.prof$method %>% factor
  
  tmp <- df.prof %>% 
    mutate(FDR = factor(FDR, levels=ct)) %>%
    mutate(dataset = factor(dataset, levels = col_ord))
  
  (pt <- tmp %>%
      ggplot(aes(x=FDR, y=value)) +
      geom_point(aes(colour=method), size=2) +
      geom_line(aes(colour=method, group=method), size=0.8) +
      geom_point(data = tmp %>% filter(method == "CB2"), aes(colour=method), size=2) +
      geom_line(data = tmp %>% filter(method == "CB2"), aes(colour=method, group=method), size=0.8) +
      facet_grid(dataset~measure) + ylim(0,1) +
      xlab("FDR") + ylab("measure") +
      #scale_color_brewer(palette = "RdYlBu") +
      theme(axis.text.x = element_text(angle=90))) +
    theme(strip.text.x = element_text(face="bold"))
}


library(tidyverse)
library(RobustRankAggreg)
e <- scan("data/Evers/essential-genes.txt", what = "character")
cb2 <- read_tsv("RRA_tmp/cb_rra_out.txt") %>% select(gene = group_id, fdr=FDR) %>% 
  mutate(essential = gene %in% e) %>%
  mutate(method = "CB2-RRA") 
mageck <- read_tsv("RRA_tmp/sample1.gene.low.txt") %>% select(gene = group_id, fdr=FDR) %>% 
  mutate(essential = gene %in% e) %>% 
  mutate(method = "MAGeCK")

read_tsv("RRA_tmp/sample1.plow.txt") %>% select(sgrna, symbol, plow = p.low) %>% 
  mutate(rank = rank(plow) / n()) %>% 
  select(-plow) %>% 
  spread(sgrna, rank) %>% column_to_rownames("symbol") %>% as.matrix -> r

# mageck_rra <- aggregateRanks(rmat = r, method = "RRA") %>%
#   remove_rownames() %>% 
#   select(gene = Name, fdr = Score) %>% 
#   mutate(fdr = p.adjust(fdr, method = "fdr")) %>% 
#   mutate(method = "MAGeCK_RRA") %>% 
#   mutate(essential = gene %in% e)


read_tsv("RRA_tmp/cb_rra.txt") %>% select(sgrna, symbol, plow) %>% 
  mutate(rank = rank(plow) / n()) %>% 
  select(-plow) %>% 
  spread(sgrna, rank) %>% column_to_rownames("symbol") %>% as.matrix -> r

cb2_rra <- aggregateRanks(rmat = r, method = "RRA") %>%
  remove_rownames() %>% 
  select(gene = Name, fdr = Score) %>% 
  mutate(fdr = p.adjust(fdr, method = "fdr")) %>% 
  mutate(method = "CB2_RRA") %>% 
  mutate(essential = gene %in% e)

cb2_org <- read_csv("results/Evers/CRISPRn-RT112/AUC/CB2.csv") %>% 
  select(gene, fdr) %>% 
  mutate(method = "CB2") %>% 
  mutate(essential = gene %in% e)

bind_rows(cb2_org, cb2, mageck) %>% 
  as.data.frame() %>% 
  lineplot_perf_cb2(c("CRISPRn-RT112")) +
  theme(legend.position = "bottom")



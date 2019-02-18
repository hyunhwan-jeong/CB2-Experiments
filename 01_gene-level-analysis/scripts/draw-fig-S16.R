library(tidyverse)
library(glue)
library(pheatmap)
library(GGally)
library(here)
library(cowplot)
source(here("utils/draw_auc_curve.R"))


draw_f1 <- function(all.df, col_ord = NULL) {
  require(tidyverse)
  require(cowplot)
  ct <- c(0.1, 0.05, 0.01, 0.005, 0.001)
  df.prof <- tibble()
  for(pset in unique(all.df$param)) {
    tmp <- all.df %>% filter(param==pset)
    for(fdr in ct) {
      #fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$category == 1))
      FP <- sum((tmp$fdr <= fdr) & (tmp$category == 0))
      FN <- sum((tmp$fdr > fdr) & (tmp$category == 1))
      
      precision <- TP / max(1,(TP+FP))
      recall <- TP / max(1,(TP+FN))
      fmeasure <- 2*(precision*recall)/(precision+recall)
      if(is.na(fmeasure)) fmeasure <- 0
      df.prof <- bind_rows(df.prof, tibble(param=pset, FDR=fdr, value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
  
  df.prof <- df.prof %>% as.data.frame

  tmp <- df.prof %>% 
    mutate(FDR = factor(FDR, levels=ct)) %>%
    mutate(param = factor(param, levels = col_ord))
  
  (pt <- tmp %>%
      ggplot(aes(x=FDR, y=value)) +
      geom_point(aes(colour=param), size=2) +
      geom_line(aes(colour=param, group=param), size=0.8) +
      ylim(0,1) +
      xlab("FDR") + ylab("F1-Score") +
      scale_color_brewer(palette = "RdYlBu") +
      theme(axis.text.x = element_text(angle=90))) +
    theme(strip.text.x = element_text(face="bold"))
}

show_tune_results <- function(method, fmt, params) {
  egene <- scan(here("data/Evers/essential-genes.txt"), what="character")
  
  df_dat <- lapply(params, 
                   function(p) read_csv(here(fmt %>% glue)) %>% 
                     mutate(param = p)) %>% bind_rows()
  
  df_dat %>% select(-stat) %>% 
    mutate(category = gene %in% egene) -> df_all
  
  title <- ggdraw() + draw_label(method, fontface='bold')
  
  p1 <- draw_f1(df_all, params)
  
  (p2 <- df_dat %>% select(-stat) %>%
      mutate(fdr = ifelse(fdr <= .Machine$double.eps, 1e-10, fdr)) %>% 
      mutate(fdr = -log10(fdr)) %>%
      spread(param, fdr) %>%
      column_to_rownames("gene") %>% ggpairs() + xlab("-log10(FDR)") + ylab("-log10(FDR)") + theme(axis.text.x = element_text(angle = 90)))
  
  plot_grid(title, p1, p2 %>% ggmatrix_gtable, ncol = 1, rel_heights = c(0.5,3,5))
}

plot_grid(
  plotlist = list(
    MAGeCK = show_tune_results("MAGeCK", "results_params_tune/mageck-{p}.txt", c(100, 1000, 10000, 100000) %>% as.integer),
    ScreenBEAM = show_tune_results("ScreenBEAM", "results_params_tune/ScreenBEAM-{p}-15000.txt", c(50, 500, 5000, 10000) %>% as.integer),
    PBNPA = show_tune_results("PBNPA", "results_params_tune/PBNPA-{p}.txt", c(10, 50, 100, 500) %>% as.integer),
    sgRSEA = show_tune_results("sgRSEA", "results_params_tune/sgRSEA-{p}.txt", c(10, 50, 100, 500) %>% as.integer),
    RIGER = show_tune_results("RIGER", "results_params_tune/RIGER-{p}.txt", c("0.1", "0.5", "1.0", "1.5", "2.0"))
  ),
  nrow = 1
)


save_plot(here("figures/fig-16.pdf"), last_plot(), base_width = 22, base_height = 10)

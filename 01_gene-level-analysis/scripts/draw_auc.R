library(tidyverse)
library(precrec)
library(ggsci)
library(gghighlight)
library(cowplot)
library(glue)
library(here)
source(here("utils/draw_auc_curve.R"))

draw_auc_figures <- function(study, dataset, methods, title, filename) {
  df_dataset <- lapply(methods, function(m) 
    read_csv("results/{study}/{dataset}/AUC/{m}.csv" %>% glue) %>%
      mutate(method = m) %>% 
      mutate(score = 1-fdr)) %>% bind_rows()  
  
  df_dataset %>% select(-fdr) %>% spread(method, score) -> df_dataset
  
  df_dataset <- df_dataset[df_dataset %>% is.na() %>% rowSums() != length(methods),]
  df_dataset[is.na(df_dataset)] <- 0
  df_dataset <- df_dataset %>% gather(method, score, -gene, -category)  
  scores <- join_scores (
    lapply(methods, function(m) 
      df_dataset %>% filter(method == m) %>% pull(score)),
    chklen = F
  )
  
  labels <- join_scores (
    lapply(methods, function(m) 
      df_dataset %>% filter(method == m) %>% pull(category)),
    chklen = F
  )
  
  mmdat <- mmdata(scores, 
                  labels, 
                  modnames = methods)
  
  (eval_ret <- evalmod(mmdat))
  
  
  p1 <- draw_auc_curve(eval_ret, "PRC") + xlab("Recall") + ylab("Precision") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("PR Curves") + theme(axis.text.x = element_text(angle = 90))
  p2 <- draw_auc_curve(eval_ret, "ROC") + xlab("1- Specificity") + ylab("Sensitivity") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("ROC Curves") + theme(axis.text.x = element_text(angle = 90))
  p_title <- ggdraw() + draw_label(title, fontface='bold')
  p3 <- plot_grid(p_title, plot_grid(p1, p2, ncol=1, labels = "auto"), rel_heights = c(0.1,1), ncol=1)
  save_plot(filename = filename, p3, base_width = 14, base_height = 8)  
} 

library(tidyverse)
library(precrec)
library(ggsci)
library(gghighlight)
library(cowplot)
library(glue)
library(here)
source(here("utils/draw_auc_curve.R"))

dataset <- "CRISPRn-A375"
methods <- c("CB2", "ScreenBEAM", "MAGeCK", "PBNPA", "RSA", "sgRSEA", "HitSelect")


scores <- join_scores (
  lapply(methods, function(m) 
    read_csv("results/Sanson/{dataset}/AUC/{m}.csv" %>% glue) %>%
      mutate(fdr = 1-fdr) %>%
      pull(fdr)),
  chklen = F
)

labels <- join_scores (
  lapply(methods, function(m) 
    read_csv("results/Sanson/{dataset}/AUC/{m}.csv" %>% glue) %>%
      pull(category)),
  chklen = F
)

mmdat <- mmdata(scores, 
                labels, 
                modnames = methods)

(eval_ret <- evalmod(mmdat))


p1 <- draw_auc_curve(eval_ret, "PRC") + xlab("Recall") + ylab("Precision") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("PR Curves")
p2 <- draw_auc_curve(eval_ret, "ROC") + xlab("1- Specificity") + ylab("Sensitivity") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("ROC Curves")
title <- ggdraw() + draw_label("Genome-wide CRISPRn library (Brunello)", fontface='bold')
p3 <- plot_grid(title, plot_grid(p1, p2, ncol=1, labels = "auto"), rel_heights = c(0.1,1), ncol=1)

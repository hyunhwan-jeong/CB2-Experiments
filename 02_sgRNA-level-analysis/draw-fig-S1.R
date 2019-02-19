library(tidyverse)
library(cowplot)
library(here)
draw_mean_variance_plot <- function(df_count, title) {
  require(tidyverse)
  require(matrixStats)
  df_norm <- df_count
  
  for(i in 1:ncol(df_norm)) {
    df_norm[,i] <- df_count[,i] / sum(df_count[,i]) * 10^6
  }
  
  means <- rowMeans(df_norm)
  vars <- rowVars(df_norm %>% as.matrix)
  
  tibble(m = log10(means), v = log10(vars)) %>% 
    ggplot(aes(x=m,y=v)) +
    geom_point(alpha=0.5, size = 0.5, color = "grey") + 
    geom_abline(slope = 1, intercept = 0, color="red") + 
    geom_smooth() +
    ggtitle(title) + xlab("log10(meanCPM)") + ylab("log01(varianceCPM)")
}

p1 <- read_csv("../01_gene-level-analysis/data/Evers/CRISPRn-RT112.csv")  %>% 
  select(sgRNA, B1,B2, B3) %>% column_to_rownames("sgRNA") %>% draw_mean_variance_plot(title = "CRISPRn-RT112 (T0)")

p2 <- read_tsv("../01_gene-level-analysis/data/Sanson/CRISPRn-A375.tsv")  %>% 
  select(sgRNA, RepA, RepB, RepC) %>% column_to_rownames("sgRNA") %>% draw_mean_variance_plot(title = "CRISPRn-A375")

p3 <- read_tsv("../02_sgRNA-level-analysis/dat/Gsk3_readcount.txt")  %>% 
  select(gRNA, Before_1, Before_2, Before_3, Before_4) %>% column_to_rownames("gRNA") %>% draw_mean_variance_plot(title = "Gsk3 (Before)")

plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO")
save_plot("figures/fig-S1.png", last_plot(), base_height = 3, base_width = 8)

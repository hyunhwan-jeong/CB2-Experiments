library(tidyverse)
library(ggsci)
library(pheatmap)
library(cowplot)
library(glue)
library(here)

draw <- function(screen) {
  methods <- c("CB2", "MAGeCK", "PBNPA", "sgRSEA", "RSA", "RIGER", "ScreenBEAM")
  methods <- c(methods, "PinAPL-Py", "HitSelect", "CRISPhieRmix")
  
  df_rank <- lapply(methods, 
                    function(m) read_csv("results/Evers/{screen}/FDR/{m}.csv" %>% glue %>% here) %>% 
                      mutate(method = m, rank = rank(stat)) %>%
                      select(method, gene, rank)) %>% bind_rows()

  df_rank %>% spread(method, rank) -> df_matrix
  df_matrix %>% select(-gene) %>% cor(use="complete.obs") %>% pheatmap(display_numbers = T, silent=F, fontsize_number = 12) %>% .$gtable -> hm
  df_matrix %>% select(-gene) %>% cor(use="complete.obs") -> tmp
  tmp[lower.tri(tmp)] %>% summary %>% print
  e <- scan("data/Evers/essential-genes.txt", what="character")
  
  df_rank %>% 
    mutate(essential = ifelse(gene %in% e, "essential", "non-essential")) %>% 
    ggplot(aes(x=method, y=rank)) +
    geom_boxplot(aes(fill=essential), width = 0.25, alpha=0.7) +
    geom_jitter(aes(fill=essential, color=essential), width = 0.25, alpha=0.5) + scale_y_reverse() +
    scale_fill_npg() + scale_color_npg() + ggtitle(screen) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) -> gg_bo
  
  if(screen == "CRISPRn-RT112")
    plot_grid(gg_bo, hm, ncol=1, labels ="AUTO")
  else
    plot_grid(gg_bo, hm, ncol=1)
}


plot_grid(
  "CRISPRn-RT112" %>% draw,
  "CRISPRn-UMUC3" %>% draw,
  "CRISPRi-RT112" %>% draw,
  nrow = 1
)

#save_plot(plot = last_plot(), filename = "figures/fig_S1.pdf", base_width = 16, base_height = 10)

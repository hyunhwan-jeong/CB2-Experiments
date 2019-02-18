library(cowplot)
library(ggsci)
library(glue)
library(here)
library(tidyverse)
plot.dot.bar <- function(gene.name,
                         dataset.name,
                         ncol) {
  methods <- c("CB2", "MAGeCK", "PBNPA", "sgRSEA", "RSA", "RIGER", "ScreenBEAM")
  methods <- c(methods, "PinAPL-Py", "HitSelect")
  
  df_fdr <- lapply(methods, 
                    function(m) 
                      read_csv("results/Evers/{dataset.name}/FDR/{m}.csv" %>% glue %>% here) %>% 
                     mutate(method = m)) %>% 
    bind_rows() %>%
    as.data.frame
  df_count <- read_csv(here("data/Evers/read_count.csv"))
  
  df_fdr$method <- df_fdr$method %>% factor
  lev <- df_fdr$method %>% levels
  lev[lev=="CB2"] <- expression("CB"^2)
  plot.bar <-
    df_fdr %>%
    filter(gene==gene.name) %>%
    mutate(FDR=-log10(fdr)) %>%
    ggplot(aes(x=method, y=FDR)) +
    geom_bar(stat = "identity", aes(fill=method)) +
    ylab("-log10(FDR)") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = -log10(0.05), color = "black", alpha=0.5) +
    theme(axis.title.x = element_blank()) +
    scale_x_discrete(labels=lev) +
    ylim(0,20) +
    scale_fill_brewer(palette = "RdYlBu") 
  
  
  plot.dots <-
    df_count %>%
    filter(gene==gene.name,
           dataset==dataset.name,
           count_type=="countPerMil") %>%
    mutate(CPM=read_count) %>%
    ggplot(aes(x=group, y=CPM)) +
    geom_dotplot(aes(fill=group, color=group),
                 binaxis='y',
                 stackdir='center',
                 stackratio=1.5,
                 dotsize=1.2) +
    facet_wrap(~sgRNA, scales = "free_y", ncol=ncol) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_color_npg() +
    scale_fill_npg() +
    ggtitle(str_c( gene.name, " (", dataset.name %>% str_replace("\\.", "-"), ")")) +
    theme(axis.title.x = element_blank())
  
  plot_grid(plot.dots, plot.bar, ncol = 1)
}

plot_dot_and_bar <- function(gene.name, file.name) {
  plots <- list()
  for(dataset in c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112")) {
    plots[[dataset]] <- plot.dot.bar(gene.name, dataset, 3)
  }
  
  plot_merged <- plot_grid(plotlist = plots, labels = "AUTO", nrow=1)
  save_plot(file.name, plot_merged, base_width = 14, base_height = 12)
}

plot_dot_and_bar("RPL5", file.name = "figures/fig2.pdf") # Fig. 2
plot_dot_and_bar("COPS8", file.name = "figures/fig-S10.pdf") # Supplemental Fig. 10
plot_dot_and_bar("RPL27", file.name = "figures/fig-S11.pdf") # Supplemental Fig. 11

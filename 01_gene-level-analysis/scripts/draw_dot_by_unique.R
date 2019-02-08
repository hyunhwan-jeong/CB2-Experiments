library(here)
library(tidyverse)
source(here("utils/load_data.R"))
source(here("utils/upset.R"))

df_evers <- load_data("Evers") 

mat_evers <- df_evers %>% select(-stat) %>%  spread(method, fdr) %>% 
  filter(essential, CB2 < 0.01, rowSums(.[,5:12]<0.01)==0) %>%  arrange(dataset)

plot.dot <- function(gene.name,
                         dataset.name,
                         ncol) {
  require(cowplot)
  require(ggsci)
  require(glue)
  require(here)
  require(tidyverse)
  
  df_count <- read_csv(here("data/Evers/read_count.csv"))
  

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
  
  plot.dots  
}

plot_dot <- function(gene.names, dataset) {
  plots <- list()
  for(g in gene.names) {
    plots[[g]] <- plot.dot(g, dataset, 3)
  }
  
  plot_merged <- plot_grid(plotlist = plots)
  #save_plot("figures/fig_dotandbar-{gene.name}.pdf" %>% glue, plot_merged, base_width = 14, base_height = 12)
  plot_merged
}

plot_dot(mat_evers %>% filter(dataset == "CRISPRn-RT112") %>% pull(gene), "CRISPRn-RT112")
plot_dot(mat_evers %>% filter(dataset == "CRISPRn-UMUC3") %>% pull(gene), "CRISPRn-UMUC3")
plot_dot(mat_evers %>% filter(dataset == "CRISPRi-RT112") %>% pull(gene), "CRISPRi-RT112")

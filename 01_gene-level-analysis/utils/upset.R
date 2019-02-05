draw_upset <- function(all_df, dataset) {
  require(here)
  require(UpSetR)
  require(glue)
  require(grid)
  for(d in dataset) {
    pdf(here("figures/UpSet-{d}.pdf" %>% glue), onefile=FALSE)
    fdat <- upset_plots[[d]] <- all_df %>% 
      filter(essential, dataset == d) %>% 
      mutate(fdr = ifelse(fdr <0.1,1,0)) %>% 
      select(method, gene, fdr) %>% spread("method", "fdr") %>% 
      column_to_rownames("gene")
    
    upset(fdat)
    grid.text(d,x = 0.65, y=0.95, gp=gpar(fontsize=20))
    dev.off() 
  }
}

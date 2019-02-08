draw_upset <- function(all_df, dataset) {
  require(here)
  require(UpSetR)
  require(glue)
  require(grid)
  vp <- list()
  for(d in dataset) {
    #pdf(here("figures/UpSet-{d}.pdf" %>% glue), onefile=FALSE, useDingbats=FALSE)
    fdat <- all_df %>% 
      filter(essential, dataset == d) %>% 
      mutate(fdr = ifelse(fdr <0.01,1,0)) %>% 
      select(method, gene, fdr) %>% spread("method", "fdr") %>% 
      column_to_rownames("gene")
    
    upset(fdat, text.scale = 2)
    #grid.text(d,x = 0.65, y=0.95, gp=gpar(fontsize=20))
    grid.edit('arrange',name='arrange2')
    vp[[d]] <- plot_grid(ggdraw() + draw_label(d, fontface = "bold"),  grid.grab(), ncol=1, rel_heights=c(0.1, 1))
    #dev.off() 
  }
  plot_grid(plotlist = vp, labels = "AUTO")
}

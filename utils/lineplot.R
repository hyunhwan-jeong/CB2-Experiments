lineplot_perf <- function(all.df, col_ord = NULL) {
  require(tidyverse)
  require(cowplot)
  ct <- c(0.1, 0.05, 0.01, 0.005, 0.001)
  df.prof <- tibble()
  for(dset in unique(all.df$dataset)) {
    for(mat in unique(all.df$method)) {
      tmp <- all.df %>% filter(dataset==dset, method==mat)
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
  }
  
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
      facet_grid(measure~dataset) + ylim(0,1) +
      xlab("FDR") + ylab("measure") +
      scale_color_brewer(palette = "RdYlBu") +
      theme(axis.text.x = element_text(angle=90))) +
    theme(strip.text.x = element_text(face="bold"))
}

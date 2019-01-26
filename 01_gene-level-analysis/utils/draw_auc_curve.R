draw_auc_curve <- function(eval_ret, ctype = "PRC") {
  df_auprc <- auc(eval_ret) %>% filter(curvetypes == ctype) %>% 
    dplyr::rename(modname = modnames) %>% 
    mutate(method = sprintf("%s (%.3f)", modname, aucs))
  
  df_pr <- eval_ret %>% as.data.frame %>% filter(type==ctype) %>% left_join(df_auprc, by="modname")
  df_pr %>% 
    ggplot(aes(x=x,y=y)) +
    geom_line(aes(color=method), size=1) +
    #geom_line(data = df_pr %>% filter(type==ctype, startsWith(method, "CB2")) , 
    #            aes(x=x,y=y,color=method), size=1) +
    theme_minimal() + scale_color_brewer(palette = "RdYlBu") 
}
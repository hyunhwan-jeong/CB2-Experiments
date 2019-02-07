genes <- df_merge %>% filter(fdr < 0.1,fdr_ts > 0.1) %>% pull(gene)
print(genes %>% head)
print(length(genes))

df_norm <- df_count 
for(j in 1:ncol(df_norm)) {
  df_norm[,j] <- log2(df_norm[,j] / sum(df_norm[,j]) * 10^6+1)
}

df_norm %>% rownames_to_column("id") %>% separate(id, c("gene", "gRNA"), extra = "merge") %>% 
  filter(gene %in% genes) %>% 
  select(-gene) %>% 
  column_to_rownames("gRNA") -> df_heatmap

p1 <- df_heatmap[rowSums(df_heatmap) > 0,] %>% 
  pheatmap(show_rownames = F, scale = "row", cell_width = 10, silent = T)

genes <- df_merge %>% filter(fdr < 0.1,fdr_ts < 0.1) %>% pull(gene)
print(genes %>% head)
print(length(genes))

df_norm <- df_count 
for(j in 1:ncol(df_norm)) {
  df_norm[,j] <- log2(df_norm[,j] / sum(df_norm[,j]) * 10^6+1)
}

df_norm %>% rownames_to_column("id") %>% separate(id, c("gene", "gRNA"), extra = "merge") %>% 
  filter(gene %in% genes) %>% 
  select(-gene) %>% 
  column_to_rownames("gRNA") -> df_heatmap

p2 <- df_heatmap[rowSums(df_heatmap) > 0,] %>% 
  pheatmap(show_rownames = F, scale = "row", cellwidth = 10, silent = T)

library(cowplot)

plot_grid(p2$gtable, p1$gtable, labels = "AUTO", rel_widths = c(1.1,1))
save_plot(filename = "figures/heatmap_Gsk3_with_CRmix.pdf", last_plot())

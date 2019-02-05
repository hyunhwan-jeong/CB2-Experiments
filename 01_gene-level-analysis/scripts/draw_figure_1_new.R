library(tidyverse)
library(here)
library(cowplot)
library(pheatmap)

generate_heatmap_legend_no_essential_2 <- function() {
  df <- data.frame(
    x = 1:5,
    y = 1:1,
    FDR = c("8", "6", "4", "2", "0"),
    #FDR = sprintf("\u2265%d", seq(8,0)),
    w = 1:5,
    stringsAsFactors = F
  )
  (col.pal <- RColorBrewer::brewer.pal(9, "OrRd"))
  #col.pal[1] <- "#FFFFFF"
  
  test <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = FDR), colour = "grey50") +
    scale_fill_manual(
      name = "-log10(FDR)",
      #values = rev(col.pal)[c(1,3,5,7,9)],
      values = rev(col.pal)[c(1,3,5,7,9)],
      limits = df$FDR
    ) + geom_text(aes(label=FDR)) +
    theme(legend.position = "none") + ylab("-log10(FDR)") +
    theme(axis.title.y = element_text(angle=0)) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) + xlab("")
  test
}

files <- Sys.glob("results/Evers/*/FDR/*")

e <- scan("data/Evers/essential-genes.txt", what ="character")
n <- scan("data/Evers/non-essential-genes.txt", what ="character")
mk_df <- function(f) {
  read_csv(f) %>% 
    mutate(
      method = basename(f) %>% str_replace(".csv", ""),
      dataset = f %>% strsplit("/")) %>% 
    mutate(dataset=dataset[[1]][3]) %>%
    mutate(essential = gene %in% e)
}

all_df <- lapply(files, mk_df) %>% bind_rows %>% filter(gene %in% c(e, n))

datasets <- c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112")

pt.f1 <- list()
ct <- c(0.1, 0.05, 0.01, 0.005, 0.001)
for(dset in datasets) {
  df.prof <- tibble()
  for(mat in unique(all_df$method)) {
    tmp <- all_df %>% filter(dataset==dset, method==mat)
    for(fdr in ct) {
      #fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$essential == T))
      FP <- sum((tmp$fdr <= fdr) & (tmp$essential == F))
      FN <- sum((tmp$fdr > fdr) & (tmp$essential == T))
      
      precision <- TP / max(1,(TP+FP))
      recall <- TP / max(1,(TP+FN))
      fmeasure <- 2*(precision*recall)/(precision+recall)
      if(is.na(fmeasure)) fmeasure <- 0
      df.prof <- bind_rows(df.prof, tibble(dataset=dset,method=mat, FDR=fdr, value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
  
  pt.f1[[dset]] <- df.prof %>% mutate(FDR = factor(FDR, levels=ct)) %>%
    ggplot( aes(x=FDR, y=value)) + geom_point(aes(colour=method), size=2) +
    geom_line(aes(colour=method, group=method), size=0.8)+
    xlab("FDR") + ylab("F1-score") + ylim(0,1) +
    #scale_color_npg() +
    scale_color_brewer(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle=90)) +
    theme(strip.text.x = element_text(face="bold")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

legend <- get_legend(pt.f1[[1]] + theme(legend.position = "bottom", legend.justification="center"))
plot_grid(plot_grid(plotlist = pt.f1, nrow=1), legend, ncol=1, rel_heights = c(4,1))


generate_heatmap_legend_no_essential <- function() {
  df <- data.frame(
    x = 1:5,
    y = 1:1,
    FDR = c("8", "6", "4", "2", "0"),
    #FDR = sprintf("\u2265%d", seq(8,0)),
    w = 1:5,
    stringsAsFactors = F
  )
  (col.pal <- RColorBrewer::brewer.pal(9, "Reds"))
  col.pal[1] <- "#FFFFFF"
  
  test <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = FDR), colour = "grey50") +
    scale_fill_manual(
      name = "-log10(FDR)",
      #values = rev(col.pal)[c(1,3,5,7,9)],
      values = rev(col.pal)[c(1,3,5,7,9)],
      limits = df$FDR
    ) + geom_text(aes(label=FDR)) +
    theme(legend.position = "none") + ylab("-log10(FDR)") +
    theme(axis.title.y = element_text(angle=0)) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) + xlab("")
  test
}

prof.level <- "gene"
order.methods <- c("CB2", "PBNPA", "sgRSEA", "HitSelect", "MAGeCK", "RIGER", "RSA", "ScreenBEAM", "PinAPL-Py", "CRISPhieRmix")
all_df$fdr[is.na(all_df$fdr)] <- 1
all_df <- all_df %>% select(-stat)
heatmap <- list()
(col.pal <- RColorBrewer::brewer.pal(9, "OrRd"))
#col.pal[1] <- "#FFFFFF"
for (dset in unique(all_df$dataset)) {
  print(dset)
  ess <- all_df %>% filter(dataset == dset) %>% select(gene, essential) %>% distinct()
  x <-
    all_df %>% filter(dataset == dset) %>% select(-dataset, -essential) %>% spread(method, fdr) %>%
    remove_rownames()
  
  x[x < 1e-8] <- 1e-8
  x[,-1] <- floor(-log10(x[,-1]))
  x <- x[order(-rowSums(x[,-1])),]
  x <- x %>% left_join(ess, by = prof.level)
  x <- x[order(-x$essential),]
  x <- x %>% remove_rownames()
  x <- x[, c(1, 3, 2, 4:ncol(x))]
  x$essential[x$essential] <- 8
  tmp <-
    x %>% mutate(Essentiality = ifelse(essential, "Essential", "Non-essential")) %>%
    select(prof.level, Essentiality) %>%
    column_to_rownames(prof.level)
  
  if(is.null(order.methods)) {
    order.methods <- 1:(ncol(x)-1)
  }
  
  hm <-
    pheatmap(
      column_to_rownames(x, prof.level) %>%
        #filter(essential > 0) %>%
        select(order.methods) %>% t,
      scale = "none",
      cluster_cols = F,
      cluster_rows = F,
      # main = dset,
      color = col.pal,
      legend = F,
      show_rownames = T,
      show_colnames = F,
      annotation_legend = F,
      annotation_col = tmp,
      #border_color = "black",
      annotation_colors = list(
        "Essentiality" = c("Essential" = "#000000", "Non-essential" = "#bdbdbd")
      ),
      cellheight = 12,
    )
  heatmap[[dset]] <- hm$gtable
}
legend_heatmap <- generate_heatmap_legend_no_essential_2()
pt.merged <- list()
for(d in datasets) {
  x <- plot_grid(ggdraw() + draw_label(d, fontface="bold", angle=90), heatmap[[d]], rel_widths = c(1,9)) + theme(plot.margin = margin(l=20, b=20))
  if(d == "CRISPRn-A375") {
    pt.merged[[d]] <- plot_grid(x, pt.f1[[d]], nrow=1, rel_widths =  c(7,3),
                                labels = "AUTO",
                                label_size = 18
    )
    
  } else {
    pt.merged[[d]] <- plot_grid(x, pt.f1[[d]], nrow=1, rel_widths =  c(7,3))
  }
}

legend_f1 <- plot_grid(
  get_legend(pt.f1[[1]] +
               theme(legend.position = "bottom", legend.title = element_blank()) +
               guides(color=guide_legend(ncol=3))
  )
)

top <- plot_grid(plotlist = pt.merged, nrow=3)
bottom <- plot_grid(NULL, plot_grid(legend_heatmap), NULL, legend_f1, nrow=1, rel_widths = c(1,7,2,3))

fig2 <- plot_grid(top,
                  bottom,
                  nrow = 2,
                  rel_heights = c(9,1))
#save_plot("figures/fig2.pdf", fig2, base_width = 14, base_height = 8)

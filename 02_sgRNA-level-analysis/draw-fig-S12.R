library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)
library(ggsci)

CRISPRn <- list()
read_delim("dat/CRISPRn_readcount.txt", delim = "\t") %>%
  rename(sgRNA = sgRNA,
         gene = Gene,
         "plasmid" = pDNA,
         "Rep A" = RepA,
         "Rep B" = RepB,
         "Rep C" = RepC) %>% unite("sgRNA", c("gene", "sgRNA")) -> CRISPRn$readcount

read_and_merge <- function(obj, screen) {
  sprintf("dat/%s_cb2_sgrna.csv", screen) %>% read_csv() -> obj$cb2_sg
  sprintf("dat/%s_mageck_sgrna.txt", screen) %>% read.table(sep="\t", stringsAsFactors = F, header = T) %>% unite("id", c("Gene", "sgrna")) -> obj$mageck_sg
  obj$cb2_sg %>% left_join(obj$mageck_sg, by = c("sgRNA"="id")) -> obj$merged_sg
  obj
}

CRISPRn %>% read_and_merge("CRISPRn") -> CRISPRn

CRISPRn$merged_sg %>% filter(str_detect(sgRNA, "^NO_CURRENT")) %>% select(p_ts, p.twosided) %>% summary

CRISPRn$merged_sg %>% filter(str_detect(sgRNA, "^NO_CURRENT")) %>% 
  select(sgRNA, logFC, CB2=p_ts, MAGeCK=p.twosided) -> df_nontarget

CRISPRn$cb2_sg %>% select(sgRNA, P = p_ts) %>% mutate(label = str_detect(sgRNA, "^NO_CURRENT") %>% as.integer, method = "CB2") -> cb2
CRISPRn$mageck_sg %>% select(id, P = p.twosided) %>% mutate(label = str_detect(id, "^NO_CURRENT") %>% as.integer, method = "MAGeCK") -> mageck

df_fdr <- tibble()
for(eFDR in c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)) {
  df_fdr <- bind_rows(
    df_fdr,
    tibble(
      "P-value" = eFDR,
      "method" = "CB2",
      "Specificity (1-oFDR)" = 1 - sum(df_nontarget$CB2 < eFDR) / nrow(df_nontarget)
    )
  )    
  df_fdr <- bind_rows(
    df_fdr,
    tibble(
      "P-value" = eFDR,
      "method" = "MAGeCK",
      "Specificity (1-oFDR)" = 1 - sum(df_nontarget$MAGeCK < eFDR) / nrow(df_nontarget)
    )
  )    
}

as.character(df_fdr$`P-value`) -> df_fdr$`P-value`

(df_fdr %>% 
    rename(`Specificity`= `Specificity (1-oFDR)`) %>%
    ggplot(aes(x=`P-value`, y=`Specificity`)) +
    geom_point(aes(color=method, shape=method), size=2.5) +
    geom_line(aes(group=method, color = method), size=1) + 
    scale_color_npg() +
    ylim(0,1) +
    NULL -> ggp)




df_nontarget %>%
  gather(method, `-log10(P-value)`, -sgRNA, -logFC) %>%
  mutate(`-log10(P-value)` = -log10(`-log10(P-value)`),
         sig = `-log10(P-value)` >= -log10(0.01)) %>%
  ggplot(aes(x=logFC, y=`-log10(P-value)`)) +
  geom_point(alpha=0.25, size=0.7, aes(color=sig)) +
  geom_hline(yintercept =-log10(0.01), color = "blue", alpha=0.5) +
  theme(legend.position = "none") +
  facet_grid(~method) + 
  scale_color_manual(values = c("grey", "red")) +
  scale_y_continuous( limits = c(0, 50)) -> ggv


plot_grid(ggp, ggv, nrow=1, rel_widths = c(1,1), labels = "AUTO") -> pg
pg
save_plot("figures/fig-S12.pdf", pg, base_height = 4, base_width = 10)
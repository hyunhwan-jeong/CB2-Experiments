library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)
library(ggsci)

Gsk3 <- list()
read_delim("dat/Gsk3_readcount.txt", delim = "\t") %>%
  rename(sgRNA = gRNA,
         gene = Gene,
         "Presort 1" = Before_1,
         "Presort 2" = Before_2,
         "Presort 3" = Before_3,
         "Presort 4" = Before_4,
         "GFP+ 1" = Pos_1,
         "GFP+ 2" = Pos_2,
         "GFP+ 3" = Pos_3,
         "GFP+ 4" = Pos_4) -> Gsk3$readcount

Tnf <- list()
read_delim("dat/Tnf_readcount.tsv", delim = "\t")  %>% 
  select(-starts_with("Pre_")) %>%
  rename("TNF-Low 1" = Low_TNF_EX1,
         "TNF-Low 2" = Low_TNF_EX2,
         "TNF-Low 3" = Low_TNF_EX3,
         "TNF-High 1" = High_TNF_EX1,
         "TNF-High 2" = High_TNF_EX2,
         "TNF-High 3" = High_TNF_EX3) -> Tnf$readcount

summary(Tnf$readcount)

read_and_merge <- function(obj, screen) {
  sprintf("dat/%s_cb2_sgRNA.csv", screen) %>% read_csv() -> obj$cb2_sg
  sprintf("dat/%s_mageck_sgrna.txt", screen) %>% read.table(sep="\t", stringsAsFactors = F, header = T) -> obj$mageck_sg
  obj$cb2_sg %>% left_join(obj$mageck_sg, by = c("sgRNA"="sgrna")) -> obj$merged_sg
  obj
}

Tnf %>% read_and_merge("Tnf") -> Tnf

normalize <- function(obj) {
  for(i in 3:ncol(obj$readcount)) {
    s <- sum(obj$readcount[,i])
    obj$readcount[,i] <- obj$readcount[,i] / s * 10^6
  }
  obj
}

Tnf %>% normalize() -> Tnf

Tnf$merged_sg %>% filter(str_detect(sgRNA, "^NonTarget")) %>% select(p_value_twosided, p.twosided) %>% summary

Tnf$merged_sg %>% filter(str_detect(sgRNA, "^NonTarget")) %>% select(sgRNA, log2FC, p_value_twosided, p.twosided) -> df_nontarget
df_nontarget %>% select(sgRNA, CB2 = p_value_twosided, MAGeCK = p.twosided) %>%
  ggplot(aes(x=CB2, y=MAGeCK)) +
  geom_point()

Tnf$merged_sg %>% filter(str_detect(sgRNA, "^NonTarget")) %>% 
  select(sgRNA, CB2=p_value_twosided, MAGeCK=p.twosided)-> df_nontarget

Tnf$cb2_sg %>% select(sgRNA, P = p_value_twosided) %>% mutate(label = str_detect(sgRNA, "^NonTarget") %>% as.integer, method = "CB2") -> cb2
Tnf$mageck_sg %>% select(sgrna, P = p.twosided) %>% mutate(label = str_detect(sgrna, "^NonTarget") %>% as.integer, method = "MAGeCK") -> mageck

df_fdr <- tibble()
for(pval in c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)) {
  df_fdr <- bind_rows(
    df_fdr,
    tibble(
      "P-value" = pval,
      "CB2" = 1 - sum(df_nontarget$CB2 < pval) / nrow(df_nontarget),
      "MAGeCK" = 1 - sum(df_nontarget$MAGeCK < pval) / nrow(df_nontarget)
    )
  )    
}

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
    NULL -> ggp)
#theme(legend.position = "bottom") -> ggp)

Tnf$merged_sg %>% filter(str_detect(sgRNA, "^NonTarget")) %>% 
  select(sgRNA, log2FC, CB2=p_value_twosided, MAGeCK=p.twosided) -> df_nontarget

df_nontarget %>%
  gather(method, `-log10(P-value)`, -sgRNA, -log2FC) %>%
  mutate(`-log10(P-value)` = -log10(`-log10(P-value)`),
         sig = `-log10(P-value)` >= -log10(0.01)) %>%
  ggplot(aes(x=log2FC, y=`-log10(P-value)`)) +
  geom_point(alpha=0.25, size=0.7, aes(color=sig)) +
  geom_hline(yintercept =-log10(0.01), color = "blue", alpha=0.5) +
  theme(legend.position = "none") +
  facet_grid(~method) + 
  scale_color_manual(values = c("grey", "red")) +
  scale_y_continuous( limits = c(0, 50)) -> ggv

plot_grid(ggp, ggv, nrow=1, rel_widths = c(1,1), labels = "AUTO") -> pg
pg
save_plot("figures/fig4.pdf", pg, base_height = 4, base_width = 10)

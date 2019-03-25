library(tidyverse)
library(here)
library(cowplot)
library(pheatmap)
source(here("utils/heatmap_fdr.R"))

files <- Sys.glob("results/Evers/*/FDR/*")

e <- scan("data/Evers/essential-genes.txt", what ="character")
mk_df <- function(f) {
  read_csv(f) %>% 
    mutate(
      method = basename(f) %>% str_replace(".csv", ""),
      dataset = f %>% strsplit("/")) %>% 
    mutate(dataset=dataset[[1]][3])
}

all_df <- lapply(files, mk_df) %>% bind_rows %>% 
  mutate(method = ifelse(method == "HitSelect", "HiTSelect", method))

datasets <- c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112")
methods <- c("CB2", "ScreenBEAM", "PBNPA", "sgRSEA", "HiTSelect", "MAGeCK", "RIGER", "RSA", "PinAPL-Py") 

heatmap_fdr(all_df %>% select(-stat), methods, datasets, e)
# NOTE: If you export the plot as a PDF file, the file will not show some unicode character properly.
#save_plot("figures/fig-S7.png", last_plot(), base_height = 12) 
save_plot("figures/fig-S7.pdf", last_plot(), base_height = 12) 
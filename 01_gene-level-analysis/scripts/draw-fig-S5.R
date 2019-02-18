library(here)
source(here("utils/lineplot.R"))

files <- Sys.glob("results/Evers/*/FDR/*")

e <- scan("data/Evers/essential-genes.txt", what ="character")
mk_df <- function(f) {
  read_csv(f) %>% 
    mutate(
      method = basename(f) %>% str_replace(".csv", ""),
      dataset = f %>% strsplit("/")) %>% 
    mutate(dataset=dataset[[1]][3]) %>%
    mutate(essential = gene %in% e)
}

all_df <- lapply(files, mk_df) %>% bind_rows

datasets <- c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112")

lineplot_perf(all_df, col_ord = datasets)
save_plot(here("figures/fig-S5.pdf"), last_plot(), base_width = 8, base_height = 6)

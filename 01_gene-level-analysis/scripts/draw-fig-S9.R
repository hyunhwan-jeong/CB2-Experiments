library(here)
source(here("utils/lineplot.R"))

files <- Sys.glob("results/Sanson/*/FDR/*")

e <- scan("data/Sanson/essential-genes.txt", what ="character")
n <- scan("data/Sanson/non-essential-genes.txt", what ="character")
mk_df <- function(f) {
  read_csv(f) %>% 
    mutate(
      method = basename(f) %>% str_replace(".csv", ""),
      dataset = f %>% strsplit("/")) %>% 
    mutate(dataset=dataset[[1]][3]) %>%
    mutate(essential = gene %in% e) %>%
    filter(gene %in% c(e, n))
}

all_df <- lapply(files, mk_df) %>% bind_rows

datasets <- c("CRISPRn-A375", "CRISPRi-A375")

lineplot_perf(all_df, col_ord = datasets)
save_plot(here("figures/fig-S9.pdf"), last_plot(), base_width = 8, base_height = 6)

library(here)
library(tidyverse)
source(here("utils/load_data.R"))
source(here("utils/upset.R"))

draw_upset(load_data("Evers"), c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112"))
save_plot(filename = here("figures/UpSet-Evers.pdf"), last_plot(), base_height = 11, base_width = 10)
draw_upset(load_data("Sanson"), c("CRISPRn-A375", "CRISPRi-A375"))
save_plot(filename = here("figures/UpSet-Sanson.pdf"), last_plot(), base_height = 6, base_width = 10)




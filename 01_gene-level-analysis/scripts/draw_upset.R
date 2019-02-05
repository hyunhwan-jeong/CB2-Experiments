library(here)
library(tidyverse)
source(here("utils/load_data.R"))
source(here("utils/upset.R"))

draw_upset(load_data("Evers"), c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112"))
draw_upset(load_data("Sanson"), c("CRISPRn-A375", "CRISPRi-A375"))

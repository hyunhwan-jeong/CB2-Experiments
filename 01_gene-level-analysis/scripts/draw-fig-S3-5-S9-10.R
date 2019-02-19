library(tidyverse)
library(here)
source(here("scripts/draw_auc.R"))
study <- "Sanson"
methods <- c("CB2", "ScreenBEAM", "MAGeCK", "PBNPA", "RSA", "RIGER", "sgRSEA", "HitSelect")

fname_dict <- list( "CRISPRn-A375" = "figures/fig-S9.pdf", 
                    "CRISPRi-A375" = "figures/fig-S10.pdf",
                    "CRISPRn-RT112" = "figures/fig-S3.pdf",
                    "CRISPRn-UMUC3" = "figures/fig-S4.pdf", 
                    "CRISPRi-RT112" = "figures/fig-S5.pdf")

for(d in c("CRISPRn-A375", "CRISPRi-A375")) {
  filename <- here(fname_dict[[d]])
  if(d == "CRISPRn-A375") {
    title <- "Genome-wide CRISPRn library(Brunello)"
  } else {
    title <- "Genome-wide CRISPRi library(Dolcetto Set A)"
  }
  print(filename)
  draw_auc_figures(study, d, methods, title, filename)
}

study <- "Evers"
methods <- c("CB2", "ScreenBEAM", "MAGeCK", "PBNPA", "RSA", "RIGER", "sgRSEA", "HitSelect", "PinAPL-Py")
for(d in c("CRISPRn-RT112", "CRISPRn-UMUC3", "CRISPRi-RT112")) {
  filename <- here(fname_dict[[d]])
  title <- "Screen from Evers et al. ({d})" %>% glue
  print(filename)
  draw_auc_figures(study, d, methods, title, filename)
}

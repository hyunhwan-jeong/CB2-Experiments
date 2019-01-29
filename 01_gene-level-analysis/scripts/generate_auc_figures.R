source(here("scripts/draw_auc.R"))
study <- "Sanson"
methods <- c("CB2", "ScreenBEAM", "MAGeCK", "PBNPA", "RSA", "RIGER", "sgRSEA", "HitSelect")
for(d in c("CRISPRn-A375", "CRISPRi-A375")) {
  filename <- here("figures/AUC-{d}.pdf") %>%  glue
  if(d == "CRISRn-A375") {
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
  filename <- here("figures/AUC-{d}.pdf") %>%  glue
  title <- "Screen from Evers et al. ({d})" %>% glue
  print(filename)
  draw_auc_figures(study, d, methods, title, filename)
}

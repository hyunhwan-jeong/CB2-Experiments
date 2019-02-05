load_data <- function(study_id) {
  require(tidyr)
  require(magrittr)
  require(dplyr)
  require(glue)

  mk_df <- function(f) {
    read_csv(f) %>% 
      mutate(
        method = basename(f) %>% str_replace(".csv", ""),
        dataset = f %>% strsplit("/")) %>% 
      mutate(dataset=dataset[[1]][3]) %>%
      filter(gene %in% c(e, n)) %>% 
      mutate(essential = gene %in% e)
  }
  
  files <- Sys.glob("results/{study_id}/*/FDR/*" %>% glue)
  
  e <- scan("data/{study_id}/essential-genes.txt" %>% glue, what ="character")
  n <- scan("data/{study_id}/non-essential-genes.txt" %>% glue, what ="character")
  
  all_df <- lapply(files, mk_df) %>% bind_rows
  
}
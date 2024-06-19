
# getting source code -----------------------------------------------------
Source_patterns <- list(All = "[.]R$|[.]r$|[.]md$|[.]Rmd$|[.]rmd$")

Source_FullPath <- Source_patterns %>%
  purrr::imap(function(x, y) {
    c(
   
      file.path(getwd())
      
    ) %>%
      list.files(pattern = x,
                 ignore.case = TRUE, full.names = TRUE, recursive = TRUE) %>% 
      stringr::str_subset("BuildSupport/", negate = TRUE)
  }) 

Source_FullPath

sRead <- Source_FullPath %>%
  purrr::imap(function(x, y) {
    x %>%
      purrr::set_names(x) %>%
      purrr::imap(~ brio::readLines(.x))
  })


# dplyr::select(c("Feature" = "symbol", "ensembl_id")) 

sd_out<- sRead$All %>%
  purrr::imap(~ .x %>% stringr::str_subset(stringr::regex("tidyverse",
                                                          ignore_case = FALSE))) %>%
  purrr::compact()

sd_out

# open all files
sd_out %>% names() %>% file.edit()


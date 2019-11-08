get_reactome_pathways <- function(uniprot) {
    httr::GET(glue('https://reactome.org/ContentService/data/mapping/UniProt/{uniprot}/pathways?species=9606')) %>% {
        if (.$status != 200) {
            tibble(reactome = NA_character_)
        } else {
            httr::content(., as = 'parsed') %>%
                map(as_tibble) %>%
                bind_rows %>%
                select(stId, name) %>%
                rowwise() %>%
                mutate(name = paste(name, collapse = '; ')) %>%
                ungroup() %>%
                distinct() %>%
                rename(reactome = stId)
        }
    }
}


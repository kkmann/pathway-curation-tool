get_participating_proteins <- function(reactome) {
    if (is.na(reactome)) {
        return(tibble(tmp = NA_character_))
    } else {
        httr::GET(glue('https://reactome.org/ContentService/data/participants/{reactome}')) %>%
            httr::content(as = 'parsed') %>%
            map(as_tibble) %>%
            bind_rows() %>%
            transmute(
                tmp = map(refEntities, ~as_tibble(.))
            ) %>%
            unnest(tmp) %>%
            select(displayName) %>%
            filter(
                str_detect(displayName, '^(UniProt:)')
            ) %>%
            transmute(
                uniprotswissprot = str_extract(displayName, '(?<=:)[[:alnum:]]+(?=( |-))')
            ) %>%
            distinct()
    }
}

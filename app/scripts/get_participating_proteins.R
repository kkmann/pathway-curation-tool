get_participating_proteins <- function(reactome) {
    if (str_detect(reactome, '^UniProt:.+$')) {
        return(tibble(
            uniprotswissprot = str_extract(reactome, '(?<=^UniProt:)[0-9A-Z]+$')
        ))
    }
    if (str_detect(reactome, '^reactome:.+$')) {
        reactome <- str_extract(reactome, '(?<=reactome:)R-HSA-[0-9]+$')
        # prepare for timeouts
        query <- function() httr::GET(glue('https://reactome.org/ContentService/data/participants/{reactome}'))
        cntr  <- 0
        while (cntr >= 0 & cntr < 10) {
            tryCatch({
                    result <- query()
                    cntr <- -1
                }, error = function(e) {
                    cntr <- cntr + 1
                    warning(e, immediate = TRUE)
                    message('waiting 5 seconds, retrying')
                    Sys.sleep(5)
                }
            )
        }
        return(
            result %>%
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
        )
    }
    return(tibble(uniprotswissprot = NA_character_))
}

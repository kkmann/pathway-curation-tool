get_reactome_pathways <- function(uniprot) {
    # prepare for timeouts
    query <- function() httr::GET(glue('https://reactome.org/ContentService/data/mapping/UniProt/{uniprot}/pathways?species=9606'))
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
    result %>% {
        if (.$status != 200) {
            tibble(reactome = sprintf('UniProt:%s', uniprot))
        } else {
            httr::content(., as = 'parsed') %>%
                map(as_tibble) %>%
                bind_rows %>%
                select(stId, name) %>%
                rowwise() %>%
                mutate(name = paste(name, collapse = '; ')) %>%
                ungroup() %>%
                distinct() %>%
                rename(reactome = stId) %>%
                mutate(reactome = sprintf('reactome:%s', reactome))
        }
    }
}

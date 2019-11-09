get_reactome_interactions <- function(file) {
        read_tsv(
            file,
            col_types = cols_only(
                `# Interactor 1 uniprot id`    = col_character(),
                `Interactor 2 uniprot id`      = col_character()
            )
        ) %>%
        # only keep uniprot
        transmute(
            interactor_A = `# Interactor 1 uniprot id`,
            interactor_B = `Interactor 2 uniprot id`
        ) %>%
        # throw out everything not representing a true protein-protein edge
        filter(
            str_detect(interactor_A, '^uniprotkb:'),
            str_detect(interactor_B, '^uniprotkb:')
        ) %>%
        distinct() %>%
        # order
        mutate(
            tmp = map2(interactor_A, interactor_B, ~sort(c(.x, .y)))
        ) %>%
        filter(map_int(tmp, length) == 2) %>%
        mutate(
            interactor_A = map_chr(tmp, ~.[[1]]),
            interactor_B = map_chr(tmp, ~.[[2]])
        ) %>%
        select(-tmp) %>%
        distinct() %>%
        mutate(
            id    = sprintf("reactome:%09i", row_number()),
            score = 1
        ) %>%
        select(id, score, everything())
}

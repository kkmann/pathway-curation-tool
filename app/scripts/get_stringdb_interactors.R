get_stringdb_interactors <- function(file, tbl_ensembl, min_combined_score) {
    read_delim(
        gzfile(file),
        delim = ' ',
        col_types = cols(
            protein1 = col_character(),
            protein2 = col_character(),
            neighborhood = col_double(),
            neighborhood_transferred = col_double(),
            fusion = col_double(),
            cooccurence = col_double(),
            homology = col_double(),
            coexpression = col_double(),
            coexpression_transferred = col_double(),
            experiments = col_double(),
            experiments_transferred = col_double(),
            database = col_double(),
            database_transferred = col_double(),
            textmining = col_double(),
            textmining_transferred = col_double(),
            combined_score = col_double()
        )
    ) %>%
    filter(
        experiments > 0 | database > 0, combined_score > min_combined_score
    ) %>%
    transmute(
        interactor_A = protein1,
        interactor_B = protein2
    ) %>%
    distinct %>%
    filter(
        interactor_A != interactor_B
    ) %>%
    mutate_all(
        ~str_extract(., '(?<=^9606.)ENSP[0-9]{11}$')
    ) %>%
    mutate(
        id = sprintf("stringdb:%09i", row_number())
    ) %>%
    pivot_longer(-id, names_to = 'interactor', values_to = 'ensembl_peptide_id') %>%
    mutate(
        interactor = str_extract(interactor, '(?<=^interactor_)[A|B]$')
    ) %>%
    left_join(
        tbl_ensembl %>% select(uniprotswissprot, ensembl_peptide_id),
        by = 'ensembl_peptide_id'
    ) %>%
    filter(
        !is.na(uniprotswissprot)
    ) %>%
    select(., -ensembl_peptide_id) %>%
    group_by(id) %>% {
        # resolve 1-p for peptide id to uniprot
        bind_rows(
            filter(
                .,
                length(unique(uniprotswissprot)) != 2 |
                length(unique(interactor)) != 2
            ) %>%
            nest(data = uniprotswissprot) %>%
            pivot_wider(names_from = interactor, values_from = data) %>%
            unnest(A) %>%
            rename(A = uniprotswissprot) %>%
            unnest(B) %>%
            rename(B = uniprotswissprot) %>%
            filter(A != B) %>%
            mutate(
                id2 = sprintf("%s.%i", id, row_number())
            ) %>%
            ungroup() %>%
            transmute(
                id = id2, A, B
            ) %>%
            pivot_longer(-id, names_to = 'interactor', values_to = 'uniprotswissprot'),

            filter(
                .,
                length(unique(uniprotswissprot)) == 2,
                length(unique(interactor)) == 2
            ) %>%
            ungroup()
        )
    } %>%
    pivot_wider(names_from = interactor, values_from = uniprotswissprot) %>%
    distinct(A, B, .keep_all = TRUE) %>%
    pivot_longer(-id, names_to = 'interactor', values_to = 'uniprotswissprot')
}

get_stringdb_interactions <- function(file, tbl_ensembl) {
    read_delim(
        gzfile(file),
        delim = ' ',
        col_types = cols_only(
            protein1       = col_character(),
            protein2       = col_character(),
            experiments    = col_double(),
            database       = col_double(),
            combined_score = col_double()
        )
    ) %>%
    filter(
        experiments > 0 | database > 0
    ) %>%
    transmute(
        interactor_A = protein1,
        interactor_B = protein2,
        score        = combined_score / 1000
    ) %>%
    filter(
        interactor_A != interactor_B
    ) %>%
    mutate_at(
        vars(starts_with('interactor')),
        ~str_extract(., '(?<=^9606.)ENSP[0-9]{11}$')
    ) %>%
    mutate(
        id = sprintf("stringdb:%09i", row_number())
    ) %>%
    left_join(
        tbl_ensembl %>% select(uniprotswissprot, ensembl_peptide_id),
        by = c(interactor_A = 'ensembl_peptide_id')
    ) %>%
    rename(uniprotswissprot_A = uniprotswissprot) %>%
    left_join(
        tbl_ensembl %>% select(uniprotswissprot, ensembl_peptide_id),
        by = c(interactor_B = 'ensembl_peptide_id')
    ) %>%
    rename(uniprotswissprot_B = uniprotswissprot) %>%
    transmute(
        id,
        score,
        interactor_A = uniprotswissprot_A,
        interactor_B = uniprotswissprot_B
    ) %>%
    mutate(
        tmp = map2(interactor_A, interactor_B, ~sort(c(.x, .y)))
    ) %>%
    filter(map_int(tmp, length) == 2) %>%
    mutate(
        interactor_A = map_chr(tmp, ~.[[1]]),
        interactor_B = map_chr(tmp, ~.[[2]])
    ) %>%
    select(-tmp) %>%
    distinct(interactor_A, interactor_B, .keep_all = TRUE)
}

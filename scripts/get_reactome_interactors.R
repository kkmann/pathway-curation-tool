get_reactome_interactors <- function(file) {
        read_tsv(
            file,
            col_types = cols(
                `# Interactor 1 uniprot id`    = col_character(),
                `Interactor 1 Ensembl gene id` = col_character(),
                `Interactor 1 Entrez Gene id`  = col_character(),
                `Interactor 2 uniprot id`      = col_character(),
                `Interactor 2 Ensembl gene id` = col_character(),
                `Interactor 2 Entrez Gene id`  = col_character(),
                `Interaction type`             = col_character(),
                `Interaction context`          = col_character(),
                `Pubmed references`            = col_character()
            )
        ) %>%
        # only keep more xhaustive ensembl gene id info
        transmute(
            interactor_A = `# Interactor 1 uniprot id`,
            interactor_B = `Interactor 2 uniprot id`,
            type         = `Interaction type`,
            context      = `Interaction context`,
            pubmed       = `Pubmed references`
        ) %>%
        # throw out everything not representing a true edge
        filter(
            !is.na(interactor_A),
            interactor_A != '-',
            !is.na(interactor_B),
            interactor_B != '-',
        ) %>%
        distinct() %>%
        # convert to interactor long format (one row per interactor)
        mutate(
            id = sprintf("reactome:%09i", row_number())
        ) %>%
        pivot_longer(
            c(interactor_A, interactor_B),
            names_to  = 'interactor',
            values_to = 'uniprotswissprot'
        ) %>%
        mutate(
            uniprotswissprot = str_extract(uniprotswissprot, '(?<=^uniprotkb:)[0-9A-Z]{6}$'),
            interactor      = str_extract(interactor, '(?<=^interactor_)[A|B]$')
        ) %>%
        select(id, interactor, uniprotswissprot, everything()) %>%
        filter(
            !is.na(uniprotswissprot)
        ) %>%
        distinct(id, uniprotswissprot, .keep_all = TRUE) %>%
        group_by(id) %>%
        filter(
            n() == 2
        ) %>%
        ungroup() %>%
        arrange(id, interactor)
}

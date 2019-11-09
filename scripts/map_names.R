map_proteins_to_genes <- function(proteins, tbl_ensembl) {
    tibble(
        uniprotswissprot = proteins
    ) %>%
    left_join(
        select(tbl_ensembl, 'uniprotswissprot', 'ensembl_gene_id'),
        by = 'uniprotswissprot'
    ) %>%
    pull(ensembl_gene_id) %>%
    unique
}

map_ensembl_to_name <- function(ensembl) {
    tibble(
        ensembl_gene_id = ensembl
    ) %>%
    left_join(
        select(tbls$ensembl, 'external_gene_name', 'ensembl_gene_id'),
        by = 'ensembl_gene_id'
    ) %>%
    pull(external_gene_name) %>%
    unique
}

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

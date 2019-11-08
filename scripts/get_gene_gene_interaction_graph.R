get_gene_gene_interaction_graph <- function(proteins, tbl_gene_gene_interactions, tbl_ensembl) {
    genes <- tibble(
            uniprotswissprot = proteins
        ) %>%
        left_join(
            tbl_ensembl, by = 'uniprotswissprot'
        ) %>%
        filter(complete.cases(.)) %>%
        pull(ensembl_gene_id) %>%
        unique()

    edges <- tbl_gene_gene_interactions %>%
        filter(gene_A %in% genes & gene_B %in% genes) %>%
        as.matrix() %>%
        t() %>%
        as.vector() %>%
        as.character()

    isolates <- genes[!genes %in% unique(edges)]

    igraph::make_graph(
        edges    = edges,
        isolates = isolates,
        directed = FALSE
    )
}

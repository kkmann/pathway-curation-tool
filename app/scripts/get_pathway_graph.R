get_pathway_graph <- function(gene_data, tbl_gene_gene_interactions, tbl_ensembl) {
    genes         <- gene_data$ensembl_gene_ids_pruned %>% unlist
    seed_genes    <- gene_data$ensembl_gene_ids_seed %>% unlist
    covered_genes <- gene_data$ensembl_gene_ids_pruned_covered %>% unlist
    tbl_interactions <- tbl_gene_gene_interactions %>%
        filter(gene_A %in% genes & gene_B %in% genes)
    edges <- tbl_interactions %>%
        as.matrix() %>%
        t() %>%
        as.character()
    isolates <- genes[!(genes %in% unique(edges))]
    gr <- igraph::make_graph(
        edges    = edges,
        isolates = isolates,
        directed = FALSE
    )
    # add external gene names
    external_gene_names <- tibble(
        ensembl_gene_id = igraph::vertex_attr(gr, name = 'name', igraph::V(gr))
    ) %>%
    left_join(
        select(tbl_ensembl, ensembl_gene_id, external_gene_name) %>% distinct(),
        by = 'ensembl_gene_id'
    )  %>%
    pull(external_gene_name)
    gr <- igraph::set_vertex_attr(gr, 'external_gene_name', value = external_gene_names)
    # add seed gene information
    is_seed_gene <- igraph::vertex_attr(gr, name = 'name', igraph::V(gr)) %in% seed_genes
    gr <- igraph::set_vertex_attr(gr, 'seed_gene', value = is_seed_gene)
    # add coverage information
    is_covered <- igraph::vertex_attr(gr, name = 'name', igraph::V(gr)) %in% covered_genes
    gr <- igraph::set_vertex_attr(gr, 'covered', value = is_covered)
    return(gr)
}

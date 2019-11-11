get_pathway_graph <- function(genes, seed_genes, tbl_interactions) {
    edges <- tbl_interactions %>%
        filter(gene_A %in% genes$ensembl_gene_id & gene_B %in% genes$ensembl_gene_id) %>%
        select(gene_A, gene_B) %>%
        as.matrix() %>%
        t() %>%
        as.character()
    isolates <- genes$ensembl_gene_id[!(genes$ensembl_gene_id %in% unique(edges))]
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
            select(tbls$ensembl, ensembl_gene_id, external_gene_name) %>% distinct(),
            by = 'ensembl_gene_id'
        )  %>%
        pull(external_gene_name)
    gr <- igraph::set_vertex_attr(gr, 'external_gene_name', value = external_gene_names)
    # add seed gene information
    is_seed_gene <- igraph::vertex_attr(gr, name = 'name', igraph::V(gr)) %in% seed_genes$ensembl_gene_id
    gr <- igraph::set_vertex_attr(gr, 'seed_gene', value = is_seed_gene)
    # add coverage information
    is_covered <- tibble(
            ensembl_gene_id = igraph::vertex_attr(gr, name = 'name', igraph::V(gr))
        ) %>% left_join(
            genes, by = 'ensembl_gene_id'
        ) %>%
        pull(covered)
    gr <- igraph::set_vertex_attr(gr, 'covered', value = is_covered)
    return(gr)
}

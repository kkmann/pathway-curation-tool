prune_genelist <- function(genes, seed_genes, k, tbl_interactions) {
    tbl_pruned_interactions <- tbl_interactions %>%
        select(gene_A, gene_B) %>%
        filter(
            gene_A %in% genes,
            gene_B %in% genes
        )
    net <- igraph::graph(
        tbl_pruned_interactions %>%
            as.matrix() %>%
            t() %>%
            as.character(),
        isolates = genes[
                !(genes %in% tbl_pruned_interactions$gene_A) &
                !(genes %in% tbl_pruned_interactions$gene_B)
            ],
        directed = FALSE
    )
    d <- igraph::distances(
        net,
        igraph::V(net),
        to = seed_genes
    )
    genes %in% rownames(d)[apply(d, 1, min) <= k]
}

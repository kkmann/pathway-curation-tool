prune <- function(uniprotswissprot, seed_proteins, k, tbl_interactions) {
    tbl_pruned_interactions <- tbl_interactions %>%
        select(interactor_A, interactor_B) %>%
        filter(
            interactor_A %in% uniprotswissprot,
            interactor_B %in% uniprotswissprot
        )
    net <- igraph::graph(
        tbl_pruned_interactions %>%
            as.matrix() %>%
            t() %>%
            as.vector(),
        isolates = uniprotswissprot[
                !(uniprotswissprot %in% tbl_pruned_interactions$interactor_A) &
                !(uniprotswissprot %in% tbl_pruned_interactions$interactor_B)
            ],
        directed = FALSE
    )
    d <- igraph::distances(
        net,
        igraph::V(net),
        to = seed_proteins
    )
    uniprotswissprot %in% rownames(d)[apply(d, 1, min) <= k]
}

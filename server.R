library(shiny)
library(glue)
library(tidyverse)

map(
    list.files('scripts', pattern = '.+\\.R$', full.names = TRUE),
    ~source(., verbose = FALSE)
)

dir.create('data', showWarnings = FALSE)
download_zenodo <- function(filename, doi = '3538474') {
    if (!file.exists(glue('data/{filename}')))
        download.file(
            glue('https://zenodo.org/record/{doi}/files/{filename}?download=1'),
            glue('data/{filename}')
        )
}
c('tbl_grch37.ensembl.org.rds', 'tbl_interactions_prot_prot.rds', 'tbl_interactions_gene_gene.rds') %>%
    map(download_zenodo)

# tbl_ensembl <- biomaRt::getBM(
#         attributes = c(
#             'external_gene_name',
#             'ensembl_gene_id',
#             'gene_biotype',
#             'uniprotswissprot',
#             'ensembl_peptide_id'
#         ),
#         mart = biomaRt::useMart(
#             "ENSEMBL_MART_ENSEMBL",
#             host    = "grch37.ensembl.org",
#             path    = "/biomart/martservice",
#             dataset = "hsapiens_gene_ensembl"
#         )
#     ) %>%
#     as_tibble %>%
#     arrange(gene_biotype, external_gene_name) %>%
#     # we are using ensembl_pepdide_id and ensembl_gene_id as primary keys, need both!
#     mutate_all(
#         ~ifelse(. == '', NA_character_, .)
#     )
# write_rds(tbl_ensembl, 'data/tbl_grch37.ensembl.org.rds', compress = 'gz')


# bind_rows(
#     get_reactome_interactions(
#         file = 'data/reactome.homo_sapiens.interactions.tab-delimited.txt'
#     ),
#     get_stringdb_interactions(
#         file               = 'data/9606.protein.links.full.v11.0.txt.gz',
#         tbl_ensembl        = tbls$ensembl
#     )
# ) %>%
# write_rds('data/tbl_interactions_prot_prot.rds', compress = 'gz')


tbls                            <- list()
tbls$ensembl                    <- read_rds('data/tbl_grch37.ensembl.org.rds')
tbls$interactions_prot_prot_all <- read_rds('data/tbl_interactions_prot_prot.rds')
tbls$interactions_gene_gene_all <- read_rds('data/tbl_interactions_gene_gene.rds')

renderDataTable_custom <- function(tbl) {
    DT::renderDataTable(
        {tbl},
        options = list(
            lengthMenu = c(5, 15, 50, 10000),
            pageLength = 10000
        ),
        escape  = FALSE)
}



server <- function(input, output) {

    tbls$interactions_gene_gene <- reactive({
        tbls$interactions_gene_gene_all %>%
            filter(score >= input$minScore)
    })

    tbls$seed_genes <- tibble(
        external_gene_name = character(0),
        ensembl_gene_id    = character(0),
        gene_biotype       = character(0),
        proteins           = list()
    )
    updateSeedGenesTable <- function() {
        output$tblSeedGenes <<- renderDataTable_custom(
            tbls$seed_genes %>%
                mutate(
                    proteins = map_chr(proteins, ~paste(., collapse = '; '))
                )
        )
    }
    updateSeedGenesTable()
    observeEvent(input$uploadSeedGenesFile, {
        tbls$seed_genes <<- read_csv(
                input$uploadSeedGenesFile$datapath[1],
                col_types = cols_only(
                    ensembl_gene_id = col_character()
                )
            ) %>%
            left_join(
                tbls$ensembl %>%
                    select(external_gene_name, ensembl_gene_id, gene_biotype, uniprotswissprot) %>%
                    filter(complete.cases(.)),
                by = 'ensembl_gene_id'
            ) %>%
            distinct() %>%
            arrange(external_gene_name) %>%
            group_by(ensembl_gene_id) %>%
            nest(proteins = uniprotswissprot) %>%
            mutate(proteins = map(proteins, ~.$uniprotswissprot))
        updateSeedGenesTable()
    })



    tbls$candidate_pathways <- tibble(
        pathway   = character(0),
        names     = character(0),
        seed_data = list(),
        proteins  = list(),
        genes     = list(),
        selected  = logical()
    )
    get_link <- function(pathway) {
        id <- str_extract(pathway, '(?<=:).+')
        if (str_detect(pathway, '^reactome:')) {
            link <- glue("https://reactome.org/content/detail/{id}")
        }
        if (str_detect(pathway, '^UniProt:')) {
            link <- glue("https://www.uniprot.org/uniprot/{id}")
        }
        glue("<a href={link}>{pathway}</a>")
    }
    queryPathways <- function() {
        withProgress(
            message = 'querying reactome pathways',
            value   = 1,
            max     = tbls$seed_genes %>%
                filter(!is.na(ensembl_gene_id)) %>%
                nrow(),
            tbl_tmp <- tbls$seed_genes %>%
                mutate(
                    pathways = map(
                        proteins,
                        function(x) {
                            incProgress(1, detail = sprintf("UniProt: %s", x))
                            map(x, get_reactome_pathways) %>%
                                bind_rows() %>%
                                distinct()
                        }
                    )
                )
        )
        withProgress(
            message = 'querying participating proteins',
            value   = 1,
            max     = nrow(tbl_tmp),
            tbls$candidate_pathways <<- tbl_tmp %>%
                unnest(pathways) %>%
                group_by(pathway) %>%
                nest(
                    seed_data = c(ensembl_gene_id, external_gene_name, gene_biotype, proteins),
                    names     = name
                ) %>%
                ungroup() %>%
                mutate(
                    seed_data = map(seed_data, ~distinct(., ensembl_gene_id, .keep_all = TRUE)),
                    names     = map(names, ~sort(unique(.$name))),
                    proteins = map(
                        pathway,
                        function(x) {
                            incProgress(1, detail = x)
                            get_participating_proteins(x)
                        }
                    ),
                    genes = map(proteins, function(x) {
                        tibble(
                            ensembl_gene_id = map_proteins_to_genes(x, tbls$ensembl)
                        ) %>%
                            left_join(
                                select(tbls$ensembl, ensembl_gene_id, external_gene_name) %>% distinct(),
                                by = 'ensembl_gene_id'
                            )
                    })
                ) %>%
                mutate(selected = FALSE) %>%
                select(pathway, names, seed_data, proteins, genes, selected)
        )
        updateCandidatePathwaysTables()
    }
    observeEvent(input$queryPathways, {queryPathways()})
    updateCandidatePathwaysTables <- function() {
        output$tbl_candidate_pathways <<- renderDataTable_custom(
            tbls$candidate_pathways %>%
                transmute(
                    selected,
                    pathway       = map_chr(pathway, get_link),
                    names         = map_chr(names, ~paste(., collapse = '; ')),
                    seed_genes    = map_chr(seed_data, ~paste(.$external_gene_name, collapse = '; ')),
                    seed_proteins = map_chr(seed_data, ~paste(.$proteins, collapse = '; ')),
                    n_proteins    = map_int(proteins, length),
                    n_genes       = map_int(genes, ~length(.$ensembl_gene_id))
                )
        )
        output$tbl_candidate_pathway_clusters <- renderTable({
            if (any(tbls$candidate_pathways$selected)) {
                tbls$candidate_pathways %>%
                    filter(selected) %>%
                    summarise(
                        pathways     = as.integer(n()),
                        genes        = bind_rows(genes) %>% pull(ensembl_gene_id) %>% unique %>% length %>% as.integer(),
                        `seed genes` = bind_rows(seed_data) %>% pull(external_gene_name) %>% unique %>% sort %>% paste(collapse = '; ')
                    )
            } else {
                tibble(
                    pathways     = 0L,
                    genes        = 0L,
                    `seed genes` = 'none'
                )
            }
        })
        refresh_pruning()
    }
    observeEvent(input$assign_cluster, {
        tbls$candidate_pathways$selected[input$tbl_candidate_pathways_rows_selected] <<-
            !tbls$candidate_pathways$selected[input$tbl_candidate_pathways_rows_selected]
        updateCandidatePathwaysTables()
    })
    output$downloadClusterAssignment <- downloadHandler(
        'cluster-assignments.csv',
        function(file) {
            tbls$candidate_pathways %>%
                filter(selected) %>%
                transmute(
                    pathway,
                    names         = map_chr(names, ~paste(., collapse = '; ')),
                    seed_genes    = map_chr(seed_data, ~paste(.$external_gene_name, collapse = '; ')),
                    seed_proteins = map_chr(seed_data, ~paste(.$proteins, collapse = '; ')),
                    n_proteins    = map_int(proteins, length),
                    n_genes       = map_int(genes, ~length(.$ensembl_gene_id))
                ) %>%
                write_csv(file)
        }
    )
    observeEvent(input$uploadClusterAssignmentFile, {
        tbls$candidate_pathways <<- tbls$candidate_pathways %>%
            select(-selected) %>%
            left_join(
                read_csv(
                    input$uploadClusterAssignmentFile$datapath[1],
                    col_types = cols_only(pathway = col_character()),
                ) %>%
                mutate(selected = TRUE),
                by = 'pathway'
            ) %>%
            select(selected, everything()) %>%
            mutate(selected = ifelse(is.na(selected), FALSE, selected))
        updateCandidatePathwaysTables()
    })

    tbls$coverage <- reactive({
        if (is.null(input$uploadGeneCoverageFile$datapath)) {
            res <- NULL
        } else {
            res <- read_csv(
                    input$uploadGeneCoverageFile$datapath[1],
                    col_types = cols(
                        ensembl_gene_id = col_character()
                    )
                ) %>%
                distinct()
        }
        return(res)
    })
    refresh_pruning <- function() {
        if (nrow(tbls$candidate_pathways %>% filter(selected)) == 0) {
            tbls$pathway_cluster <- tibble(
                pathways      = character(0),
                names         = character(0),
                seed_genes    = character(0),
                genes         = character(0),
                genes_pruned  = character(0),
                igraph_pruned = character(0)
            )
        } else {
            k <- input$pruning_distance
            covered <- function(ensembl_gene_ids) {
                if (is.null(tbls$coverage())) {
                    return(ensembl_gene_ids == ensembl_gene_ids)
                } else {
                    return(ensembl_gene_ids %in% tbls$coverage()$ensembl_gene_id)
                }
            }
            tbls$pathway_cluster <<- tbls$candidate_pathways %>%
                filter(selected) %>%
                select(-selected) %>%
                summarize(
                    pathways   = list(pathway),
                    names      = list(unlist(names) %>% unique %>% sort),
                    seed_genes = list(seed_data %>% bind_rows %>% pull(ensembl_gene_id) %>% unique %>% sort),
                    genes      = list(genes %>% bind_rows %>% distinct %>% pull(ensembl_gene_id))
                ) %>%
                mutate(
                    genes_pruned = list(
                        tibble(
                            ensembl_gene_id = genes[[1]][prune_genelist(genes[[1]], seed_genes[[1]], k, tbl_interactions = tbls$interactions_gene_gene())],
                            covered         = covered(ensembl_gene_id)
                        )
                    ),
                    seed_genes = list(
                        tibble(
                            ensembl_gene_id = seed_genes[[1]],
                            covered         = covered(ensembl_gene_id)
                        )
                    ),
                    genes = list(
                        tibble(
                            ensembl_gene_id = genes[[1]],
                            covered         = covered(ensembl_gene_id)
                        )
                    ),
                    igraph_pruned = list(get_pathway_graph(
                        genes_pruned[[1]],
                        seed_genes[[1]],
                        tbl_interactions = tbls$interactions_gene_gene(),
                        tbl_ensembl =  tbls$ensembl
                    ))
                )
        }
        output$tbl_pruned_pw_cluster <- DT::renderDataTable({
                if (nrow(tbls$candidate_pathways %>% filter(selected)) == 0) {
                    tbls$pathway_cluster
                } else {
                    tbls$pathway_cluster %>%
                        transmute(
                            pathways             = length(pathways[[1]]),
                            genes                = sprintf('%i (%5.1f%% covered)', nrow(genes[[1]]), 100*nrow(genes[[1]] %>% filter(covered))/nrow(genes[[1]])),
                            pruned               = sprintf('%i (%5.1f%% covered)', nrow(genes_pruned[[1]]), 100*nrow(genes_pruned[[1]] %>% filter(covered))/nrow(genes_pruned[[1]])),
                            `seed genes covered` = paste(
                                seed_genes[[1]] %>%
                                    filter(covered) %>%
                                    pull(ensembl_gene_id) %>%
                                    map_ensembl_to_name(tbls$ensembl),
                                collapse = '; '
                            ),
                            `seed genes not covered` = paste(
                                seed_genes[[1]] %>%
                                    filter(!covered) %>%
                                    pull(ensembl_gene_id) %>%
                                    map_ensembl_to_name(tbls$ensembl),
                                collapse = '; '
                            )
                        )
                }
            },
            options = list(
                lengthMenu = c(10, 25, 100),
                pageLength = 25
            )
        )
    }
    observeEvent(input$pruning_distance, {refresh_pruning()})
    observeEvent(input$intersectGREX, {refresh_pruning()})

    pathway_cluster_plotter <- reactive({
        if (is.null(tbls$pathway_cluster$igraph_pruned[[1]])) return(NULL)
        withProgress({
                gr  <- tbls$pathway_cluster$igraph_pruned[[1]]
                ggr <- tidygraph::as_tbl_graph(gr) %>%
                    tidygraph::activate(nodes) %>%
                    mutate(
                        alpha = ifelse(covered, '1', '0.33'),
                        size  = ifelse(seed_gene, '6', '1')
                    ) %>%
                    {if (input$plotNonCovered) . else filter(., covered)} %>%
                    mutate(
                        label = if (input$ensemblId == 'external_gene_name') external_gene_name else name
                    )
                p <- ggraph::ggraph(ggr) +
                    ggraph::geom_node_point(aes(alpha = alpha, size = size), color = 'black') +
                    ggraph::geom_edge_link(alpha = .05) +
                    ggraph::geom_node_text(aes(alpha = alpha, label = label), size = input$labelSize, colour = 'white', repel = TRUE) +
                    scale_alpha_manual(values = c('0.33' = 0.33, '1' = 1), guide = 'none') +
                    scale_size_manual(values = c('6' = 3*input$pointSize, '1' = input$pointSize), guide = 'none') +
                    ggraph::theme_graph(background = '#dddddd')
                incProgress(1, message = 'done')
                return(p)
            },
            value = 0, message = 'plotting...'
        )
    })
    output$plotPathwayCluster <- renderPlot({
            pathway_cluster_plotter()
        },
        height = function() 72*input$plotHeight,
        width  = function() 72*input$plotWidth,
        res    = 72
    )
    output$downloadPlot <- downloadHandler(
        'pathway_cluster_graph.pdf',
        function(file) ggsave(filename = file, pathway_cluster_plotter(), device = 'pdf', height = input$plotHeight, width  = input$plotWidth)
    )

    output$downloadGraph <- downloadHandler(
        'pathway_cluster_graph.rds',
        function(file) write_rds(tbls$pathway_cluster$igraph_pruned[[1]], path = file, compress = 'gz')
    )

    updateCandidatePathwaysTables()

}

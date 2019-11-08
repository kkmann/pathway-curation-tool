library(shiny)
library(glue)
library(tidyverse)

map(
    list.files('scripts', pattern = '.+\\.R$', full.names = TRUE),
    ~source(., verbose = FALSE)
)

dir.create('data', showWarnings = FALSE)
download_zenodo <- function(filename, doi = '3532192') {
    if (!file.exists(glue('data/{filename}')))
        download.file(
            glue('https://zenodo.org/record/{doi}/files/{filename}?download=1'),
            glue('data/{filename}')
        )
}
c('tbl_ensembl.rds', 'tbl_gene_gene_interactions.rds', 'tbl_interactions.rds') %>%
    map(download_zenodo)

# ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# write_rds(tbls$ensembl, 'data/tbl_ensembl.rds')
# tbls$ensembl <- biomaRt::getBM(
#         attributes = c(
#             'external_gene_name',
#             'ensembl_gene_id',
#             'gene_biotype',
#             'uniprotswissprot',
#             'ensembl_peptide_id'
#         ),
#         mart       = ensembl
#     ) %>%
#     as_tibble %>%
#     arrange(gene_biotype, external_gene_name) %>%
#     # we are using ensembl_pepdide_id and ensembl_gene_id as primary keys, need both!
#     filter(
#         ensembl_peptide_id != '',
#         ensembl_gene_id != ''
#     ) %>%
#     mutate_all(
#         ~ifelse(. == '', NA_character_, .)
#     )

# tbls$interactions <- bind_rows(
#         get_reactome_interactors(
#             file = 'data/reactome.homo_sapiens.interactions.tab-delimited.txt'
#         ),
#         get_stringdb_interactors(
#             file               = 'data/9606.protein.links.full.v11.0.txt.gz',
#             tbl_ensembl        = tbls$ensembl,
#             min_combined_score = 750
#         )
#     ) %>%
#     pivot_wider(
#         names_from   = interactor,
#         names_prefix = 'interactor_',
#         values_from  = uniprotswissprot
#     ) %>%
#     distinct(interactor_A, interactor_B, .keep_all = TRUE) %>%
#     select(id, interactor_A, interactor_B, everything())
# write_rds(tbls$interactions, 'data/tbl_interactions.rds')

# tbls$gene_gene_interactions <- tbls$interactions %>%
#     select(interactor_A, interactor_B) %>%
#     left_join(
#         select(tbls$ensembl, ensembl_gene_id, uniprotswissprot) %>%
#             filter(complete.cases(.)) %>%
#             distinct(),
#         by = c(interactor_A = 'uniprotswissprot')
#     ) %>%
#     filter(complete.cases(.)) %>%
#     rename(gene_A = ensembl_gene_id) %>%
#     left_join(
#         select(tbls$ensembl, ensembl_gene_id, uniprotswissprot) %>%
#             filter(complete.cases(.)) %>%
#             distinct(),
#         by = c(interactor_B = 'uniprotswissprot')
#     ) %>%
#     rename(gene_B = ensembl_gene_id) %>%
#     select(gene_A, gene_B) %>%
#     distinct() %>%
#     filter(gene_A != gene_B) %>%
#     rowwise() %>%
#     mutate(
#         gene_A1 = sort(c(gene_A, gene_B))[1],
#         gene_B1 = sort(c(gene_A, gene_B))[2]
#     ) %>%
#     ungroup() %>%
#     transmute(
#         gene_A = gene_A1,
#         gene_B = gene_B1
#     ) %>%
#     distinct()
# write_rds(tbls$gene_gene_interactions, 'data/tbl_gene_gene_interactions.rds')

tbls                        <- list()
tbls$ensembl                <- read_rds('data/tbl_ensembl.rds')
tbls$interactions           <- read_rds('data/tbl_interactions.rds')
tbls$gene_gene_interactions <- read_rds('data/tbl_gene_gene_interactions.rds')



server <- function(input, output) {

    tbls$candidate_pathways_reactome <- c()
    tbls$seed_genes <- reactive({
        if (is.null(input$uploadSeedGenesFile$datapath)) {
            res <- NULL
        } else {
            res <- input$uploadSeedGenesFile$datapath[1] %>%
                read_delim(
                    delim   = ' ',
                    trim_ws = TRUE,
                    col_types = cols(
                        gene_name = col_character(),
                        reference = col_character()
                    )
                ) %>%
                distinct() %>%
                select(-reference) %>%
                left_join(
                    tbls$ensembl %>%
                        select(external_gene_name, ensembl_gene_id, gene_biotype, uniprotswissprot) %>%
                        filter(complete.cases(.)),
                    by = c(gene_name = 'external_gene_name')
                ) %>%
                distinct() %>%
                arrange(gene_name) %>%
                rename(external_gene_name = gene_name)
        }
        return(res)
    })
    observeEvent(input$uploadSeedGenesFile, {
        output$tblSeedGenes <- DT::renderDataTable({
                tbls$seed_genes() %>% filter(!is.na(ensembl_gene_id))
            },
            options = list(
                lengthMenu   = c(5, 10, 25, 100),
                pageLength   = 25
            )
        )
        output$tblSeedGenesNotFound <- DT::renderDataTable({
                tbls$seed_genes() %>% filter(is.na(ensembl_gene_id))
            },
            options = list(
                lengthMenu   = c(5, 10, 25, 100),
                pageLength   = 5
            )
        )
    })
    output$downloadSeedGenes <- downloadHandler(
        'seed-genes.txt',
        function(file) {
            read_delim(
                input$uploadSeedGenesFile$datapath[1],
                delim   = ' ',
                trim_ws = TRUE,
                col_types = cols(
                    gene_name        = col_character(),
                    source_reference = col_character()
                )
            ) %>%
            write_delim(
                path    = file,
                delim   = ' '
            )
        }
    )

    update_candidate_pathways_tables <- function() {
        output$tbl_candidate_pathways <- DT::renderDataTable({
                tbls$candidate_pathways_reactome %>%
                    transmute(
                        cluster,
                        reactome,
                        names         = map_chr(names, ~paste(., collapse = '; ')),
                        seed_genes    = map_chr(seed_genes, ~paste(sort(unique(.)), collapse = '; ')),
                        seed_proteins = map_chr(seed_proteins, ~paste(sort(unique(.)), collapse = '; ')),
                        n_proteins    = map_int(proteins, ~unique(.) %>% length),
                        n_genes       = map_int(ensembl_gene_ids, ~unique(.) %>% length),
                    )
            },
            options = list(
                lengthMenu   = c(5, 10, 25, 50, 100, 10000),
                pageLength   = 10000
            )
        )
        output$tbl_candidate_pathway_clusters <- renderTable({
            tbls$candidate_pathways_reactome %>%
                filter(cluster != 'not assigned') %>%
                group_by(cluster) %>%
                summarise(
                    pathways     = n(),
                    genes        = unlist(ensembl_gene_ids) %>% unique %>% length,
                    `seed genes` = unlist(seed_genes) %>% unique %>% sort %>% paste(collapse = '; ')
                )
        })
    }
    observeEvent(input$queryPathways, {
        withProgress(
            message = 'querying reactome pathways',
            value   = 1,
            max     = tbls$seed_genes() %>%
                filter(!is.na(ensembl_gene_id)) %>%
                nrow(),
            tmp1 <- tbls$seed_genes() %>%
                filter(!is.na(ensembl_gene_id)) %>%
                select(external_gene_name, uniprotswissprot) %>%
                mutate(
                    reactome_pathway = map(
                        uniprotswissprot,
                        function(x) {
                            res <- get_reactome_pathways(x)
                            incProgress(1, detail = sprintf("UniProt: %s", x))
                            return(res)
                        }
                    )
                ) %>%
                unnest(reactome_pathway) %>%
                select(reactome, name, everything()) %>%
                arrange(reactome) %>%
                distinct()
        )
        withProgress(
            message = 'querying participating proteins',
            value   = 1,
            max     = nrow(tmp1),
            tbls$candidate_pathways_reactome <<- tmp1 %>%
                group_by(reactome) %>%
                nest(
                    names         = name,
                    seed_genes    = external_gene_name,
                    seed_proteins = uniprotswissprot
                ) %>%
                ungroup() %>%
                mutate_at(vars(-reactome),
                    ~map(., ~sort(unique(.[[1]], )))
                ) %>%
                mutate(
                    proteins = map(
                        reactome,
                        function(x) {
                            res <- get_participating_proteins(x)
                            incProgress(1, detail = x)
                            return(res[[1]])
                        }
                    )
                ) %>%
                mutate(
                    ensembl_gene_ids_seed = map(seed_proteins, ~map_proteins_to_genes(., tbls$ensembl)),
                    ensembl_gene_ids      = map(proteins, ~map_proteins_to_genes(., tbls$ensembl))
                ) %>%
                mutate(cluster = 'not assigned')
        )
        update_candidate_pathways_tables()
    })
    observeEvent(input$assign_cluster, {
        tbls$candidate_pathways_reactome$cluster[input$tbl_candidate_pathways_rows_selected] <<-
            input$cluster_name
        update_candidate_pathways_tables()
    })
    output$downloadClusterAssignment <- downloadHandler(
        'cluster-assignments.txt',
        function(file) {
            tbls$candidate_pathways_reactome %>%
                select(reactome, cluster) %>%
                write_delim(
                    path  = file,
                    delim = ' '
                )
        }
    )
    observeEvent(input$uploadClusterAssignmentFile, {
        tbls$candidate_pathways_reactome <<- tbls$candidate_pathways_reactome %>%
            select(-cluster) %>%
            left_join(
                input$uploadClusterAssignmentFile$datapath[1] %>%
                    read_delim(
                        delim   = ' ',
                        trim_ws = TRUE,
                        col_types = cols(
                            reactome = col_character(),
                            cluster  = col_character()
                        )
                    ),
                by = 'reactome'
            )
        update_candidate_pathways_tables()
    })

    refresh_pruning <- function() {
        k <- if (input$prune) input$pruning_distance else Inf
        tbls$pruned_pathways <- tbls$candidate_pathways_reactome %>%
            filter(cluster != 'not assigned') %>%
            group_by(cluster) %>%
            summarize_if(
                is.list,
                ~list(unlist(.) %>% unique)
            ) %>%
            mutate(
                pruned_proteins = map2(
                    proteins, seed_proteins,
                    ~.x[prune(.x, .y, k, tbls$interactions)]
                ),
                ensembl_gene_ids_pruned = map(pruned_proteins, ~map_proteins_to_genes(., tbls$ensembl)),
                igraph = map2(ensembl_gene_ids_pruned, ensembl_gene_ids_seed,
                              ~get_pathway_graph(.x, .y, tbls$gene_gene_interactions, tbls$ensembl))
            )
        output$tbl_pruned_pw_cluster <- DT::renderDataTable({
                tbls$pruned_pathways %>%
                    transmute(
                        cluster,
                        n_genes           = map_int(ensembl_gene_ids_pruned, length),
                        fraction_retained = map2_chr(ensembl_gene_ids, ensembl_gene_ids_pruned, ~sprintf('%5.1f%%', 100*length(.y)/length(.x))),
                        seed_genes        = map_chr(seed_genes, ~paste(sort(unique(unlist(.))), collapse = '; '))
                    )
            },
            options = list(
                lengthMenu = c(10, 25, 100),
                pageLength = 25
            )
        )
        clusters <- tbls$pruned_pathways$cluster %>% unique %>% sort
        output$plotSelector <- renderUI({
            selectInput("selectPlot", "select pathway cluster", clusters, clusters[1])
        })
        output$plotPathwayCluster <- renderPlot({
                tmp <- tbls$pruned_pathways_gene_space %>%
                    filter(cluster == clusters[1]) %>%
                    pull(igraph)
                if (length(tmp) != 1) stop('did not find unique pathway cluster graph')
                gr         <- tmp[[1]]
                n_vertices <- igraph::V(gr) %>% length()
                ggr        <- tidygraph::as_tbl_graph(gr)
                ggr %>%
                    activate(nodes) %>%
                    mutate(
                        bla = map2_chr(external_gene_name, name, ~paste(c(.x, .y), collapse = '\n\r'))
                    ) %>%
                    ggraph::ggraph() +
                    ggraph::geom_edge_link(alpha = .1, size = .33) +
                    ggraph::geom_node_point(aes(color = seed_gene), size = 3) +
                    ggraph::geom_node_text(aes(label = bla), size = 1.5, colour = 'white', vjust = 0.4)
                # ggsave('test.pdf', width = sqrt(n_vertices), height = sqrt(n_vertices))
            },
            height = 72*sqrt(2*10), width = 72*sqrt(2*10), res = 72
        )
    }
    observeEvent(input$refreshPruning, {refresh_pruning()})
    output$plotSelector <- renderUI({selectInput("selectPlot", "select pathway cluster", c())})
}

library(shiny)
library(glue)
library(tidyverse)

map(
    list.files('scripts', pattern = '.+\\.R$', full.names = TRUE),
    ~source(., verbose = FALSE)
)

dir.create('data', showWarnings = FALSE)
download_zenodo <- function(filename, doi = '3533634') {
    if (!file.exists(glue('data/{filename}')))
        download.file(
            glue('https://zenodo.org/record/{doi}/files/{filename}?download=1'),
            glue('data/{filename}')
        )
}
c('tbl_ensembl.rds', 'tbl_interactions_gene_gene.rds', 'tbl_interactions_prot_prot.rds') %>%
    map(download_zenodo)

# ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# write_rds(tbls$ensembl, 'data/tbl_ensembl.rds', compress = 'gz')
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
#
# tbls$interactions_gene_gene <- tbls$interactions_prot_prot %>%
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
#     select(-interactor_A, -interactor_B) %>%
#     distinct() %>%
#     filter(gene_A != gene_B) %>%
#     mutate(
#         tmp = map2(gene_A, gene_B, ~sort(c(.x, .y)))
#     ) %>%
#     filter(map_int(tmp, length) == 2) %>%
#     mutate(
#         gene_A = map_chr(tmp, ~.[[1]]),
#         gene_B = map_chr(tmp, ~.[[2]])
#     ) %>%
#     select(-tmp) %>%
#     distinct()
# write_rds(tbls$interactions_gene_gene, 'data/tbl_interactions_gene_gene.rds', compress = 'gz')

tbls                        <- list()
tbls$ensembl                <- read_rds('data/tbl_ensembl.rds')
tbls$interactions_prot_prot <- read_rds('data/tbl_interactions_prot_prot.rds')
tbls$interactions_gene_gene <- read_rds('data/tbl_interactions_gene_gene.rds')


renderDataTable_custom <- function(tbl) DT::renderDataTable({tbl}, options = list(lengthMenu = c(5, 15, 50, 10000), pageLength = 10000), escape = FALSE)


server <- function(input, output) {

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
    output$downloadSeedGenes <- downloadHandler(
        'seed-genes.csv',
        function(file) write_csv(tbls$seed_genes, path = file)
    )



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
    observeEvent(input$queryPathways, {
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
                            get_reactome_pathways(x)
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
    })
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
                        pathways     = n(),
                        genes        = bind_rows(genes) %>% pull(ensembl_gene_id) %>% unique %>% length,
                        `seed genes` = bind_rows(seed_data) %>% pull(external_gene_name) %>% unique %>% sort %>% paste(collapse = '; ')
                    )
            } else {
                tibble(
                    pathways     = 0,
                    genes        = 0,
                    `seed genes` = 'none'
                )
            }

        })
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
            res <- input$uploadGeneCoverageFile$datapath[1] %>%
                read_delim(
                    delim   = ' ',
                    trim_ws = TRUE,
                    col_types = cols(
                        ensembl_gene_id = col_character()
                    )
                ) %>%
                distinct()
        }
        return(res)
    })
    refresh_pruning <- function() {
        k <- if (input$prune) input$pruning_distance else Inf
        covered <- function(ensembl_gene_ids) {
            if (is.null(tbls$coverage())) {
                return(ensembl_gene_ids == ensembl_gene_ids)
            } else {
                return(ensembl_gene_ids %in% tbls$coverage()$ensembl_gene_id)
            }
        }
        # browser()
        tbls$pruned_pathways <- tbls$candidate_pathways %>%
            filter(cluster != 'not assigned') %>%
            mutate(reactome = as.list(reactome)) %>%
            group_by(cluster) %>%
            summarize_all(~list(unlist(.) %>% unique)) %>%
            mutate(
                pruned_proteins = map2(
                    proteins, seed_proteins,
                    ~.x[prune(.x, .y, k, tbls$interactions_prot_prot)]
                ),
                ensembl_gene_ids_pruned = map(pruned_proteins, ~map_proteins_to_genes(., tbls$ensembl)),
                ensembl_gene_ids_seed_covered = map(ensembl_gene_ids_seed, ~.[covered(.)]),
                ensembl_gene_ids_pruned_covered = map(ensembl_gene_ids_pruned, ~.[covered(.)])
            ) %>%
            nest(
                gene_data = starts_with('ensembl')
            ) %>%
            mutate(
                igraph_pruned = map(gene_data, ~get_pathway_graph(., tbls$interactions_gene_gene, tbls$ensembl))
            ) %>%
            select(cluster, reactome, names, gene_data, igraph_pruned)
        output$tbl_pruned_pw_cluster <- DT::renderDataTable({
                tbls$pruned_pathways %>%
                    transmute(
                        cluster,
                        n_genes            = map_int(gene_data, ~length(.$ensembl_gene_ids[[1]])),
                        n_pruned           = map_int(gene_data, ~length(.$ensembl_gene_ids_pruned[[1]])),
                        n_pruned_covered   = map_int(gene_data, ~length(.$ensembl_gene_ids_pruned_covered[[1]])),
                        seed_genes         = map_chr(gene_data, ~paste(map_ensembl_to_name(.$ensembl_gene_ids_seed[[1]]), collapse = '; ')),
                        seed_genes_covered = map_chr(gene_data, ~paste(map_ensembl_to_name(.$ensembl_gene_ids_seed_covered[[1]]), collapse = '; '))
                    )
            },
            options = list(
                lengthMenu = c(10, 25, 100),
                pageLength = 25
            )
        )
        clusters <- tbls$pruned_pathways$cluster %>% sort
        output$plotSelector <- renderUI({
            selectInput("selectPlot", "select pathway cluster", clusters, clusters[1])
        })
        output$plotPathwayCluster <- renderPlot({
                browser()
                tmp <- tbls$pruned_pathways %>%
                    filter(cluster == clusters[1]) %>%
                    pull(igraph_pruned)
                if (length(tmp) != 1) stop('did not find unique pathway cluster graph')
                gr  <- tmp[[1]]
                ggr <- tidygraph::as_tbl_graph(gr)
                ggr %>%
                    filter(covered) %>%
                    mutate(
                        bla = map2_chr(external_gene_name, name, ~paste(c(.x, .y), collapse = '\n\r'))
                    ) %>%
                    ggraph::ggraph() +
                    ggraph::geom_edge_link(alpha = .1) +
                    ggraph::geom_node_point(aes(color = seed_gene), size = 3) +
                    ggraph::geom_node_text(aes(label = bla), size = 1.5, colour = 'white', vjust = 0.4)
                ggsave('test.pdf', width = 20, height = 20)
            },
            height = 72*sqrt(2*10), width = 72*sqrt(2*10), res = 72
        )
    }
    observeEvent(input$refreshPruning, {refresh_pruning()})
    output$plotSelector <- renderUI({selectInput("selectPlot", "select pathway cluster", c())})
}

library(shiny)

ui <- navbarPage("Pathway Curation Tool",
    tabPanel("Info",
        includeMarkdown('readme.md')
    ),
    tabPanel("Upload Seed Genes",
        sidebarLayout(
            sidebarPanel(
                includeMarkdown('markdown/seed_genes_sidebar.md'),
                fileInput('uploadSeedGenesFile', 'upload seed-genes file',
                    accept = c('.xlsx')
                )
            ),
            mainPanel(
                includeMarkdown('markdown/seed_genes_header.md'),
                br(), br(),
                DT::dataTableOutput('tblSeedGenes')
            )
        )
    ),
    tabPanel("Define Pathway Cluster",
        sidebarLayout(
            sidebarPanel(
                includeMarkdown('markdown/define_pathway_cluster_sidebar_query.md'),
                actionButton('queryPathways', 'query reactome.org', width = '100%'),
                br(), hr(),
                includeMarkdown('markdown/define_pathway_cluster_sidebar.md'),
                actionButton('assign_cluster', 'assign/unassign selected', width = '100%'),
                br(), br(),
                downloadButton('downloadClusterAssignment', 'download assignment'),
                br(), br(),
                fileInput('uploadClusterAssignmentFile', 'upload assignment',
                          accept = c('.csv')
                ),
                hr(),
                h3('Summary'),
                tableOutput('tbl_candidate_pathway_clusters')
            ),
            mainPanel(
                includeMarkdown('markdown/define_pathway_cluster_header.md'),
                br(), br(),
                DT::dataTableOutput('tbl_candidate_pathways')
            )
        )
    ),
    tabPanel("Prune Pathway Cluster",
        sidebarLayout(
            sidebarPanel(
                withMathJax(includeMarkdown('markdown/prune_sidebar.md')),
                numericInput('pruning_distance', label = 'k',
                            min = 0, max = NA, value = 2
                ),
                numericInput('minScore', 'minimumScore', min = 0, max = 1, step = .01, value = .9),
                hr(),
                includeMarkdown('markdown/prune_sidebar2.md'),
                fileInput('uploadGeneCoverageFile', 'upload gene coverage file',
                          accept = c('.csv')
                ),
                actionButton('intersectGREX', 'apply', width = '100%'),
                hr(),
                downloadButton("downloadGraph", "download graph data"),
            ),
            mainPanel(
                includeMarkdown('markdown/prune.md'),
                br(), br(),
                DT::dataTableOutput('tbl_pruned_pw_cluster')
            )
        )
    ),
    tabPanel("Plot",
         sidebarPanel(
             uiOutput('plotSelector'),
             numericInput('plotWidth', label = 'width (in)',
                          min = 1, max = 50, value = 15
             ),
             numericInput('plotHeight', label = 'height (in)',
                          min = 1, max = 50, value = 10
             ),
             numericInput('labelSize', label = 'label size',
                          min = 1, max = 20, value = 5
             ),
             numericInput('pointSize', label = 'point size',
                          min = 1, max = 50, value = 5
             ),
             checkboxInput('plotNonCovered', label = 'plot non-covered', value = FALSE, width = '100%'),
             radioButtons('ensemblId', '', choiceNames = c('gene name', 'ensembl ID'), choiceValues = c('external_gene_name', 'ensembl_gene_id'), selected = 'external_gene_name', width = '100%'),
             hr(),
             downloadButton('downloadPlot', 'download current plot')
         ),
         mainPanel(
             plotOutput('plotPathwayCluster')
         )
     )
)

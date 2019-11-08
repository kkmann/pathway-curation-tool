library(shiny)

ui <- navbarPage("Pathway Curation Tool",
    tabPanel("Info"),
    tabPanel("Seed Genes",
        sidebarLayout(
            sidebarPanel(
                fileInput('uploadSeedGenesFile', 'upload seed-genes file',
                    accept = c('.txt')
                ),
                downloadButton('downloadSeedGenes', 'adsas')
            ),
            mainPanel(
                h1('Seed-genes not found in ensembl'),
                DT::dataTableOutput('tblSeedGenesNotFound'),
                h1('Seed-genes found in ensembl'),
                DT::dataTableOutput('tblSeedGenes')
            )
        )
    ),
    tabPanel("Pathway Curation",
        sidebarLayout(
            sidebarPanel(
                actionButton('queryPathways', 'query pathways', width = '100%'),
                textInput('cluster_name', label = 'Cluster name',
                          value = 'not assigned', placeholder = 'delete assignment'
                ),
                actionButton('assign_cluster', 'assign selected')
            ),
            mainPanel(
                h1('Reactome pathways'),
                DT::dataTableOutput('tbl_candidate_pathways'),
                h1('Gene names not found in reactome pathways'),
                DT::dataTableOutput('tbl_candidate_pathways_not_found'),
                h1('Aggregation by pathway cluster'),
                DT::dataTableOutput('tbl_candidate_pathway_clusters')
            )
        )
    ),
    tabPanel("Pathway Pruning",
        sidebarLayout(
            sidebarPanel(
                numericInput('pruning_distance', label = 'max indirect interactions',
                            min = 0, max = NA, value = 1
                ),
                checkboxInput('prune', label = 'prune', value = 0)
            ),
            mainPanel(
                h1('Gene names that could not be mapped to pathways'),
                DT::dataTableOutput('tbl_pruned_pw_cluster')
            )
        )
    ),
    tabPanel("Output",
         sidebarPanel(
             downloadButton("downloadData", "Download Data")
         ),
         mainPanel()
     )
)

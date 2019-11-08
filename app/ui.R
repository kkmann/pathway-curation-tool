library(shiny)

ui <- navbarPage("Pathway Curation Tool",
    tabPanel("Info"),
    tabPanel("Seed Genes",
        sidebarLayout(
            sidebarPanel(
                fileInput('uploadSeedGenesFile', 'upload seed-genes file',
                    accept = c('.txt')
                ),
                downloadButton('downloadSeedGenes', 'download seed genes')
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
                actionButton('assign_cluster', 'assign selected', width = '100%'),
                downloadButton('downloadClusterAssignment', 'download cluster assignments'),
                fileInput('uploadClusterAssignmentFile', 'upload cluster assignment file',
                          accept = c('.txt')
                ),
                h1('cluster summary'),
                tableOutput('tbl_candidate_pathway_clusters')
            ),
            mainPanel(
                DT::dataTableOutput('tbl_candidate_pathways')
            )
        )
    ),
    tabPanel("Pathway Pruning",
        sidebarLayout(
            sidebarPanel(
                actionButton('refreshPruning', 'refresh', width = '100%'),
                numericInput('pruning_distance', label = 'max indirect interactions',
                            min = 0, max = NA, value = 1
                ),
                checkboxInput('prune', label = 'prune', value = 0),
                fileInput('uploadGeneCoverageFile', 'upload gene coverage file',
                          accept = c('.txt')
                ),
            ),
            mainPanel(
                h1('Pruned Pathways'),
                DT::dataTableOutput('tbl_pruned_pw_cluster')
            )
        )
    ),
    tabPanel("Output",
         sidebarPanel(
             uiOutput('plotSelector'),
             downloadButton("downloadData", "Download Data")
         ),
         mainPanel(
             plotOutput('plotPathwayCluster')
         )
     )
)

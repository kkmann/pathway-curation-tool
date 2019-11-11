# Candidate Gene Pathways

**tl;dr:** elicit candidate gene networks based on seed-genes.

## Details

PrediXcan allows the estimation of genetically regulated gene expression (GREX)
levels by leveraging GTEx donor data.
Since the vast majority of the genes for which GREX can be estimate are 
protein coding genes, this tool is restricted to the curation of gene-networks 
containing soly *coding* genes.
We do so in steps:

1. Upload a .csv file with ensembl IDs of seed-genes. 
These can be genes for which prior associations with the outcome of interest 
have been established or which are known to play an important role in pathays
*a prior* cosidered to be of interest for outcome. The .csv file must contain a
column `ensemble_gene_id` with the ensembl IDs of the seed-genes 
(without version sugffix!).

2. Query biological pathways from reactome and hand-curate the list of 
pathways. The gene-network consists of the union of all protein coding
genes in the selcted reactome pathways.

3. Prune the gene network using protein-protein interaction data from reactome
and string db.
To reduce the size of very large gene sets, we suggest to prune by imposing 
a maximal distance of indirect associations with one of the seed genes in the
network. I.e., only genes will be retained that are at most k-edges from a
seed gene. Edges represent protein-protein interactions.

4. Intersect gene set with available GREX data. We then augment the pathways 
with information about the actually available GREX data.

5. Plot and download.

# Candidate Gene Pathways

**tl;dr:** elicit candidate gene networks based on seed-genes.

PrediXcan allows the estimation of genetically regulated gene expression (GREX)
levels by leveraging GTEx donor data.
Since the vast majority of the genes for which GREX can be estimate are 
protein coding genes, this tool is restricted to the curation of gene-networks 
containing soly *coding* genes.

## 1. Define Seed-Genes

Upload a .csv file with ensembl IDs of seed-genes. 
These can be genes for which prior associations with the outcome of interest 
have been established or which are known to play an important role in pathays
*a prior* of interest. 
The .csv file must contain a column ensemble_gene_id with the ensembl IDs of 
the seed-genes (without version sugffix!).

## 2. Query reacrome.org pathways

Query biological pathways from reactome.org and hand-curate the list of 
pathways to be included in the pathway cluster. 
The pathway cluster is then the set of (coding) genes in the union of the 
selected pathways. 
Note that for seed-genes that cannot be mapped to reactome pathways, 
we allow direct inclusion as pseudo-pathway only condaining the respective
coded gene (denoted as UniProt:*** instead of reactome:***).

## 3. Pruning

Prune the gene network using protein-protein interaction data from reactome
and string db.
To reduce the size of very large gene sets, we suggest to prune by imposing 
a maximal distance in terms of interactions with one of the seed-genes in the
cluster.
I.e., only genes that are at most k-edges from a seed-gene will be retained. 
Edges represent protein-protein interactions.
We also allow the annotation of the resulting gene-networks with coverage
data (list of ensembl IDs of genes for which GREX is actually available).

## 4. Plotting

Plot the final pathway cluster gene network graph to verify curation result.

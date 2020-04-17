# GO Enrichment Analyses

This repository contains a collection of R and python scripts that may be useful for conducting gene ontology (GO) enrichment analyses.

The script ```00scripts/GOannotation_2_map.py``` transforms a many-to-many file containing gene IDs in the first column with GO IDs in the second column into a map file that can be read by the ```00scripts/RunTopGO.R``` script which does the enrichment testing.

## Methods employed

The enrichment tests are based on the methods implemented in the [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) R package available on BioConductor.

Briefly, a Fisher's Exact Test is calculated for each GO term using the "weight01" algorithm of Alexa et al. ([2006](http://doi.org/10.1093/bioinformatics/btl140)).

Only GO terms with more than 5 annotated genes in the reference set are considered.

### Running the script

The script requires 3 input files:
* A gene ID to GO annotation map file with a separate line for each gene (e.g. geneID1<tab>GOID1, GOID2, etc.)
* A list of the genes in the reference set
* A list of genes to test for enrichment

Example:
```
00scripts/RunTopGO.R 01data/geneID_2_GO_map.txt 01data/ref_gene_set.txt 01data/test_gene_set.txt
```

Results are saved into the `02results` directory using the name of the test gene set input file

### Disclaimer

These scripts are provided with no guarantee of usefulness to anyone, including the author.


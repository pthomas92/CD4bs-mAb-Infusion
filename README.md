# Code availability for publication of Thomas *et al.,* 2024

## Data

Data provided in the {input path to data after deposition}:
- Paired VH and VK sequencing data from antigen specific B cells
- Unpaired bulk VH sequencing data from IgG enriched plasma cells
- Affinity obtained per well

### Data provided in the following Mendeley data link:
- Thomas, Peter (2024), “2024-Bulk_SC BCR Sequencing Data”, Mendeley Data, V1, doi: 10.17632/5k5thsgxw7.1

Sequence data is provided in Changeo format<sup>1</sup> `CLONE` column contains lineage definition and is consistent between antigen specific VH and bulk VH sequences:

<sup>1</sup>: Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv359

## B cell network generation

R code in the scripts/ directory describes the process by which B cell network plots were generated in the following publication:
- High affinity mAb infusion can enhance maximum affinity maturation during HIV Env immunization; doi: 10.1016/j.isci.2024.109495
- Contact peter.thomas.19@ucl.ac.uk (author) with any queries

Required packages:
- `tidyverse` (script uses v1.3.2)
- `igraph` (script uses v1.3.5)
- `RColorBrewer` (script uses v1.1.3)
- `ape` (script uses v5.6.2)

By default, the script will check if a version of the above packages has been installed, and install the latest version if none exists. Functionality should be retained across later package versions but this is not guaranteed.

### Outline
The purpose of this script is to allow network visualisation of BCR sequences, using the `igraph` R package. Briefly:
1. A user-defined data frame is loaded into the session and partitioned based on previously defined lineages definitions.
2. A levenshtein distance matrix of the nucleotide VDJ sequence is constructed per lineage.
3. Neighbour joining trees are created per lineage, which defines the edges for each node.
4. Edgelists are generated for each trees, and concatenated to create a single edgelist for the starting data.
5. Metadata is added to allow colouring of nodes by a categorical variable.
6. Network is plotted

### Script Requirements
- Must have at least 3 unique sequences within the cluster for creation of the nj tree. If this is not the case, it's recommended to deduplicate sequences and change the size of the duplicated nodes to reflect number of duplicates.

## Customising Colour
Script is set up to colour by IGHV gene family, however no variables are hard-coded. These variables can therefore be easily altered to create new colouring sets, providing they're categorical. Mapping a continuous colouring variable requires definition of a colouring gradient (for example: https://gist.github.com/cshukai/6911981)

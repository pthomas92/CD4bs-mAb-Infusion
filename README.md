# B Cell Network Generation

This repository describes the process by which B cell network plots were generated in the following publication:
- input publication details when possible
- Contact peter.thomas.19@ucl.ac.uk with any queries

Required packages:
- `tidyverse` (script uses v1.3.2)
- `igraph` (script uses v1.3.5)
- `RColorBrewer` (script uses v1.1.3)
- `ape` (script uses v5.6.2)

By default, the script will check if a version of the above packages has been installed, and install the latest version if none exists. Functionality should be retained across later package versions but this is not guaranteed.

## Outline
The purpose of this script is to allow network visualisation of BCR sequences, using the `igraph` R package. Briefly:
1. A user-defined data frame is loaded into the session and partitioned based on previously defined lineages definitions.
2. A levenshtein distance matrix of the nucleotide VDJ sequence is constructed per lineage.
3. Neighbour joining trees are created per lineage, which defines the edges for each node.
4. Edgelists are generated for each trees, and concatenated to create a single edgelist for the starting data.
5. Metadata is added to allow colouring of nodes by a categorical variable.
6. Network is plotted

## Script Requirements
- Must have at least 3 unique sequences within the cluster for creation of the nj tree. If this is not the case, it's recommended to deduplicate sequences and change the size of the duplicated nodes to reflect number of duplicates.

## Customising Colour
Script is set up to colour by IGHV gene family, however no variables are hard-coded. These variables can therefore be customised to easily create new colouring sets.

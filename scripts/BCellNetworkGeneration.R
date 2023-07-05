####################################################################################################################
#                                                                                                                  #
# Script written by Peter Thomas, Laura McCoy's Lab, UCL IIT, 2023                                                 #
#                                                                                                                  #
# Manuscript currently under review                                                                                #
#                                                                                                                  #
# Outline:                                                                                                         #
#   For each defined group (i.e., immunisation group/epitope specificity), it partitions sequences into            #
#   lineages, and then creates a neighbour joining tree to capture relationships between nodes, which              #
#   it then plots. Intermediate nodes created by the NJ tree are shrunk until not visible.                         #
#                                                                                                                  #
# Requirements:                                                                                                    #
#   A source file (can be either csv or tsv/tab, specify with 'file_sep' argument) with:                           #
#     * a unique sequence identifier column (name input to 'seq_id' variable)                                      #
#     * a lineage definitions column (name input to 'clone_' variable)                                             #
#     * a sequence for arranging positions within lineages (name input to 'sequence_' variable)                    #
#     * germline gene definitions column (name input to 'v_column')                                                #
#     * percentage alignment of v sequence to germline (e.g., column 'V_IDENTITY in Change-O output)               #
#                                                                                                                  #
# Outputs:                                                                                                         #
#   A single network graphic coloured as desired.                                                                  #
#   Supports colouring by any categorical variable, however new palette will need to be defined if ≥ 8 groups      #
#                                                                                                                  #
####################################################################################################################


##### Define functions for data processing #####
load_installPackages <- function(){
  
  ### vector for packages needed
  requiredPackages <- c('tidyverse', 'igraph', 'RColorBrewer', 'ape')
  versions <- c('1.3.2', '1.3.5', '1.1.3', '5.6.2')
  
  ### get all packages installed
  packages <- installed.packages()
  
  ### find packages not installed
  packages_to_install <- requiredPackages %in% packages == F
  
  ### loop through vector and install ones marked as T
  for(i in 1:length(packages_to_install)){
    
    if(packages_to_install[i] == T){
      
      cat(paste('>Package', requiredPackages[i], 'not installed. Installing...'))
      
      install.packages(requiredPackages[i])
      
    }
    
  }
  
  cat('>Loading packages...\n')
  
  ### loop through packages and load all sequentially
  for(i in 1:length(packages_to_install)){
    
    library(requiredPackages[i], character.only = TRUE)
    
  }
  
  message(paste('>Original script uses ',
                paste(paste(requiredPackages, paste('v', versions, sep = '')), collapse = ', '),
                '. Compatibility is not assessed for other package versions.', sep = ''))
  
  cat('>Packages loaded\n')
  
}

makeNetworks <- function(clone_list = data,
                         id_col = 'SEQUENCE_ID',
                         sequence_col = 'SEQUENCE_VDJ',
                         doublet_order = 'V_IDENTITY',
                         germ_distance = 'V_IDENTITY'){
  
  network_data <- list()
  
  ### loop through each clone list item
  for(clone in 1:length(clone_list)){
    
    ### if lineage is expanded
    if(nrow(clone_list[[clone]]) >= 3){
      
      ### create a distance matrix across the VDJ nucleotide sequence
      mat = adist(clone_list[[clone]][,sequence_col])
      ### rename columns and rows to the sequence id
      colnames(mat) = rownames(mat) = clone_list[[clone]][,id_col]
      
      g = nj(mat)
      brlen = compute.brlen(g)
      g = as.igraph(g)
      
      el = as.data.frame(get.edgelist(g))
      colnames(el) = c('source', 'dest')
      el$weight = brlen$edge.length
      
      el$source <- ifelse(grepl('Node', el$source),
                          paste(unique(clone_list[[clone]]$CLONE),
                                el$source, sep = '_'), el$source)
      
      el$dest <- ifelse(grepl('Node', el$dest),
                          paste(unique(clone_list[[clone]]$CLONE),
                                el$dest, sep = '_'), el$dest)
      
    } else if(nrow(clone_list[[clone]]) == 2) {
      
      first = which.min(clone_list[[clone]][,doublet_order])
      last = seq(1, nrow(clone_list[[clone]]))[-first]
      
      el = data.frame(source = clone_list[[clone]][,id_col][first],
                      dest = clone_list[[clone]][,id_col][last],
                      weight = (clone_list[[clone]][,doublet_order][last] - clone_list[[clone]][,doublet_order][first]))
      
    } else {
      
      el = data.frame(source = clone_list[[clone]][,id_col],
                      dest = clone_list[[clone]][,id_col],
                      weight = NA)
      
    }
    
    ### package into output list
    network_data[[clone]] <- el
    
  }
  
  ### join list into a dataframe + reference list (was in a hurry writing it so duplicated this)
  network_data <- do.call('rbind', network_data)
  data = do.call('rbind', clone_list)
  
  ### create igraph object and remove self-linkages
  g <- igraph::graph_from_data_frame(network_data)
  g <- igraph::simplify(g)
  
  V(g)$germ_dist <- NA
  V(g)$sample <- NA
  
  ### create a new vertex attribute for the v identity of each vertex
  for(i in 1:length(V(g))){
      
      if(grepl('Node', V(g)$name[i]) == F){
        
        V(g)$germ_dist[i] <- data[grep(paste(V(g)$name[i], '$', sep = ''), data[, id_col]), germ_distance]
        V(g)$size[i] <- 6
        
      } else {
        
        V(g)$size[i] <- 0.001
        
      }
    
  }
  
  return(g)
  
}

addMetaData <- function(network_graph,
                        metadata,
                        column_add,
                        name_add,
                        seq_id = seq_id){
  
  naming_vector = rep(NA, length(V(network_graph)$name))
  
  for(i in 1:length(naming_vector)){
    
    if(grepl('Node', V(network_graph)$name[i])){
      
      next
      
    } else {
      
      naming_vector[i] <- metadata[V(network_graph)$name[i] == metadata[,seq_id], column_add]
      
    }
    
  }
  
  network_graph <- set.vertex.attribute(network_graph, name = name_add, value = naming_vector)
  
  return(network_graph)
  
}

setColours <- function(vertex_characteristic,
                       g,
                       palette_ = cbbPalette){
  
  ### assign the names of the colour palette to be unique elements in the desired colouring column
  names(palette_) <- 
    sort(unique(
      get.vertex.attribute(g, vertex_characteristic)[!is.na(get.vertex.attribute(g,vertex_characteristic))]
      ))
  
  ### if vertex attribute is not NA (i.e., it's not an intermediate node in the nj tree), assign the colour named above
  colours_ = sapply(get.vertex.attribute(g, vertex_characteristic), function(x, cp = palette_){
    return(ifelse(!is.na(x), unlist(cp[names(cp) == x]), '#FFFFFF'))
  })
  
  ### assign new colours to colours attribute
  V(g)$colour = colours_
  
  ### set legend information (colours and text), removing NAs
  col_ = unique(colours_[order(names(colours_))])
  legend_ = unique(names(colours_)[order(names(colours_))])
  legend_ = legend_[!is.na(legend_)]
  col_ = col_[1:length(legend_)]
  
  return(list(graph = g,
              legend_ = list(legend_colours = col_,
                             legend_label = legend_)))
  
}

#####

##### load packages and set constants for the script #####

load_installPackages()
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

##### Variables used within the paper #####
path_to_file = '/path/to/file.csv'
file_sep = ','
seq_id = 'SEQUENCE_ID'
sequence_ = 'SEQUENCE_VDJ'
clone_ = 'CLONE'
v_column = 'GERMLINE_V_CALL'

### read data
sequence_data <- read.csv(path_to_file, sep = file_sep) %>%
  mutate(V_FAM = gsub('^.*?(IG.V\\d).*?$', '\\1', get(v_column), ignore.case = T))

#####

##### split sequences into lineages and create networks #####
g <- makeNetworks(split(sequence_data, sequence_data[,clone_]),
                  id_col = seq_id,
                  sequence_col = 'JUNCTION',
                  doublet_order = 'V_IDENTITY',
                  germ_distance = 'V_IDENTITY')

#####

##### generate graph layouts with graph optimisation (adds time, but improves visualisation) #####
set.seed(1234)
layout_ = layout_with_graphopt(graph = g, niter = 50000) # should be set to ≥ 20000 to resolve lineages properly

#####

##### add metadata to network #####

### metadata addition is currently through a loop (so a bit slow)
### but code was failing to save output correctly in the function when using 'set.vertex.attribute'
### will change later

### default colouring requires v family to be added
g = addMetaData(network_graph = g, metadata = sequence_data, column_add = 'V_FAM', name_add = 'v_fam')

### example of adding additional data
g = addMetaData(network_graph = g, metadata = sequence_data, column_add = 'SAMPLE_ID', name_add = 'immunisation')

#####

##### link the colour palette to nodes #####
plot_vfam = setColours(vertex_characteristic = 'v_fam',
                       g = g,
                       palette_ = cbbPalette)

vfam_graph = plot_vfam$graph
vfam_legend = plot_vfam$legend_

plot_imms = setColours(vertex_characteristic = 'immunisation',
                       g = g,
                       palette_ = cbbPalette)

imms_graph = plot_imms$graph
imms_legend = plot_imms$legend_

#####

##### plot networks #####
par(xpd = T)

plot.igraph(vfam_graph, layout = layout_,
            vertex.size = V(vfam_graph)$size, edge.arrow.size = 0.01,
            vertex.label = NA, vertex.color = V(vfam_graph)$colour)

legend('bottom', inset = c(0, -0.12),
       legend = vfam_legend$legend_label, col = vfam_legend$legend_colours,
       pch = 20, horiz = T,
       cex = 0.75, pt.cex = 1)

plot.igraph(imms_graph, layout = layout_,
            vertex.size = V(imms_graph)$size, edge.arrow.size = 0.01,
            vertex.label = NA, vertex.color = V(imms_graph)$colour)

legend('bottom', inset = c(0, -0.12),
       legend = imms_legend$legend_label, col = imms_legend$legend_colours,
       pch = 20, horiz = T,
       cex = 0.75, pt.cex = 1)

#####

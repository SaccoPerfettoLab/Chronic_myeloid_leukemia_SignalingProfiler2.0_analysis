
library(readxl)
library(tidyverse)
library(SignalingProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(writexl)
library(readxl)
library(tidyverse)
library(tidyr)
library (igraph)
library(data.table)




get_df_uni2info_2 = function(id_input, batch_size = 400) {

  # i = 18801
  # id_input = uniprot_ids$UNIPROT

  id_input = id_input[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id_input)]

  # batch_size = 10
  header_df_uni2seq_fin <- c("Entry", "Reviewed", "Entry Name", "Protein names", "Gene Names (primary)", "Organism", "Length", "Sequence")
  df_uni2seq_fin <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df_uni2seq_fin) <- header_df_uni2seq_fin

  for (i in seq(from= 1, to= length(id_input)-(batch_size-1), by = batch_size)){
    print (i)

    id= unique(id_input[i:(i+(batch_size-1))])

    query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
    url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence&format=tsv&query=',query_test)
    result <- httr::GET(url_test)
    as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
    df_uni2seq <- read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE, show_col_types =  F)
    df_uni2seq_fin <- rbind(df_uni2seq_fin, df_uni2seq )
  }

  id= unique(id_input[i:length(id_input)])
  # id = id[grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', id)]
  # id = id[!grepl('CHEBI:|SIGNOR-|NP_|CID:|SID:|PRO_|DB0|_9606', id)]
  query_test= paste0('accession%3A', paste0( id, collapse ='%20OR%20accession%3A'))
  url_test=paste0('https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_primary%2Corganism_name%2Clength%2Csequence&format=tsv&query=',query_test)
  result <- httr::GET(url_test)
  as.character(httr::content(result, "text"))-> file_uni# automatically parses JSON
  df_uni2seq <- read_delim(file_uni, delim = '\t',skip_empty_rows = TRUE,show_col_types = F)
  df_uni2seq_fin <- rbind(df_uni2seq_fin,df_uni2seq )
  df_uni2seq_fin%>%
    dplyr::select(Entry, Sequence, `Gene Names (primary)`)-> df_uni2seq_fin_sel

  df_uni2seq_fin <- df_uni2seq_fin %>% distinct()
  return(df_uni2seq_fin_sel)
}


# Functional circuits
# starts from phenotypes and retrieves upper nodes in k steps
# to a set of desired starting nodes
functional_circuits_new <-
  function (SP_object, start_nodes, phenotypes, k) {

    # SP_object = final_sp_visualization
    # start_nodes = c('BCR_ABL')
    # phenotypes = target2_df$gene_name
    # k = k

    # SP_object = final_sp_visualization
    # start_nodes = 'BCR_ABL'
    # phenotype = target1_df$gene_name
    # k = 20

    # SP_object = final_sp_visualization
    # start_nodes = 'BCR_ABL'
    # phenotypes = c('CELL CYCLE BLOCK', 'G1/S TRANSITION')
    # k = k
    #
    SP_graph <- SP_object$igraph_network

    #phenotypes <- target1_df$gene_name

    for (i in c(1:length(phenotypes))) {
      phenotype <- phenotypes[i]
      to_nodes <- V(SP_graph)$name[V(SP_graph)$name %in% start_nodes]

      all_paths <- igraph::all_simple_paths(
        SP_graph,
        from = str_replace_all(string = phenotype, pattern = "[[:space:]\\\\/]", '_'),
        to = to_nodes,
        mode = "in",
        cutoff = k
      )

      all_paths_nodes <- unique(names(unlist(all_paths)))

      if (length(all_paths_nodes) == 0) {
        warning(paste0("No path of length ", k, " have been found for ",
                       phenotype))
      }

      if (i == 1) {
        final_nodes <- all_paths_nodes
      } else {
        final_nodes <- c(final_nodes, all_paths_nodes)
      }
    }

    pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)
    pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                                 what = c("edges"))


    pheno_circuit_no_incoming <-
      igraph::graph_from_data_frame(d = pheno_circuit_edges,
                                    vertices = igraph::as_data_frame(pheno_circuit, what = c("vertices")))

    return(pheno_circuit_no_incoming)
  }

# This function removes nodes in the graph with a carnival_activiti != |100|
# and its interactors if the removed node are the only interactors
remove_nodes_and_interactors <- function(graph) {

  circuits_nod <- as_data_frame(graph, what = 'vertices')

  circuits_nod %>% filter(!carnival_activity %in% c(100, -100)) -> to_remove

  nodes_to_remove = V(graph) [V(graph)$name %in% to_remove$name]

  neighbors_to_remove <- c()
  for (node in nodes_to_remove) {
    neighbors <- igraph::neighbors(graph = graph, v = V(graph) [ node])

    for (neighbor in neighbors) {
      if(length( neighbors(graph, V(graph)[neighbor])) != 0){
        if (all( neighbors(graph, V(graph)[neighbor])) %in% nodes_to_remove) {
          neighbors_to_remove <- c(neighbors_to_remove, neighbor)
        }
      }
    }
    neighbors_to_remove <- c(neighbors_to_remove, node)
  }

  graph_clean <- delete.vertices(graph, neighbors_to_remove)
  return(graph_clean)
}

# This function takes as input a graph
# and a set of paths and returns a dataframe with the edges
filter_paths_raw <- function(graph, all_paths, causality, mode = 'out'){

  keep_vector <- c()
  #i_path = 4
  for(i_path in 1:length(all_paths)){

    path <- all_paths[[i_path]]

    # Transform a path in a vector of nodes id
    path_nodes_ids <- as.vector(path)

    # Transform the nodes ids in pairs of nodes id identifying edges
    path_pairs <- rep(path_nodes_ids, times = c(1, rep(2, length(path_nodes_ids)-2), 1))

    if(mode == 'in'){path_pairs <- rev(path_pairs)}

    # Extract the edges ids
    edges_in_path <- get.edge.ids(graph = graph,vp = path_pairs)

    # Extract the attributes of the edges
    edge_attributes <- as.numeric(as.vector(E(graph)$sign[edges_in_path]))

    # If the product is inhibitory keep the path

    if(!is.null(causality)){
      # This chunk filter paths for global causal effect
      if(causality == '-1'){
        keep <- ifelse(prod(edge_attributes) < 0, TRUE, FALSE)
      }else if(causality == '1'){
        keep <- ifelse(prod(edge_attributes) > 0, TRUE, FALSE)
      }

      keep_vector <- c(keep_vector, keep)
    }
  }
  filtered_paths <- all_paths[keep_vector]

  return(filtered_paths)
}


# Functional circuits
# starts from phenotypes and retrieves upper nodes in k steps
# to a set of desired starting nodes
functional_circuits_new_coherent <-
  function (SP_object, start_nodes, phenotypes, k) {

    SP_object = final_sp_visualization
    start_nodes = 'BCR_ABL'
    phenotype = c('DNA REPAIR')
    k = 8


    SP_graph <- SP_object$igraph_network

    for (i in c(1:length(phenotypes))) {
      phenotype <- phenotypes[i]
      to_nodes <- V(SP_graph)$name[V(SP_graph)$name %in% start_nodes]

      all_paths <- igraph::all_simple_paths(
        SP_graph,
        from = phenotype,
        to = to_nodes,
        mode = "in",
        cutoff = k
      )


      all_paths_nodes <- unique(names(unlist(all_paths_filtered)))

      if (length(all_paths_nodes) == 0) {
        warning(paste0("No path of length ", k, " have been found for ",
                       phenotype))
      }

      if (i == 1) {
        final_nodes <- all_paths_nodes
      } else {
        final_nodes <- c(final_nodes, all_paths_nodes)
      }
    }

    pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)
    pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                                 what = c("edges"))


    pheno_circuit_no_incoming <-
      igraph::graph_from_data_frame(d = pheno_circuit_edges,
                                    vertices = igraph::as_data_frame(pheno_circuit, what = c("vertices")))


    RCy3::createNetworkFromIgraph(pheno_circuit_no_incoming,
                                  title = 'DNA REPAIR')

    # data_path <- system.file("extdata", "SP_pheno_layout.xml", package = "SignalingProfiler")
    # RCy3::importVisualStyles(filename = data_path)

    RCy3::setVisualStyle('SP_pheno_layout')
    return(pheno_circuit_no_incoming)
  }

pheno_to_start_circuit <- function(SP_object, start_nodes, phenotypes, k, start_to_top = FALSE) {

  # SP_object = opt1$sp_object_phenotypes
  # start_nodes = start_nodes
  # phenotype = phenotypes[i_pheno]
  # k = k_vector[i_pheno]

  SP_graph <- SP_object$igraph_network

  for (i in c(1:length(phenotypes))) {
    phenotype <- phenotypes[i]
    to_nodes <- V(SP_graph)$name[V(SP_graph)$name %in% start_nodes]

    all_paths <- igraph::all_simple_paths(
      SP_graph,
      from = stringr::str_replace_all(string = phenotype, pattern = "[[:space:]\\\\/]", '_'),
      to = to_nodes,
      mode = "in",
      cutoff = k
    )

    all_paths_nodes <- unique(names(unlist(all_paths)))

    if (length(all_paths_nodes) == 0) {
      warning(paste0("No path of length ", k, " have been found for ",
                     phenotype))
    }

    if (i == 1) {
      final_nodes <- all_paths_nodes
    } else {
      final_nodes <- c(final_nodes, all_paths_nodes)
    }
  }

  pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)
  pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                               what = c("edges"))
  if(start_to_top == TRUE){
    pheno_circuit_edges <- pheno_circuit_edges %>% dplyr::filter(!to %in% start_nodes)
  }

  pheno_circuit_no_incoming <-
    igraph::graph_from_data_frame(d = pheno_circuit_edges,
                                  vertices = igraph::as_data_frame(pheno_circuit,
                                                                   what = c("vertices")))

  return(pheno_circuit_no_incoming)
}

pheno_to_start_circuit <- function(SP_object, start_nodes, phenotypes, k, start_to_top = FALSE) {

  # SP_object = opt1$sp_object_phenotypes
  # start_nodes = start_nodes
  # phenotype = phenotypes[i_pheno]
  # k = k_vector[i_pheno]

  SP_graph <- SP_object$igraph_network

  for (i in c(1:length(phenotypes))) {
    phenotype <- phenotypes[i]
    to_nodes <- V(SP_graph)$name[V(SP_graph)$name %in% start_nodes]

    all_paths <- igraph::all_simple_paths(
      SP_graph,
      from = stringr::str_replace_all(string = phenotype, pattern = "[[:space:]\\\\/]", '_'),
      to = to_nodes,
      mode = "in",
      cutoff = k
    )

    all_paths_nodes <- unique(names(unlist(all_paths)))

    if (length(all_paths_nodes) == 0) {
      warning(paste0("No path of length ", k, " have been found for ",
                     phenotype))
    }

    if (i == 1) {
      final_nodes <- all_paths_nodes
    } else {
      final_nodes <- c(final_nodes, all_paths_nodes)
    }
  }

  pheno_circuit <- igraph::induced_subgraph(SP_graph, final_nodes)
  pheno_circuit_edges <- igraph::as_data_frame(pheno_circuit,
                                               what = c("edges"))
  if(start_to_top == TRUE){
    pheno_circuit_edges <- pheno_circuit_edges %>% dplyr::filter(!to %in% start_nodes)
  }

  pheno_circuit_no_incoming <-
    igraph::graph_from_data_frame(d = pheno_circuit_edges,
                                  vertices = igraph::as_data_frame(pheno_circuit,
                                                                   what = c("vertices")))

  return(pheno_circuit_no_incoming)
}




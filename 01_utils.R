# Useful functions from the private repository related to the article
# Clustering and Structural Robustness in Causal Diagrams (Tikka et al., 2023)
# https://github.com/santikka/causaleffect 

library(igraph)
library(R.utils)

children <- function(x, g, v = igraph::V(g)) {
  ch_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "out"))
  v[ch_ind]$name
}

parents <- function(x, g, v = igraph::V(g)) {
  pa_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "in"))
  v[pa_ind]$name
}

descendants <- function(x, g, v = igraph::V(g)) {
  de_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "out"))
  v[de_ind]$name
}

ancestors <- function(x, g, v = igraph::V(g)) {
  an_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "in"))
  v[an_ind]$name
}

neighbors_ <- function(x, g, v = igraph::V(g)) {
  ne_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "all"))
  v[ne_ind]$name
}

connected <- function(x, g, v = igraph::V(g)) {
  co_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "all"))
  v[co_ind]$name
}

uu <- function(x) {
  if (length(x)) unique(unlist(x))
  else character(0)
}

edge_subgraph <- function(g, incoming, outgoing) {
  # Setting from and to to NULL to satisfy CRAN if we end up making a package
  # R thinks these are global bindings, but they are igraph-operators for edges
  .to <- .from <- NULL
  e <- igraph::E(g)
  e_inc <- e[.to(incoming)]
  e_out <- e[.from(outgoing)]
  igraph::subgraph.edges(g, e[setdiff(e, union(e_inc, e_out))], delete.vertices = FALSE)
}

# Convert an igraph graph using causaleffect syntax into a dag
# with explicit latent variables
to_dag <- function(g) {
  out <- g
  unobs_edges <- which(igraph::edge.attributes(g)$description == "U")
  if (length(unobs_edges)) {
    e <- igraph::get.edges(g, unobs_edges)
    e <- e[e[ ,1] > e[ ,2], , drop = FALSE]
    e_len <- nrow(e)
    new_nodes <- paste0("U[", 1:e_len, "]")
    g <- igraph::set.vertex.attribute(g, name = "description", value = "")
    g <- g + igraph::vertices(new_nodes, description = rep("U", e_len))
    v <- igraph::get.vertex.attribute(g, "name")
    g <- g + igraph::edges(c(rbind(new_nodes, v[e[ ,1]]), rbind(new_nodes, v[e[ ,2]])))
    obs_edges <- setdiff(igraph::E(g), igraph::E(g)[unobs_edges])
    out <- igraph::subgraph.edges(g, igraph::E(g)[obs_edges], delete.vertices = FALSE)
  }
  out
}

ancestors_unsrt <- function(node, G) {
  an.ind <- unique(unlist(igraph::neighborhood(G, order = igraph::vcount(G), nodes = node, mode = "in")))
  an <- igraph::V(G)[an.ind]$name
  return(an)
}

parents_unsrt <- function(node, G.obs) {
  pa.ind <- unique(unlist(igraph::neighborhood(G.obs, order = 1, nodes = node, mode = "in")))
  pa <- igraph::V(G.obs)[pa.ind]$name
  return(pa)
}

children_unsrt <- function(node, G) {
  ch.ind <- unique(unlist(igraph::neighborhood(G, order = 1, nodes = node, mode = "out")))
  ch <- igraph::V(G)[ch.ind]$name
  return(ch)
}

# Function to generate the subgraph G[A; A->B, B->A]

modify_clustered <- function(g, A, pa_X, ch_Y) {
  all_vars <- unique(c(A, pa_X, ch_Y))
  
  # Create subgraph with only these vertices
  sub_g <- induced_subgraph(g, all_vars)
  
  # Remove unwanted edges:
  # - For pa_X: remove edges that don't go TO A
  # - For ch_Y: remove edges that don't come FROM A
  
  edges_to_remove <- E(sub_g)[
    (tail_of(sub_g, E(sub_g))$name %in% pa_X & !head_of(sub_g, E(sub_g))$name %in% A) |
      (head_of(sub_g, E(sub_g))$name %in% ch_Y & !tail_of(sub_g, E(sub_g))$name %in% A)
  ]
  
  sub_g <- delete_edges(sub_g, edges_to_remove)
  return(sub_g)
}

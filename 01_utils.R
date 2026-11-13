# Some useful helper functions. 

library(igraph)
library(R.utils)

children <- function(x, g, v = igraph::V(g)) {
  ch_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "out"))
  v[ch_ind]$name
}

children_lst <- function(x, g, v = igraph::V(g), self = TRUE) {
  ch <- igraph::neighborhood(g, order = 1, nodes = x, mindist = 1 - self, mode = "out")
  setNames(lapply(ch, function(i) v[i]$name), v[x]$name)
}

parents <- function(x, g, v = igraph::V(g)) {
  pa_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "in"))
  v[pa_ind]$name
}

parents_lst <- function(x, g, v = igraph::V(g), self = TRUE) {
  pa <- igraph::neighborhood(g, order = 1, nodes = x, mindist = 1 - self, mode = "in")
  setNames(lapply(pa, function(i) v[i]$name), v[x]$name)
}

descendants <- function(x, g, v = igraph::V(g)) {
  de_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "out"))
  v[de_ind]$name
}

descendants_lst <- function(x, g, v = igraph::V(g), names = TRUE, unobs = character(0)) {
  de <- igraph::neighborhood(g, order = length(v), nodes = x, mode = "out")
  if (!names) {
    return(de)
  }
  de <- setNames(lapply(de, function(i) v[i]$name), v[x]$name)
  lapply(de, setdiff, unobs)
}

ancestors <- function(x, g, v = igraph::V(g)) {
  an_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "in"))
  v[an_ind]$name
}

ancestors_lst <- function(x, g, v = igraph::V(g), names = TRUE, unobs = character(0)) {
  an <- igraph::neighborhood(g, order = length(v), nodes = x, mode = "in")
  if (!names) {
    return(an)
  }
  an <- setNames(lapply(an, function(i) v[i]$name), v[x]$name)
  lapply(an, setdiff, unobs)
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
  if (length(x)) unique(unlist(x, use.names = FALSE))
  else character(0)
}

edge_subgraph <- function(g, e, incoming, outgoing) {
  .to <- .from <- NULL
  e_inc <- e[.to(incoming)]
  e_out <- e[.from(outgoing)]
  igraph::subgraph_from_edges(g, e[setdiff(e, union(e_inc, e_out))], delete.vertices = FALSE)
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

filter_unobs_parents <- function(pa_list, unobs) {
  from <- unlist(pa_list, use.names = FALSE)
  if (!length(from)) return(pa_list)
  grp <- rep.int(seq_along(pa_list), lengths(pa_list))
  is_u <- from %in% unobs
  uf <- from[is_u]
  keep <- !is_u | !(from %in% uf[duplicated(uf)])
  out <- split(from[keep], factor(grp[keep], levels = seq_along(pa_list)))
  names(out) <- names(pa_list)
  out
}

# Function to generate the subgraph G[A->B, B->C]

modify_clustered <- function(g, A, pa_X, ch_Y) {
  all_vars <- unique(c(A, pa_X, ch_Y))
  
  # Create subgraph with only these vertices
  sub_g <- induced_subgraph(g, all_vars)
  
  # Remove unwanted edges:
  # - For pa_X: remove edges that don't go TO A
  # - For ch_Y: remove edges that don't come FROM A
  e <- E(sub_g)
  tails <- tail_of(sub_g, e)$name
  heads <- head_of(sub_g, e)$name
  
  edges_to_remove <- e[
    !((tails %in% A & heads %in% A) |
        (tails %in% pa_X & heads %in% A) |
        (tails %in% A & heads %in% ch_Y))
  ]
  
  sub_g <- delete_edges(sub_g, edges_to_remove)
  return(sub_g)
}

modify_clustered_alt <- function(g, edges, heads, tails, A, pa_X, ch_Y) {
  all_vars <- unique(c(A, pa_X, ch_Y))
  edges_to_keep <- edges[
    (tails %in% pa_X & heads %in% A) |
      (tails %in% A & heads %in% ch_Y) |
      (heads %in% A & tails %in% A)
  ]
  sub_g <- subgraph_from_edges(g, edges_to_keep, delete.vertices = TRUE)
  return(sub_g)
}

# Check paths between parents and children of A that lie entirely in A.
check_paths <- function(A, pa_A, ch_A, rank) {
  ord <- order(rank[A])
  A <- A[ord]
  pa_A <- pa_A[ord]
  ch_A <- ch_A[ord]
  n_A <- length(A)
  
  # Edges into A: from -> A[to]; tails outside A are the parents P
  from <- unlist(pa_A, use.names = FALSE)
  to <- rep.int(seq_len(n_A), lengths(pa_A))
  from_A <- match(from, A)
  ent <- is.na(from_A)
  P <- unique(from[ent])
  n_P <- length(P)
  
  # Edges out of A: A[cfrom] -> cto; heads outside A are the children C
  cto <- unlist(ch_A, use.names = FALSE)
  cfrom <- rep.int(seq_len(n_A), lengths(ch_A))
  out <- is.na(match(cto, A))
  C <- unique(cto[out])
  n_C <- length(C)
  
  if (n_P == 0 || n_C == 0) return(TRUE)
  
  # Reverse if there are fewer children than parents
  if (n_C < n_P){
    return(Recall(A, pa_A = ch_A, ch_A = pa_A, rank = -rank))
  }
    
  # S[, i]: parents that reach A[i] by a path whose interior lies in A
  S <- matrix(FALSE, n_P, n_A)
  S[cbind(match(from[ent], P), to[ent])] <- TRUE
  par_in <- split(from_A[!ent], factor(to[!ent], levels = seq_len(n_A)))
  for (i in seq_len(n_A)) {
    pj <- par_in[[i]]
    if (length(pj) == 1L) {
      S[, i] <- S[, i] | S[, pj]
    } else if (length(pj) > 1L) {
      S[, i] <- S[, i] | (rowSums(S[, pj, drop = FALSE]) > 0)
    }
  }
  
  # Each child must be reached by every parent other than itself
  self <- match(C, P)
  exits <- split(cfrom[out], factor(match(cto[out], C), levels = seq_len(n_C)))
  for (j in seq_len(n_C)) {
    ex <- exits[[j]]
    r <- if (length(ex) == 1L) S[, ex] else rowSums(S[, ex, drop = FALSE]) > 0
    if (!is.na(self[j])) r[self[j]] <- TRUE
    if (!all(r)) return(FALSE)
  }
  return(TRUE)
}

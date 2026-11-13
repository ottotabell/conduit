# Implementation of FindConduits algorithm

source("01_utils.R")

# Expands the set of possible IPS candidates
# (line 4 of the algorithm)
# Initially C_set one of the vertices in graph g.
expandchildren <- function(C_set, pa, ch) {
  l_C_set <- length(C_set)
  # The maximum size of the IPS set is length of C_set squared
  new_C_set <- vector("list", l_C_set^2)
  new_C_set[1:l_C_set] <- C_set
  
  # To keep a tab on the size of the candidate list
  current_ind <- l_C_set
  
  for (i in 1:l_C_set) {
    upd <- TRUE
    new_c <- new_C_set[[i]]
    # Update the list while Ch(Pa(X)) \neq X,
    # then move on to the next unhandled item on the list
    while(upd) {
      old_c <- new_c
      new_par <- setdiff(uu(pa[old_c]), c(old_c))
      new_c <- setdiff(uu(ch[new_par]), new_par)
      if (setequal(new_c, old_c)) {
        upd <- FALSE
      } else {
        current_ind <- current_ind + 1
        new_C_set[[current_ind]] <- new_c
      }
    }
  }
  return(new_C_set[1:current_ind])
}

# Expands the set of possible IFS candidates
# (line 5 of the algorithm)
# The logic is analogous to expandchildren 
expandparents <- function(P_set, pa, ch) {
  l_P_set <- length(P_set)
  new_P_set <- vector("list", l_P_set^2)
  new_P_set[1:l_P_set] <- P_set
  
  current_ind <- l_P_set
  
  for (i in 1:l_P_set) {
    upd <- TRUE
    new_p <- P_set[[i]]
    while(upd) {
      old_p <- new_p
      new_ch <- setdiff(uu(ch[old_p]), old_p)
      new_p <- setdiff(uu(pa[new_ch]), c(new_ch))
      if (setequal(new_p, old_p)) {
        upd <- FALSE
      } else {
        current_ind <- current_ind + 1
        new_P_set[[current_ind]] <- new_p
      }
    }
  }
  return(new_P_set[1:current_ind])
}

# Algorithm for finding every (restricted) conduit in a given causal graph.
# MANDATORY INPUTS
# g: the DAG in the igraph format
# OPTIONAL INPUTS
# prohibit: the set of vertices prohibited to be clustered, i.e., V \cup U \setminus R
# unobs: the set of unobserved variables that are in g
# size: minimum size of the conduits returned
# empty_rec: if TRUE, the output contains also conduits with empty receiver sets
# empty_emi: if TRUE, the output contains also conduits with empty emitter sets
find_conduits <- function(g, prohibit = character(0), unobs = character(0), size = 2,
                                    empty_rec = TRUE, empty_emi = TRUE) {
  
  # Create the restriction set
  nodes <- igraph::V(g)
  n <- length(nodes)
  if (length(prohibit)) {
    restrict <- setdiff(igraph::V(g)$name, c(prohibit))
  } else {
    restrict <- igraph::V(g)$name
  }

  nodes_L <- igraph::V(g)

  observed <- setdiff(nodes_L$name, unobs)
  
  # Precomputing parents, children, ancestors and descendants
  pa_plus <- parents_lst(nodes_L, g, nodes_L, self = FALSE)
  pa <- lapply(pa_plus, \(x) setdiff(x, unobs))
  ch <- children_lst(nodes_L, g, nodes_L, self = FALSE)
  an <- ancestors_lst(nodes_L, g, nodes_L, unobs = unobs)
  de <- descendants_lst(nodes_L, g, nodes_L, unobs = unobs)
  
  # Topological order of the vertices using rank for fast reversal
  topo <- igraph::topo_sort(g, mode = "out")$name
  rk <- setNames(seq_along(topo), topo)
  
  # Initializing the iterative parental sibling (IPS) candidates
  C_set <- as.list(names(nodes_L))
  names(C_set) <- NULL
  
  # Initializing the iterative filial sibling (IFS) candidates
  P_set <- as.list(names(nodes_L))
  names(P_set) <- NULL
  
  # Expanding the IPS iteratively to create the full IPS set for all other k
  C_set <- expandchildren(C_set, pa_plus, ch)
  
  # Expanding the IFS iteratively to create the full IFS set for all other l
  P_set <- expandparents(P_set, pa_plus, ch)
  
  # Explicitly add the empty set
  C_set <- c(C_set, list(character(0)))
  P_set <- c(P_set, list(character(0)))
  
  # Remove potential duplicates
  C_set <- C_set[!duplicated(C_set)]
  P_set <- P_set[!duplicated(P_set)]
  
  conduits <- vector("list", length(C_set) * length(P_set))
  current_ind <- 0L
  A <- NULL
  
  # Check potential IPS and IFS pairs
  for (i in seq_along(C_set)) {
    X_orig <- C_set[[i]]
    if (!empty_rec & !length(X_orig)) next
    for (j in seq_along(P_set)) {
      Y_orig <- P_set[[j]]
      if (!empty_emi & !length(Y_orig)) next
      if (!length(Y_orig) & !length(X_orig)) next
      # Construct A according to the formation conditions (C1)-(C3)
      if (!length(X_orig)) A <- uu(an[Y_orig])
      if (!length(Y_orig)) A <- uu(de[X_orig]) 
      if (length(X_orig) && length(Y_orig)) {
        A <- intersect(
          unlist(an[Y_orig], use.names = FALSE), 
          unlist(de[X_orig], use.names = FALSE)
        )
      }
      # Check if A obeys the restriction rule and if it is a conduit,
      # add it to the list
      if (length(A) >= size && all(A %in% restrict)) {
        # Keep only those unobserved variables with 
        # a single child in A in the parent list
        pa_A <- filter_unobs_parents(pa_plus[A], unobs)
        ch_A <- ch[A]
        if (check_paths(A, pa_A, ch_A, rk)) {
          current_ind <- current_ind + 1L
          conduits[[current_ind]] <- A
        }
      }
    }
  }
  if (current_ind == 0) return(list())
  conduits <- conduits[1:current_ind]
  # Deduplicate ignoring element order: build a canonical key per conduit
  keys <- vapply(conduits, function(x) paste(sort(x), collapse = "\r"), character(1))
  conduits[!duplicated(keys)]
}

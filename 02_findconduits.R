# Implementation of FindConduits algorithm

source("01_utils.R")

# Algorithm 1
find_conduits <- function(g, prohibit = character(0), singletons = FALSE,
                                    empty_rec = TRUE, empty_emi = TRUE) {
  
  # Create the restriction set
  nodes <- igraph::V(g)
  n <- length(nodes)
  if (length(prohibit)) {
    restrict <- setdiff(igraph::V(g)$name, prohibit)
  } else {
    restrict <- igraph::V(g)$name
  }
  
  # Creating the subgraph L
  L <- modify_clustered(g, restrict, setdiff(parents(restrict, g), restrict), setdiff(children(restrict, g), restrict))
  
  nodes_L <- igraph::V(L)
  
  pa <- setNames(lapply(nodes_L, function(x) parents(x, L, nodes_L)[-1]), nodes_L$name)
  ch <- setNames(lapply(nodes_L, function(x) children(x, L, nodes_L)[-1]), nodes_L$name)
  an <- setNames(lapply(nodes_L, function(x) ancestors(x, L, nodes_L)), nodes_L$name)
  de <- setNames(lapply(nodes_L, function(x) descendants(x, L, nodes_L)), nodes_L$name)
  
  # Initializing the iterative parental sibling (IPS) candidates (i.e. IPS_\cL^1)
  C_set <- ch[sapply(ch, length) > 0]
  names(C_set) <- NULL
  
  # Initializing the iterative filial sibling (IFS) candidates
  P_set <- pa[sapply(pa, length) > 0]
  names(P_set) <- NULL
  
  # Add singleton components separately if requested
  if (singletons) {
    tc <- lapply(nodes$name, function(x) {
      y <- list(nodes = x, receivers = character(0), emitters = character(0))
      if (length(pa[[x]])) y$receivers <- x
      if (length(ch[[x]])) y$emitters <- x
      y
    })
  } else {
    tc <- list()
  }
  
  # Add the full graph if no restrictions
  if (length(restrict) == n) {
    tc <- c(tc, list(list(nodes = nodes$name, receivers = character(0), emitters = character(0))))
  }
  
  l_C_set <- length(C_set)
  l_P_set <- length(P_set)
  
  # Expanding the IPS iteratively to create the full IPS set for all other k
  for (i in 1:l_C_set) {
    upd <- TRUE
    new_c <- C_set[[i]]
    while(upd) {
      old_c <- new_c
      new_par <- setdiff(parents(old_c, L), old_c)
      new_c <- setdiff(children(new_par, L), new_par)
      if (setequal(new_c, old_c)) {
        upd <- FALSE
      } else {
        C_set <- append(C_set, list(new_c))
      }
    }
  }
  
  # Expanding the IFS iteratively to create the full IFS set for all other l
  for (i in 1:l_P_set) {
    upd <- TRUE
    new_p <- P_set[[i]]
    while(upd) {
      old_p <- new_p
      new_ch <- setdiff(children(old_p, L), old_p)
      new_p <- setdiff(parents(new_ch, L), new_ch)
      if (setequal(new_p, old_p)) {
        upd <- FALSE
      } else {
        P_set <- append(P_set, list(new_p))
      }
    }
  }
  
  # Explicitly add the empty set
  C_set <- c(C_set, list(character(0)))
  P_set <- c(P_set, list(character(0)))
  
  # Remove potential duplicates
  C_set <- C_set[!duplicated(C_set)]
  P_set <- P_set[!duplicated(P_set)]
  
  # Check potential IPS and IFS pairs
  for (i in seq_along(C_set)) {
    X_orig <- C_set[[i]]
    if (!empty_rec & !length(X_orig)) next
    for (j in seq_along(P_set)) {
      Y_orig <- P_set[[j]]
      if (!empty_emi & !length(Y_orig)) next
      if (!length(Y_orig) & !length(X_orig)) next
      an_de_XY_orig <- intersect(uu(an[Y_orig]), uu(de[X_orig]))
      if (length(Y_orig)) {
        # Restrict only if the other candidate is not the empty set
        X <- intersect(X_orig, an_de_XY_orig)
      } else {
        X <- X_orig
      }
      if (length(X_orig)) {
        # Restrict only if the other candidate is not the empty set
        Y <- intersect(Y_orig, an_de_XY_orig)
      } else {
        Y <- Y_orig
      }
      XY <- union(X, Y)
      n_X <- length(X)
      n_Y <- length(Y)
      n_XY <- length(XY)
      # Do not group the entire graph
      if (n_XY == n || n_XY == 0) next
      # If there are IPS candidates, they must have parents
      if (n_X && !all(sapply(pa[X], length))) next
      # If there are IFS candidates, they must have children
      if (n_Y && !all(sapply(ch[Y], length))) next
      if (n_XY) {
        # Construct A according to the formation conditions
        if (length(X) == 0) an_de_XY_orig <- unique(ancestors(Y, L))
        if (length(Y) == 0) an_de_XY_orig <- unique(descendants(X, L))
        g_ne <- induced_subgraph(L, an_de_XY_orig)
        A <- names(V(g_ne))
        # Check if A obeys the restriction rule and if it is a conduit,
        # add it to the list
        if (length(A) > 1 && all(A %in% restrict)) {
            tc <- c(tc, is_conduit(A, L))
          }    
        }
    }
  }
  unique(tc)
}

 is_conduit <- function(A, g) {
  
  #if (length(X == 0) & length(Y == 0)) return(NULL)
  
  X <- get_receivers(A, g)
  Y <- get_emitters(A, g)
  
  pa_X <- setdiff(parents(A, g), A)
  ch_Y <- setdiff(children(A, g), A)

  # Create new G*
  g_star <- modify_clustered(g, A, pa_X, ch_Y)
  
  # Find all descendants of the parents and ancestors of children in G*
  de_pa_X <- setNames(lapply(pa_X, function(x) descendants(x, g_star)), pa_X)
  an_ch_Y <- setNames(lapply(ch_Y, function(x) ancestors(x, g_star)), ch_Y)
  
  # Verify the conditions when Re or Em are empty sets
  if (length(X) == 0) {
    common_anc <- an_ch_Y[[1]]
    for (l in seq_along(an_ch_Y)) {
      common_anc <- intersect(common_anc, an_ch_Y[[l]])
      if (length(common_anc) == 0) return(NULL)
    }
  }
  if (length(Y) == 0) {
    common_de <- de_pa_X[[1]]
    for (i in seq_along(de_pa_X)) {
      common_de <- intersect(common_de, de_pa_X[[i]])
      if (length(common_de) == 0) return(NULL)
    }
  }
  
  # Check whether each parent is among the ancestors of each child,
  # i.e. whether there exists a directed path from each parent to 
  # each child.
  if (length(X) > 0 & length(Y) > 0) {
    for (p in pa_X) {
      for (c in ch_Y) {
        if (!(p %in% unlist(an_ch_Y[c]))) return(NULL)
      }
    }
  }
  
  # Confirm that there are no receivers that should not be emitters
  # and that there are no emitters that should not be receivers
  XY <- intersect(X, Y)
  ex_X <- setdiff(X, XY)
  ex_Y <- setdiff(Y, XY)
  X_leak <- setdiff(uu(children(ex_X, g)), A)
  Y_leak <- setdiff(uu(parents(ex_Y, g)), A)
  if (!length(X_leak) && !length(Y_leak)) {
    list(
      list(
        vertices = A,
        receivers = X,
        emitters = Y
      )
    )
  } else {
    NULL
  }
}

# Returns the set of receivers for a given set
get_receivers <- function(t, g) {
  rec <- c()
  outside_parents <- setdiff(parents_unsrt(t, g), t)
  for (i in 1:length(t)) {
    if (any(parents_unsrt(t[i], g) %in% outside_parents)) rec <- c(rec, t[i])
  }
  return(rec)
}

# Returns the set of emitters for a given set
get_emitters <- function(t, g) {
  emi <- c()
  outside_children <- setdiff(children_unsrt(t, g), t)
  for (i in 1:length(t)) {
    if (any(children_unsrt(t[i], g) %in% outside_children)) emi <- c(emi, t[i])
  }
  return(emi)
}

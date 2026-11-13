# Simulation: identification-invariant clusters found by FindConduits and
# FindTrClust in graphs with latent confounders
#
# Change d for different bidirected arc generating probabilities!

source("02_findconduits.R")
source("03_findtrcomp.R")
source("04_checkid.R")
load("networks_new.RData")

# --- Grab array index and settings from command line ---
args       <- commandArgs(trailingOnly = TRUE)
task_id    <- as.integer(args[1])
m          <- 1000
graph_name <- names(networks)[task_id]
g          <- networks[[task_id]]

set.seed(261008 + task_id)

message(sprintf("Task %d: running on graph '%s'",
                task_id, graph_name))


# --- Generating the latent confounders and the excluded set ---

# Adds a bidirected arc between each pair of observed vertices with probability
# p. Each arc is an explicit unobserved vertex with two children.
add_bidirected <- function(g) {
  obs <- igraph::V(g)$name
  pairs <- utils::combn(obs, 2)
  d <- 0
  p <- d / (length(obs))
  pick <- pairs[, runif(ncol(pairs)) < p, drop = FALSE]
  k <- ncol(pick)
  if (k == 0) return(list(g = g, unobs = character(0)))
  unobs <- sprintf("U_%d", seq_len(k))
  if (any(unobs %in% obs)) stop("Names of unobserved vertices clash with the graph.")
  el <- rbind(cbind(unobs, pick[1, ]), cbind(unobs, pick[2, ]))
  g  <- igraph::add_vertices(g, k, name = unobs)
  g  <- igraph::add_edges(g, as.vector(t(el)))
  list(g = g, unobs = unobs)
}

# Samples the family D: all leaf vertices (vertices without children) together
# with k = max(2, min(4, round(|V| / 25))) vertices drawn uniformly at random
# from the proper ancestors of the leaves. (In a DAG, every vertex that is not a
# leaf is an ancestor of some leaf.) If the leaves have fewer than k ancestors,
# all of them are included. `g` is the graph without the unobserved vertices.
sample_family <- function(g) {
  obs <- igraph::V(g)$name
  leaves <- obs[igraph::degree(g, mode = "out") == 0]
  anc <- setdiff(unique(unlist(lapply(leaves, function(v)
    names(igraph::subcomponent(g, v, mode = "in"))))), leaves)
  k <- min(max(2, min(4, round(length(obs) / 25))), length(anc))
  c(anc[sample.int(length(anc), k)], leaves)
}


# --- Certification ---

# Certifies a cluster with check_id(). Returns the result of check_id()
# ("plain", "c-isolated", "c-coherent", "not coherent" or "undetermined"), or
# "timeout" if check_id() timed out or errored.
certify <- function(g, clust, unobs, timeout_sec) {
  tryCatch(
    R.utils::withTimeout(
      check_id(g, clust, unobs, coherent = TRUE, check_c_coherence = TRUE),
      timeout = timeout_sec, onTimeout = "error"
    ),
    error = function(e) "timeout"
  )
}

# Certifies a list of candidates within a total time budget. Candidates left
# when the budget runs out are marked "unchecked".
certify_all <- function(g, cand, unobs, timeout_sec, budget_sec) {
  type <- rep("unchecked", length(cand))
  t_start <- proc.time()["elapsed"]
  for (i in seq_along(cand)) {
    if (proc.time()["elapsed"] - t_start > budget_sec) break
    type[i] <- certify(g, cand[[i]], unobs, timeout_sec)
  }
  type
}


# --- Clustering and graph size ---

# Greedy clustering: the largest cluster first, then the next largest that does
# not intersect the clusters chosen so far
greedy_clusters <- function(cl) {
  cl     <- cl[order(lengths(cl), decreasing = TRUE)]
  chosen <- list()
  used   <- character(0)
  for (x in cl) {
    if (!any(x %in% used)) {
      chosen <- c(chosen, list(x))
      used   <- c(used, x)
    }
  }
  chosen
}

# Clusters a set of disjoint clusters at once. An unobserved vertex with both
# children in the same cluster ends up with a single child (an error term of the
# cluster vertex).
cluster_graph <- function(g, clusters) {
  nodes <- igraph::V(g)$name
  map   <- stats::setNames(nodes, nodes)
  for (i in seq_along(clusters)) map[clusters[[i]]] <- sprintf("CL_%d", i)
  el <- igraph::as_edgelist(g)
  el <- unique(matrix(map[el], ncol = 2))
  el <- el[el[, 1] != el[, 2], , drop = FALSE]
  igraph::graph_from_data_frame(as.data.frame(el), directed = TRUE,
                                vertices = data.frame(name = unique(map)))
}

# Bidirected arcs as a two-column matrix of observed vertices (one row per arc,
# parallel arcs merged)
bidirected_arcs <- function(g, unobs) {
  el <- igraph::as_edgelist(g)
  el <- el[el[, 1] %in% unobs, , drop = FALSE]
  ch <- split(el[, 2], el[, 1])
  ch <- ch[lengths(ch) == 2]
  if (!length(ch)) return(matrix(character(0), 0, 2))
  arcs <- t(vapply(ch, sort, character(2)))
  unique(arcs)
}

# Numbers of observed vertices, directed edges and bidirected arcs
graph_size <- function(g, unobs) {
  el <- igraph::as_edgelist(g)
  c(n_obs = length(setdiff(igraph::V(g)$name, unobs)),
    n_dir = sum(!(el[, 1] %in% unobs)),
    n_bi  = nrow(bidirected_arcs(g, unobs)))
}

# --- Summaries ---

# Summary of one method: candidates, certification, greedy clustering and the
# reduced graph. `type` is the certification result of each candidate.
method_summary <- function(g, unobs, cand, type) {
  na_res <- c(candidates = NA, plain = NA, c_isolated = NA, c_coherent = NA,
              certified = NA, unchecked = NA, n_clusterings = NA, n_clustered = NA,
              reduced_obs = NA, reduced_dir = NA, reduced_bi = NA)
  if (is.null(cand)) return(na_res)
  certified <- cand[type %in% c("plain", "c-isolated", "c-coherent")]
  chosen    <- greedy_clusters(certified)
  g_cl      <- cluster_graph(g, chosen)
  size_cl   <- graph_size(g_cl, unobs)
  c(candidates    = length(cand),
    plain         = sum(type == "plain"),
    c_isolated    = sum(type == "c-isolated"),
    c_coherent    = sum(type == "c-coherent"),
    certified     = length(certified),
    unchecked     = sum(type %in% c("unchecked", "timeout")),
    n_clusterings = length(chosen),
    n_clustered   = length(unlist(chosen)),
    reduced_obs   = unname(size_cl["n_obs"]),
    reduced_dir   = unname(size_cl["n_dir"]),
    reduced_bi    = unname(size_cl["n_bi"]))
}

with_suffix <- function(x, suffix) as.list(stats::setNames(x, paste0(names(x), "_", suffix)))


# --- One repetition ---

run_replicate <- function(g, timeout_sec, cert_timeout_sec, cert_budget_sec) {
  lat   <- add_bidirected(g)
  g_u   <- lat$g
  unobs <- lat$unobs
  
  family <- sample_family(g)
  excl   <- c(family, unobs)
  
  size_g <- graph_size(g_u, unobs)
  
  # Time find_conduits
  t_start_cond <- proc.time()["elapsed"]
  cond <- tryCatch(
    R.utils::withTimeout(
      find_conduits(g_u, prohibit = excl, unobs = unobs),
      timeout = timeout_sec, onTimeout = "error"
    ),
    error = function(e) {
      message(sprintf("find_conduits timed out or errored: %s", conditionMessage(e)))
      NULL
    }
  )
  elapsed_cond <- proc.time()["elapsed"] - t_start_cond
  
  # Time find_trclust
  t_start_trclust <- proc.time()["elapsed"]
  trclust <- tryCatch(
    R.utils::withTimeout(
      {
        trcomp <- find_transit_components(g_u, prohibit = excl)
        lapply(find_transit_clusters(g_u, trcomp), function(x) x$vertices)
      },
      timeout = timeout_sec, onTimeout = "error"
    ),
    error = function(e) {
      message(sprintf("find_trclust timed out or errored: %s", conditionMessage(e)))
      NULL
    }
  )
  elapsed_trclust <- proc.time()["elapsed"] - t_start_trclust
  
  # Clusters may not contain family vertices or unobserved vertices. Neither
  # algorithm should return such clusters; the number dropped is recorded.
  drop_bad <- function(cl) {
    if (is.null(cl)) return(NULL)
    cl[!vapply(cl, function(x) any(x %in% excl), logical(1))]
  }
  dropped_cond    <- if (is.null(cond))    NA else length(cond)
  dropped_trclust <- if (is.null(trclust)) NA else length(trclust)
  cond    <- drop_bad(cond)
  trclust <- drop_bad(trclust)
  dropped_cond    <- dropped_cond    - if (is.null(cond))    0 else length(cond)
  dropped_trclust <- dropped_trclust - if (is.null(trclust)) 0 else length(trclust)
  
  # Certification
  t_start <- proc.time()["elapsed"]
  type_cond <- if (is.null(cond)) NULL else
    certify_all(g_u, cond, unobs, cert_timeout_sec, cert_budget_sec)
  cert_elapsed_cond <- proc.time()["elapsed"] - t_start
  
  t_start <- proc.time()["elapsed"]
  type_trclust <- if (is.null(trclust)) NULL else
    certify_all(g_u, trclust, unobs, cert_timeout_sec, cert_budget_sec)
  cert_elapsed_trclust <- proc.time()["elapsed"] - t_start
  
  c(
    list(
      n_nodes              = unname(size_g["n_obs"]),
      n_dir                = unname(size_g["n_dir"]),
      n_bi                 = unname(size_g["n_bi"]),
      n_family             = length(family),
      family               = paste(family, collapse = ", "),
      elapsed_cond         = unname(elapsed_cond),
      elapsed_trclust      = unname(elapsed_trclust),
      cert_elapsed_cond    = unname(cert_elapsed_cond),
      cert_elapsed_trclust = unname(cert_elapsed_trclust),
      dropped_cond         = dropped_cond,
      dropped_trclust      = dropped_trclust
    ),
    with_suffix(method_summary(g_u, unobs, cond,    type_cond),    "cond"),
    with_suffix(method_summary(g_u, unobs, trclust, type_trclust), "trclust")
  )
}


# --- Simulation ---

simulate_id <- function(g, m, update_every = NULL, timeout_sec = 900,
                        cert_timeout_sec = 60, cert_budget_sec = 1800,
                        save_every = 10, save_prefix = "results") {
  results <- vector("list", m)
  for (iter in seq_len(m)) {
    results[[iter]] <- c(list(rep = iter),
                         run_replicate(g, timeout_sec, cert_timeout_sec, cert_budget_sec))
    
    # Periodic checkpoint save (overwritten each time)
    if (iter %% save_every == 0) {
      checkpoint <- dplyr::bind_rows(results[seq_len(iter)])
      save_path  <- sprintf("%s_checkpoint.RData", save_prefix)
      save(checkpoint, file = save_path)
      message(sprintf("Checkpoint saved: %s (%d repetitions)", save_path, iter))
    }
    if (!is.null(update_every) && iter %% update_every == 0) {
      message(sprintf("Iteration %d / %d complete", iter, m))
    }
  }
  dplyr::bind_rows(results)
}


results <- simulate_id(
  g,
  m                = m,
  update_every     = 10,
  timeout_sec      = 1800,
  cert_timeout_sec = 60,
  cert_budget_sec  = 1800,
  save_every       = 10,
  save_prefix      = sprintf("/folder_name_here", graph_name)
)

# Final save
save(results, file = sprintf("/folder_name_here", graph_name))
message(sprintf("Final results saved for graph '%s'", graph_name))
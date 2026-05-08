# Running the simulation of Section 6.1

source("02_findtrclust.R")
source("03_findtrcomp.R")
load("networks.RData")

# --- Grab array index from command line ---
args       <- commandArgs(trailingOnly = TRUE)
task_id    <- as.integer(args[1])
graph_name <- names(networks)[task_id]
g          <- networks[[task_id]]

set.seed(260504 + task_id)

message(sprintf("Task %d: running on graph '%s'", task_id, graph_name))

clustering_results <- function(g, cand) {
  clustered  <- c()
  clustering <- 0
  candidates <- length(cand)
  if (candidates == 0) return(c(0, 0, 0))
  
  sorted <- cand[order(
    sapply(cand, function(x) length(x$vertices)),
    decreasing = TRUE
  )]
  
  clustered  <- sorted[[1]]$vertices
  clustering <- 1
  
  if (candidates > 1) {
    for (i in 2:length(sorted)) {
      if (!length(intersect(clustered, sorted[[i]]$vertices))) {
        clustered  <- c(clustered, sorted[[i]]$vertices)
        clustering <- clustering + 1
      }
    }
  }
  return(c(length(clustered), clustering, candidates))
}

simulate_clustering <- function(g, n, update_every = NULL, timeout_sec = 900,
                                save_every = 10, save_prefix = "results") {
  nodes   <- V(g)
  n_nodes <- length(nodes)
  results <- vector("list", n)
  
  for (iter in seq_len(n)) {
    rand_x <- NULL
    for (i in seq_len(1000)) {
      rand_y <- sample(nodes, 1)
      rand_y <- names(rand_y)
      anc_y  <- setdiff(ancestors(rand_y, g), rand_y)
      if (!length(anc_y)) next
      rand_x <- sample(anc_y, 1)
      break
    }
    if (is.null(rand_x)) stop("Failed to find a valid rand_x after 1000 attempts")
    
    prohibit <- c(rand_x, rand_y)
    
    # Time find_conduits
    t_start_cond <- proc.time()["elapsed"]
    cond <- tryCatch(
      R.utils::withTimeout(
        find_conduits(g, prohibit = prohibit),
        timeout   = timeout_sec,
        onTimeout = "error"
      ),
      error = function(e) {
        message(sprintf("find_conduits timed out or errored on iteration %d: %s", iter, conditionMessage(e)))
        NULL
      }
    )
    elapsed_cond <- proc.time()["elapsed"] - t_start_cond
    
    # Time find_trclust
    t_start_trclust <- proc.time()["elapsed"]
    trclust <- tryCatch(
      R.utils::withTimeout(
        {
          trcomp <- find_transit_components(g, prohibit = prohibit)
          find_transit_clusters(g, trcomp)
        },
        timeout   = timeout_sec,
        onTimeout = "error"
      ),
      error = function(e) {
        message(sprintf("find_trclust timed out or errored on iteration %d: %s", iter, conditionMessage(e)))
        NULL
      }
    )
    elapsed_trclust <- proc.time()["elapsed"] - t_start_trclust
    
    findcond_res    <- if (!is.null(cond))    clustering_results(g, cond)    else c(NA, NA, NA)
    findtrclust_res <- if (!is.null(trclust)) clustering_results(g, trclust) else c(NA, NA, NA)
    
    results[[iter]] <- list(
      n_nodes                = n_nodes,
      prohibit               = paste0(prohibit[1], ", ", prohibit[2]),
      n_clustered_cond       = findcond_res[1],
      n_clusterings_cond     = findcond_res[2],
      reduced_size_cond      = n_nodes - findcond_res[1] + findcond_res[2],
      candidates_cond        = findcond_res[3],
      elapsed_cond           = elapsed_cond,
      n_clustered_trclust    = findtrclust_res[1],
      n_clusterings_trclust  = findtrclust_res[2],
      reduced_size_trclust   = n_nodes - findtrclust_res[1] + findtrclust_res[2],
      candidates_trclust     = findtrclust_res[3],
      elapsed_trclust        = elapsed_trclust
    )
    
    # Periodic checkpoint save
    if (iter %% save_every == 0) {
      checkpoint <- dplyr::bind_rows(results[seq_len(iter)])
      save_path  <- sprintf("%s_checkpoint.RData", save_prefix)
      save(checkpoint, file = save_path)
      message(sprintf("Checkpoint saved: %s", save_path))
    }
    
    if (!is.null(update_every) && iter %% update_every == 0) {
      message(sprintf("Iteration %d / %d complete", iter, n))
    }
  }
  
  dplyr::bind_rows(results)
}

results <- simulate_clustering(
  g,
  n            = 1000,
  update_every = 10,
  timeout_sec  = 900,
  save_every   = 10,
  save_prefix  = sprintf("/scratch/jkarvane/conduit/%s", graph_name)
)

# Final save
save(results, file = sprintf("/scratch/jkarvane/conduit/%s_final.RData", graph_name))
message(sprintf("Final results saved for graph '%s'", graph_name))

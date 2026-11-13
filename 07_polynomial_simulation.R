# Simulation: running time of FindConduits with respect to graph size

source("06_generate_graph.R")

# --- Settings ---

sizes         <- seq(25, 150, by = 25)
densities     <- c(sparse = 2, dense = 4) 
n_tasks       <- 100                       # number of array tasks
reps_per_task <- 10                        # passes per task
# Repetitions per size and density: n_tasks * reps_per_task = 1000
timeout_sec   <- 1800                      # time limit of one FindConduits run
out_dir       <- "/folder_name_here"

combos <- expand.grid(n = sizes, d = unname(densities))

# --- Grab array index from command line ---

args    <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
if (is.na(task_id) || task_id < 1 || task_id > n_tasks) {
  stop("task_id must be between 1 and ", n_tasks)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save_prefix     <- file.path(out_dir, sprintf("timing_task%03d", task_id))
checkpoint_file <- sprintf("%s_checkpoint.RData", save_prefix)
node            <- unname(Sys.info()["nodename"])

message(sprintf("Task %d on node %s", task_id, node))


# --- One repetition ---

run_instance <- function(n, d, timeout_sec) {
  g <- random_dag(n, d)
  t0 <- proc.time()["elapsed"]
  cond <- tryCatch(
    R.utils::withTimeout(find_conduits(g),
                         timeout = timeout_sec, onTimeout = "error"),
    error = function(e) {
      message(sprintf("find_conduits timed out or errored: %s", conditionMessage(e)))
      NULL
    })
  elapsed <- unname(proc.time()["elapsed"] - t0)
  data.frame(n = n, d = d, m = igraph::ecount(g),
             timed_out = is.null(cond),
             time      = if (is.null(cond)) NA else elapsed,
             conduits  = if (is.null(cond)) NA else length(cond))
}


# --- Simulation ---

# Resume from the checkpoint if the task was interrupted earlier. The
# checkpoint only contains complete passes.
if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  results    <- checkpoint
  first_pass <- max(results$pass) + 1
  message(sprintf("Resuming from checkpoint: %d passes done", first_pass - 1))
} else {
  results    <- NULL
  first_pass <- 1
}

# Warm-up run, not recorded, so that the first recorded instance does not
# include the one-off cost of loading packages and compiling functions
invisible(find_conduits(random_dag(10, densities[1])))

if (first_pass <= reps_per_task) {
  for (pass in first_pass:reps_per_task) {
    # A seed per task and pass, so that a resumed task produces the same graphs
    set.seed(261102 + 1000 * task_id + pass)
    ord <- sample.int(nrow(combos))
    
    res_pass <- lapply(seq_along(ord), function(k) {
      data.frame(task = task_id, pass = pass, position = k, node = node,
                 run_instance(combos$n[ord[k]], combos$d[ord[k]], timeout_sec))
    })
    results <- dplyr::bind_rows(results, res_pass)
    
    # Checkpoint after every complete pass (overwritten each time)
    checkpoint <- results
    save(checkpoint, file = checkpoint_file)
    message(sprintf("Pass %d / %d complete", pass, reps_per_task))
  }
}

# Final save
save(results, file = sprintf("%s_final.RData", save_prefix))
message(sprintf("Final results saved for task %d", task_id))
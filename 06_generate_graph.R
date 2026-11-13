source("01_utils.R")
source("02_findconduits.R")

# Generates a random igraph DAG
# Inputs:
# n = graph size
# if d = 4 --> dense graphs. Otherwise edge generating probability is d / n
random_dag <- function(n, d = 3) {
  
  # Create a pool of variable names and randomly sample n of them
  variable_pool <- paste0("z", 1:n)
  vars <- sample(variable_pool, n)
  
  p <- d / n
  
  if (d == 4) p <- 0.4
  
  # Force the last variable to be "y" (outcome variable)
  vars[n] <- "y"
  
  # Randomly select position for "x" (treatment variable) between 2 and n-1
  x_ind <- sample(2:(n - 1), 1)
  vars[x_ind] <- "x"
  
  # Initialize edge vectors (v = from, w = to)
  v <- character(0)
  w <- character(0)
  
  # Phase 1: Generate edges from x_ind to n with probability p
  # Working backwards from y (position n) down to x (position x_ind)
  for (i in (n - 1):x_ind) {
    j <- n
    while(j > i) {
      u <- runif(1)
      if (u < p) {
        # Add edge from vars[i] to vars[j]
        v <- c(v, vars[i])
        w <- c(w, vars[j])
      }
      j <- j - 1
    }
  } 
  
  # Phase 2: Generate edges from position 1 to x_ind-1 with probability p
  # These are variables that come before x in topological order
  for (i in (x_ind - 1):1) {
    j <- n
    while(j > i) {
      u <- runif(1)
      if (u < p) {
        # Add edge from vars[i] to vars[j]
        v <- c(v, vars[i])
        w <- c(w, vars[j])
      }
      j <- j - 1
    }
  } 
  
  g0 <- graph_from_data_frame(data.frame(from = v, to = w), directed = TRUE,
                              vertices = data.frame(name = vars))
  comps <- split(seq_len(n), components(g0, mode = "weak")$membership)
  if (length(comps) > 1) {
    for (k in 2:length(comps)) {
      a <- comps[[k]][sample.int(length(comps[[k]]), 1)]
      prev <- unlist(comps[1:(k - 1)])
      b <- prev[sample.int(length(prev), 1)]
      v <- c(v, vars[min(a, b)])
      w <- c(w, vars[max(a, b)])
    }
  }
  
  # Create the DAG from collected edges
  edges <- cbind(v, w)
  dag <- graph_from_edgelist(edges, directed = TRUE)
  
  return(dag)
}

# Running the insurance example of Section 4.3

library(bnlearn)
library(igraph)
library(causaleffect)

source("02_findconduits.R")
source("03_findtrcomp.R")

# Function that adds an unobserved confounder between two variables
# in an igraph object

add_confounder <- function(g, var1, var2, unobs = character(0)) {
  node_names <- V(g)$name
  missing <- setdiff(c(var1, var2), node_names)
  if (length(missing) > 0) {
    stop("Variable(s) not found in graph: ", paste(missing, collapse = ", "))
  }
  if (var1 == var2) stop("The two variables must be different.")
  if (any(c(var1, var2) %in% unobs)) stop("Confounders can only be added between observed variables.")
  
  # Name the new unobserved vertex U_k with the first k not already in use
  k <- 1
  while (sprintf("U_%d", k) %in% node_names) k <- k + 1
  u <- sprintf("U_%d", k)
  
  # Add the bidirected arc as an explicit unobserved vertex with two children
  g <- add_vertices(g, 1, name = u)
  g <- add_edges(g, c(u, var1, u, var2))
  
  g
}
#############
# INSURANCE #
#############

load("networks_new.RData")
g <- networks$insurance
plot(g)


# Adding confounders present in the graph and pruning non-ancestors of the
# responses
g <- add_confounder(g, "RiskAversion", "ThisCarCost")
g <- add_confounder(g, "RiskAversion", "SocioEcon")
g_prune <- induced_subgraph(g, ancestors(c("ThisCarCost", "OtherCarCost"), g))

# case i) RiskAversion as treatment

# Discovering the conduits
con_restricted <- find_conduits(g_prune, prohibit = c("ThisCarCost", "OtherCarCost", "RiskAversion", "SocioEcon"), unobs = c("U_1", "U_2"))

con <- find_conduits(g_prune, unobs = c("U_1", "U_2"))

# Constructing the clustered graph

add_bidirected_causaleffect <- function(g, var1, var2) {
  missing <- setdiff(c(var1, var2), igraph::V(g)$name)
  if (length(missing) > 0) {
    stop("Variable(s) not found in graph: ", paste(missing, collapse = ", "))
  }
  if (var1 == var2) stop("The two variables must be different.")
  igraph::add_edges(g, c(var1, var2, var2, var1), description = "U")
}

g_clustered <- graph_from_edgelist(rbind(c("a", "m"), c("a", "s"), c("s", "m"), c("s", "r"), c("a", "r"),
                                         c("r", "m"), c("m", "o"), c("m", "t"), c("e", "m")))

g_clustered <- add_bidirected_causaleffect(g_clustered, "r", "t")
g_clustered <- add_bidirected_causaleffect(g_clustered, "r", "s")

# Effects of risk aversion

causaleffect::causal.effect(c("t", "o"), "r", G = g_clustered, prune = TRUE)
causaleffect::causal.effect(c("t"), "r", G = g_clustered, prune = TRUE)
causaleffect::causal.effect(c("o"), "r", G = g_clustered, prune = TRUE)

# Effects of socioeconomic status

causaleffect::causal.effect(c("t", "o"), "s", G = g_clustered, prune = TRUE)
causaleffect::causal.effect(c("t"), "s", G = g_clustered, prune = TRUE)
causaleffect::causal.effect(c("o"), "s", G = g_clustered, prune = TRUE)

# Joint effect of risk aversion and socioeconomic status

causaleffect::causal.effect(c("t", "o"), c("r", "s"), G = g_clustered, prune = TRUE)

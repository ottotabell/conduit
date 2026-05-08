# Running the insurance example of Section 6.2

library(bnlearn)
library(igraph)
library(stringr)
library(causaleffect)

source("02_findconduits.R")
source("03_findtrcomp.R")

# Function that adds an unobserved confounder between two variables
# in an igraph object

add_confounder <- function(g, var1, var2, confounder_name = NULL) {
  node_names <- V(g)$name
  missing <- setdiff(c(var1, var2), node_names)
  if (length(missing) > 0) {
    stop("Variable(s) not found in graph: ", paste(missing, collapse = ", "))
  }
  
  # Add bidirected edge as two opposite directed edges with description = "U"
  g <- add_edges(
    g,
    c(var1, var2, var2, var1),
    description = c("U", "U")
  )
  
  return(g)
}

#############
# INSURANCE #
#############

load("networks.RData")
g <- networks$insurance
plot(g)


# Adding confounders present in the graph and pruning non-ancestors of the
# responses
g <- add_confounder(g, "RiskAversion", "ThisCarCost", "U2")
g <- add_confounder(g, "RiskAversion", "SocioEcon", "U3")
g_prune <- induced_subgraph(g, ancestors(c("ThisCarCost", "OtherCarCost"), g))

# case i) RiskAversion as treatment

# Discovering the conduits
con <- find_conduits(g_prune, prohibit = c("ThisCarCost", "OtherCarCost", "RiskAversion"))

# Causal effect of the original graph
causaleffect::causal.effect(c("ThisCarCost", "OtherCarCost"), "RiskAversion", G = g_prune, prune = T)

# Constructing the clustered graph
g_clustered <- graph_from_edgelist(rbind(c("a", "r"), c("a", "m"), c("a", "s"), c("s", "r"), c("s", "m"),
                                     c("r", "m"), c("m", "o"), c("m", "t")))
g_clustered <- add_confounder(g_clustered, "r", "t")
g_clustered <- add_confounder(g_clustered, "r", "s")

# Causal effect for the clustered graph
causaleffect::causal.effect(c("t", "o"), "r", G = g_clustered, prune = T)

# case ii) SocioEcon as treatment

con2 <- find_conduits(g_prune, prohibit = c("ThisCarCost", "OtherCarCost", "SocioEcon"))
con2

g_clustered <- graph_from_edgelist(rbind(c("a", "m"), c("a", "s"), c("s", "m"),
                                         c("m", "o"), c("m", "t")))

g_clustered <- add_confounder(g_clustered, "m", "t")
g_clustered <- add_confounder(g_clustered, "m", "s")

causaleffect::causal.effect(c("ThisCarCost", "OtherCarCost"), "SocioEcon", G = g_prune, prune = T)
causaleffect::causal.effect(c("t", "o"), "s", G = g_clustered, prune = T)
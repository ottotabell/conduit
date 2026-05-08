# Analyzing the results of the simulation study

# files <- list.files("results", pattern = "_final\\.RData$", full.names = TRUE)

library(tidyverse)

all_results <- dplyr::bind_rows(lapply(files, function(f) {
  e <- new.env()
  load(f, envir = e)
  graph_name <- sub("_final\\.RData$", "", basename(f))
  dplyr::mutate(e$results, graph = graph_name)
}))

condres <- all_results %>% group_by(graph) %>% 
       summarise(
         size            = median(n_nodes, na.rm = T),
         reduced         = median(n_nodes - reduced_size_cond, na.rm = T),
         reduced_q1      = quantile(n_nodes - reduced_size_cond, 0.25, na.rm = T),
         reduced_q3      = quantile(n_nodes - reduced_size_cond, 0.75, na.rm = T),
         n_conduits      = median(candidates_cond, na.rm = T),
         n_conduits_q1   = quantile(candidates_cond, 0.25, na.rm = T),
         n_conduits_q3   = quantile(candidates_cond, 0.75, na.rm = T),
         time            = median(elapsed_cond, na.rm = T),
         time_q1         = quantile(elapsed_cond, 0.25, na.rm = T),
         time_q3         = quantile(elapsed_cond, 0.75, na.rm = T)
       )

trclustres <- all_results %>% group_by(graph) %>% 
  summarise(
    size            = median(n_nodes, na.rm = T),
    reduced         = median(n_nodes - reduced_size_trclust, na.rm = T),
    reduced_q1      = quantile(n_nodes - reduced_size_trclust, 0.25, na.rm = T),
    reduced_q3      = quantile(n_nodes - reduced_size_trclust, 0.75, na.rm = T),
    n_conduits      = median(candidates_trclust, na.rm = T),
    n_conduits_q1   = quantile(candidates_trclust, 0.25, na.rm = T),
    n_conduits_q3   = quantile(candidates_trclust, 0.75, na.rm = T),
    time            = median(elapsed_trclust, na.rm = T),
    time_q1         = quantile(elapsed_trclust, 0.25, na.rm = T),
    time_q3         = quantile(elapsed_trclust, 0.75, na.rm = T)
  )

Supplementary code for the paper
"Efficient and Practical Clustering of Variables in Causal Graphs"
===========================================================================
 
This folder contains an R implementation of the FindConduits algorithm
(Algorithm 1 of the paper), the algorithms it is compared with, and the code
for the simulations and the example of Section 4.
  
FILES
-----
 
  01_utils.R
      Helper functions for igraph graphs: parents, children, ancestors and
      descendants of vertex sets, and the path check for condition (A) of
      Definition 6 (check_paths). Sourced by the other files.
 
  02_findconduits.R
      find_conduits(): the FindConduits algorithm (Algorithm 1).
 
  03_findtrcomp.R
      find_transit_components() and find_transit_clusters(): the algorithm
      of Tikka et al. (2023) for finding transit clusters (FindTrClust in
      the paper).
 
  04_checkid.R
      check_id(): checks the sufficient conditions of Nykänen et al. (2026)
      for universal identification invariance (Appendix A).
 
  05_simulation00.R
      Simulation of Section 4.2 and Appendix E: FindConduits and FindTrClust
      in graphs of the Bayesian Network Repository with random bidirected
      edges.
 
  06_generate_graph.R
      random_dag(): random connected DAGs for Section 4.1.
 
  07_polynomial_simulation.R
      Simulation of Section 4.1: running time of FindConduits and the number
      of conduits in random DAGs.
 
  analyze_timing.R
      Summarizes the results of 07_polynomial_simulation.R and draws
      Figure 2.
 
  analyze_results.R
      Summarizes the results of 05_simulation00.R and prints the LaTeX code
      of Tables 1-4.
 
  insurance.R
      The insurance claim costs example of Section 4.3 (Figure 3).
 
  networks_new.RData
      A named list `networks` of DAGs from the Bayesian Network Repository
      (https://www.bnlearn.com/bnrepository/) as igraph objects.
 
  
REFERENCES
----------
 
Nykänen et al. (2026). Invariance of clustering operations in causal effect
identification.
 
Tikka et al. (2023). Clustering and structural robustness in causal
diagrams. Journal of Machine Learning Research, 24(195):1-32.
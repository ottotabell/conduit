# Analyzing the results of the simulation study

library(tidyverse)

# Folders of the results for each density s of the latent confounders
result_dirs <- c("0" = "results_00", "0.1" = "results_01", "0.2" = "results_20")

# How each value of s is described in the tables
arc_label <- c("0" = "none", "0.1" = "few", "0.2" = "moderately")
arc_text  <- function(s_value) arc_label[[sprintf("%g", s_value)]]

# Number of largest graphs (by |V|) shown in the main-text table
n_largest <- 5

# --- Reading the results ---

read_results <- function(dir, s) {
  files <- list.files(dir, pattern = "_final\\.RData$", full.names = TRUE)
  if (!length(files)) warning("No result files found in ", dir)
  dplyr::bind_rows(lapply(files, function(f) {
    e <- new.env()
    load(f, envir = e)
    graph_name <- sub("_final\\.RData$", "", basename(f))
    dplyr::mutate(e$results, graph = graph_name, s = s)
  }))
}

all_results <- dplyr::bind_rows(lapply(names(result_dirs), function(s)
  read_results(result_dirs[[s]], as.numeric(s))))

all_results$graph <- substr(all_results$graph, 1, nchar(all_results$graph) - 3)

all_results$graph <- ifelse(all_results$graph == "mehra-original", "mehra", all_results$graph)

all_results <- all_results %>% filter(graph != "survey")

# --- One row per graph, s, algorithm and instance ---

# Instances where the algorithm hit the time limit have no results; their
# running times are excluded as well.
long <- bind_rows(
  all_results %>% transmute(
    graph, s, n_nodes, r_size = n_nodes - n_family, method = "cond",
    missing   = is.na(candidates_cond),
    reduction = n_nodes - reduced_obs_cond,
    all_clusters = candidates_cond,
    clusters  = certified_cond,
    time      = ifelse(missing, NA, elapsed_cond)
  ),
  all_results %>% transmute(
    graph, s, n_nodes, r_size = n_nodes - n_family, method = "trclust",
    missing   = is.na(candidates_trclust),
    reduction = n_nodes - reduced_obs_trclust,
    all_clusters = candidates_trclust,
    clusters  = certified_trclust,
    time      = ifelse(missing, NA, elapsed_trclust)
  )
)

# Median and quartiles of each variable per graph, s and algorithm
res <- long %>%
  group_by(graph, s, method) %>%
  summarise(
    n         = n(),
    v_size    = median(n_nodes),
    size      = median(r_size),
    size_min  = min(r_size),
    size_max  = max(r_size),
    n_missing = sum(missing),
    across(c(reduction, all_clusters, clusters, time),
           list(med = ~ median(.x, na.rm = TRUE),
                q1  = ~ quantile(.x, 0.25, na.rm = TRUE, names = FALSE),
                q3  = ~ quantile(.x, 0.75, na.rm = TRUE, names = FALSE))),
    .groups = "drop"
  ) %>%
  arrange(graph, s, factor(method, levels = c("cond", "trclust")))


# --- Helpers for the LaTeX tables ---

# Format "median (Q1, Q3)" strings; "--" if there are no results
fmt <- function(med, q1, q3, digits = 1) {
  if (is.na(med)) return("--")
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
          med, q1, q3)
}

alg_name <- c(cond = "\\FindConduits{}", trclust = "\\FindTrClustA{}")

# Number of instances (repetitions), stated once in the caption
instances_text <- function(res) {
  n_rep <- unique(res$n)
  if (length(n_rep) > 1) warning("Graphs have different numbers of instances: ",
                                 paste(sort(n_rep), collapse = ", "))
  if (length(n_rep) == 1) n_rep else paste(range(n_rep), collapse = "--")
}

# The size of the restriction set R should be the same in every instance of a
# graph; warn if it varies
check_r_size <- function(res) {
  varies <- unique(res$graph[res$size_min != res$size_max])
  if (length(varies)) warning("The size of R varies between instances in: ",
                              paste(varies, collapse = ", "), "; the median is shown.")
}

# Footnote mark and text for instances that exceeded the time limit
missing_note <- function(r, notes, time_limit_min) {
  if (r$n_missing == 0) return(list(mark = "", notes = notes))
  letter <- letters[length(notes) + 1]
  list(
    mark  = sprintf("\\textsuperscript{%s}", letter),
    notes = c(notes, sprintf(
      "\\textsuperscript{%s} %d instance%s exceeded the %d-minute time limit and %s excluded from the reported statistics.",
      letter, r$n_missing, ifelse(r$n_missing == 1, "", "s"), time_limit_min,
      ifelse(r$n_missing == 1, "was", "were")))
  )
}

# Graph cell: name, |V| and |R|, on three lines or, if `one_line` is TRUE,
# with |V| and |R| on the same line
graph_cell <- function(gr, sub, n_rows, one_line = FALSE) {
  sep <- if (one_line) ", " else "\\\\"
  sprintf("\\multirow{%d}{*}{\\makecell[l]{%s\\\\$|\\+V| = %d$%s$|\\+R \\cap \\+V| = %d$}}",
          n_rows, str_to_upper(gr), round(sub$v_size[1]), sep, round(sub$size[1]))
}

# Cells in the order all clusters, invariant clusters, reduction and, if
# `time` is TRUE, running time
cells <- function(r, time = FALSE) {
  x <- c(fmt(r$all_clusters_med, r$all_clusters_q1, r$all_clusters_q3, digits = 0),
         fmt(r$clusters_med,     r$clusters_q1,     r$clusters_q3,     digits = 0),
         fmt(r$reduction_med,    r$reduction_q1,    r$reduction_q3,    digits = 0))
  if (time) x <- c(x, fmt(r$time_med, r$time_q1, r$time_q3, digits = 1))
  x
}

# Column headers in the same order
header_cols <- function(time = FALSE) {
  paste(
    "& \\makecell{\\textbf{All}\\\\\\textbf{clusters}}",
    "& \\makecell{\\textbf{Invariant}\\\\\\textbf{clusters}}",
    "& \\textbf{Reduction}",
    if (time) "& \\textbf{Time (s)}" else "",
    "\\\\")
}

# Description of the columns for the captions
caption_common <- function(time = FALSE) {
  paste(
    "The table reports the median number of all clusters found,",
    "the median number of clusters among them that are universally identification invariant according to \\citet{nykanen2026invariance},",
    if (time) "the median reduction in the number of observed vertices obtained by clustering the invariant clusters, and the median running time of the algorithm." else
      "and the median reduction in the number of observed vertices obtained by clustering the invariant clusters.",
    "Interquartile ranges (Q1, Q3) are shown in parentheses.")
}

# Description of the densities of the bidirected arcs for the captions
arcs_caption <- "Each pair of vertices was connected by a bidirected edge with probability $s/|\\+V|$, where $s = 0$ (no edges), $s = 0.1$ (few edges) or $s = 0.2$ (moderately many edges)."


# --- Appendix tables: all graphs, one table per s ---

generate_latex_table <- function(res, s_value, time_limit_min = 30) {
  
  check_r_size(res)
  notes <- character(0)
  rows  <- character(0)
  
  for (gr in unique(res$graph)) {
    sub <- res %>% filter(graph == gr)
    rows <- c(rows, graph_cell(gr, sub, 2, one_line = TRUE))
    for (i in seq_len(nrow(sub))) {
      r  <- sub[i, ]
      mn <- missing_note(r, notes, time_limit_min)
      notes <- mn$notes
      x <- cells(r, time = TRUE)
      rows <- c(rows, sprintf("& %s%s & %s \\\\",
                              alg_name[[r$method]], mn$mark, paste(x, collapse = " & ")))
    }
    rows <- c(rows, "\\hline", "")
  }
  
  header <- c(
    "\\begin{table}[!b]",
    "\\centering",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{|l|lcccc|}",
    "\\hline",
    paste("\\textbf{Graph} & \\textbf{Algorithm}", header_cols(time = TRUE)),
    "\\hline",
    ""
  )
  
  caption <- paste(
    sprintf("\\caption{Comparison of \\FindConduits{} and \\FindTrClustA{} in all graphs with %s bidirected edges ($s = %g$).",
            c(none = "no", few = "few", moderately = "moderately many")[[arc_text(s_value)]], s_value),
    sprintf("Each graph was studied in %s instances (i.e., repetitions of the simulation), and each graph name is followed by the number of observed vertices $|\\+V|$ and the number of observed vertices in the restriction set $|\\+R \\cap \\+V|$.",
            instances_text(res)),
    caption_common(time = TRUE),
    paste(notes, collapse = " "))
  caption <- paste0(trimws(caption), "}")
  
  footer <- c(
    "\\end{tabular}",
    caption,
    sprintf("\\label{tab_results_s%s}", gsub("\\.", "", sprintf("%g", s_value))),
    "\\end{table}"
  )
  
  cat(paste(c(header, rows, footer), collapse = "\n"), "\n\n")
}

for (s_value in sort(unique(res$s))) {
  generate_latex_table(res %>% filter(s == s_value), s_value)
}


# --- Main-text table: the largest graphs, all values of s ---

generate_combined_table <- function(res, graphs, time_limit_min = 30) {
  
  check_r_size(res)
  notes    <- character(0)
  rows     <- character(0)
  s_values <- sort(unique(res$s))
  methods  <- c("cond", "trclust")
  n_rows   <- length(s_values) * length(methods)
  
  for (gr in graphs) {
    sub <- res %>% filter(graph == gr)
    for (j in seq_along(s_values)) {
      for (k in seq_along(methods)) {
        r <- sub %>% filter(s == s_values[j], method == methods[k])
        # Graph cell on the first row of the graph, s cell on the first row of
        # each value of s
        g_cell <- if (j == 1 && k == 1) graph_cell(gr, sub, n_rows) else ""
        s_cell <- if (k == 1) sprintf("\\multirow{%d}{*}{%s}", length(methods), arc_text(s_values[j])) else ""
        if (!nrow(r)) {
          x <- c("--", "--", "--"); mark <- ""
        } else {
          mn    <- missing_note(r, notes, time_limit_min)
          notes <- mn$notes
          mark  <- mn$mark
          x <- cells(r)
        }
        rows <- c(rows, sprintf("%s & %s & %s%s & %s \\\\", g_cell, s_cell,
                                alg_name[[methods[k]]], mark, paste(x, collapse = " & ")))
      }
      # Line between the values of s, and below the graph after the last one
      if (j < length(s_values)) rows <- c(rows, "\\cline{2-6}")
    }
    rows <- c(rows, "\\hline", "")
  }
  
  header <- c(
    "\\begin{table}[!t]",
    "\\footnotesize",
    "\\centering",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{|l|c|lccc|}",
    "\\hline",
    paste("\\textbf{Graph} & \\makecell{\\textbf{Bidirected}\\\\\\textbf{edges}} & \\textbf{Algorithm}", header_cols()),
    "\\hline",
    ""
  )
  
  caption <- paste(
    sprintf("\\caption{\\small Comparison of \\FindConduits{} and \\FindTrClustA{} in the %d largest graphs with different numbers of bidirected edges (none, few or moderately relative to the graph size; see Appendix~\\ref{sec_results} for details).", length(graphs)),
    sprintf("Each graph was studied in %s instances for each number of edges, and each graph name is followed by the number of observed vertices $|\\+V|$ and the number of observed vertices in the restriction set $|\\+R \\cap \\+V|$.",
            instances_text(res)),
    caption_common(),
    "",
    paste(notes, collapse = " "))
  caption <- paste0(trimws(caption), "}")
  
  footer <- c(
    "\\end{tabular}",
    caption,
    "\\label{tab_results}",
    "\\end{table}"
  )
  
  cat(paste(c(header, rows, footer), collapse = "\n"), "\n\n")
}

# The largest graphs by the number of vertices, from largest to smallest
largest <- res %>%
  distinct(graph, v_size) %>%
  arrange(desc(v_size)) %>%
  slice_head(n = n_largest) %>%
  pull(graph)

generate_combined_table(res %>% filter(graph %in% largest), largest)


# --- Maximum running times, for reporting in the text ---

long %>%
  group_by(graph, s, method) %>%
  summarise(max_time = if (all(is.na(time))) NA else max(time, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = max_time, names_prefix = "max_time_") %>%
  print(n = Inf)

all_results %>% group_by(graph, s) %>% summarise(
  any_reduction_cond = mean(n_clustered_cond > 0, na.rm = TRUE),
  cond_larger        = mean(reduced_obs_cond <  reduced_obs_trclust, na.rm = TRUE),
  tie                = mean(reduced_obs_cond == reduced_obs_trclust, na.rm = TRUE),
  cert_rate_cond     = sum(certified_cond, na.rm = TRUE) / sum(candidates_cond, na.rm = TRUE),
  plain_share_cond   = sum(plain_cond, na.rm = TRUE) / sum(certified_cond, na.rm = TRUE),
  cert_rate_trclust  = sum(certified_trclust, na.rm = TRUE) / sum(candidates_trclust, na.rm = TRUE),
  .groups = "drop"
) %>% print(n = Inf)


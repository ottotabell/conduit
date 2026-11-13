# Analyzing the running times of FindConduits (simulation 07_timing.R)

library(tidyverse)

# Folder of the result files of 07_timing.R
result_dir <- "timing"

# How each value of d is described in the table and the figure
density_label <- c("2" = "sparse", "4" = "dense")
# The value of d for a label
density_label_value <- function(lbl) names(density_label)[density_label == lbl]

# --- Reading the results ---

files <- list.files(result_dir, pattern = "_final\\.RData$", full.names = TRUE)
if (!length(files)) stop("No result files found in ", result_dir)
timing <- bind_rows(lapply(files, function(f) {
  e <- new.env()
  load(f, envir = e)
  e$results
})) %>%
  mutate(density = factor(density_label[as.character(d)], levels = density_label))

# --- Summary per graph size and density ---

# Instances that hit the time limit have no running time; they are counted
# separately and excluded from the statistics.
summ <- timing %>%
  group_by(density, d, n) %>% 
  summarise(
    instances  = n(),
    timeouts   = sum(timed_out),
    edges      = median(m),
    time_med   = median(time, na.rm = TRUE),
    time_q1    = quantile(time, 0.25, na.rm = TRUE, names = FALSE),
    time_q3    = quantile(time, 0.75, na.rm = TRUE, names = FALSE),
    time_max   = if (all(is.na(time))) NA else max(time, na.rm = TRUE),
    cond_med   = median(conduits, na.rm = TRUE),
    cond_q1    = quantile(conduits, 0.25, na.rm = TRUE, names = FALSE),
    cond_q3    = quantile(conduits, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  arrange(density, n)

print(summ, n = Inf)

summ$density <- ifelse(summ$density == "sparse", "Sparse", "Dense")

# --- Figure ---

p <- ggplot(summ, aes(x = n, y = time_med,
                      colour = density, fill = density,
                      linetype = density)) +
  geom_line() +
  geom_ribbon(aes(ymin = time_q1, ymax = time_q3),
              alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_point(size = 1) +
  # scale_x_log10(breaks = c(10, 20, 30, 50, 70, 100)) +
  scale_y_log10() +
  labs(x = "Number of vertices", y = "Time (s, log scale)",
       colour = "Graphs", fill = "Graphs", linetype = "Graphs") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.text.position = "bottom",
    legend.key.height = unit(4, "pt"),         # shorter key box around the line
    legend.key.width  = unit(14, "pt"),        # length of the line sample
    legend.text = element_text(size = 8, hjust = 0.5,
                               margin = margin(t = 1)),  # gap between line and name
    legend.key.spacing.x = unit(8, "pt"),      # gap between Dense and Sparse
    legend.margin = margin(2, 2, 2, 2),        # padding around the whole legend
    legend.title = element_blank(),
    axis.text = element_text(size = 7, face = "bold"),
    axis.title = element_text(size = 8, face = "bold")
  )

ggsave("timing.pdf", p, width = 2.75, height = 2)

p_conduits <- ggplot(summ, aes(x = n, y = cond_med,
                      colour = density, fill = density,
                      linetype = density)) +
  geom_line() +
  geom_ribbon(aes(ymin = cond_q1, ymax = cond_q3),
              alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_point(size = 1) +
  # scale_x_log10(breaks = c(10, 20, 30, 50, 70, 100)) +
  scale_y_log10() +
  labs(x = "Number of vertices", y = "Number of conduits (log scale)",
       colour = "Graphs", fill = "Graphs", linetype = "Graphs") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.text.position = "bottom",
    legend.key.height = unit(4, "pt"),         # shorter key box around the line
    legend.key.width  = unit(14, "pt"),        # length of the line sample
    legend.text = element_text(size = 8, hjust = 0.5,
                               margin = margin(t = 1)),  # gap between line and name
    legend.key.spacing.x = unit(4, "pt"),      # gap between Dense and Sparse
    legend.margin = margin(2, 2, 2, 2),        # padding around the whole legend
    legend.title = element_blank(),
    axis.text = element_text(size = 7, face = "bold"),
    axis.title = element_text(size = 8, face = "bold")
  )

ggsave("conduits.pdf", p_conduits, width = 2.75, height = 2)



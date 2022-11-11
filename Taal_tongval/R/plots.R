library(tidyverse)
library(ggrepel)
library(showtext) # for a fun font
font_add_google("Schoolbell", "bell")
showtext_auto()

# Functions -------------------------------------------------------------------
tsne_plot <- function(coords, seed = 7) {
  # Scatterplot with t-sne coordinates, colored by cluster (and negative silhouette
  # in gray). Instead of a legend the name of the medoid is indicated in its
  # approximate position (and the same color).
  ggplot(coords, aes(x = x, y = y, color = final_cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_label_repel(data = filter(coords, cw == medoid_label),
                     aes(label = cw),
                     seed = seed, size = 10) +
    theme_void() +
    scale_color_viridis_d(na.value = "gray", guide = "none") +
    coord_fixed()
}
add_point <- function(p, lemma_data, add_x = 1, add_y = 1, curvature = -0.2,
                      hjust = 1, vjust = -0.1, font_size = 10) {
  p + annotate("point", x = lemma_data$x, y = lemma_data$y, shape = 1, size = 5) +
    annotate("curve",
             x = lemma_data$x + add_x, xend = lemma_data$x, 
             y = lemma_data$y + add_y, yend = lemma_data$y,
             size = 1, arrow = arrow(length = unit(0.15, "inches")),
             curvature = curvature) +
    annotate("text",
             x = lemma_data$x + add_x, y = lemma_data$y + add_y,
             label = lemma_data$cw, hjust = hjust, vjust = vjust, family = "bell", size = font_size)
}

clusters_barplot <- function(tokens, column) {
  # Barplot with the number of tokens covered by each cluster, with different
  # shades of gray for each variant.
  ggplot(tokens, aes(x = fct_infreq({{ column }}), fill = Variant)) +
    geom_bar() +
    scale_fill_grey(start = 0.4, end = 0.6) +
    labs(x = "Medoid", y = "N of tokens") +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.9)) +
    coord_flip()
}
varimp_dotplot <- function(varimp_data) {
  varimp_data %>% 
    mutate(
      Predictor = str_remove(Predictor, ".bin"),
      Slot = str_extract(Predictor, "^[RT][a-z]+"),
      Slot = if_else (is.na(Slot), "Traditional", Slot)) %>% 
    ggplot(aes(x = Varimp, y = reorder(Predictor, Varimp),
               color = Slot)) +
    geom_point(size = 3) +
    scale_color_viridis_d(begin = 0.2, end = 0.8, option = "plasma", direction = -1) +
    labs(x = "Variable importance", y = "Predictor") +
    geom_vline(xintercept = 0,
               color = "gray30", linetype = 2) +
    theme_minimal(base_size = 25) +
    theme(legend.position = c(0.9, 0.2))
}

library(tidyverse)
library(semcloud)
library(Rtsne)
library(ggrepel)
library(ggparty)
library(showtext) # for a fun font
font_add_google("Schoolbell", "bell")
showtext_auto()

data_folder <- here::here("slides", "BKL_22", "data")
model_name <- list.dirs(data_folder, recursive = FALSE, full.names = FALSE)

# Read data ====================================================================
find_file <- function(pattern, subfolder = NA) {
  folder <- file.path(data_folder, model_name)
  if (!is.na(subfolder)) {
    folder <- file.path(folder, subfolder)
  }
  fname <- dir(folder, pattern = pattern, full.names = TRUE)
  if (!file.exists(fname)) {
    stop("File doesn't exist")
  }
  fname
}

## Filenames -------------------------------------------------------------------
fnames <- list(
  tokens = file.path(data_folder, "dat_us.tsv"),
  rds = find_file("*.rds"),
  wwmx_t = find_file("themes"),
  wwmx_r = find_file("recipients"),
  clusters_t = find_file("themes", "clusters"),
  clusters_r = find_file("recipients", "clusters")
)
## Load data -------------------------------------------------------------------
### Random forests ----
tree_data <- readRDS(fnames$rds)
varimp_bin <- partykit::varimp(tree_data$forest_bin)
tree_df <- as_tibble(ggparty(tree_data$tree)$data) %>%
  filter(is.na(splitvar)) %>%
  select(id, nodesize, starts_with("nodedata")) %>%
  unnest(starts_with("nodedata")) %>%
  select(-nodedata_fitted_nodes) %>% 
  distinct()

### Tokens with clustering ----
recipient_clusters <- read_tsv(fnames$clusters_r, show_col_types = FALSE)
theme_clusters <- read_tsv(fnames$clusters_t, show_col_types = FALSE) 
tokens <- read_tsv(fnames$tokens, show_col_types = FALSE) %>% 
  filter(Theme.inlist, Recipient.inlist) %>% 
  select(Token.ID, Context, Variant = Response.variable,
         ends_with("lemma"), ends_with("type")) %>% 
  mutate(Variant = fct_recode(Variant, ditransitive = "D", prepositional = "P")) %>%
  left_join(
    rename_with(recipient_clusters, ~ paste0("Recipient.", .x)),
    by = c("Recipient.lemma" = "Recipient.type_word")
    ) %>% 
  left_join(
    rename_with(theme_clusters, ~ paste0("Theme.", .x)),
    by = c("Theme.lemma" = "Theme.type_word")
  ) %>% 
  mutate(
    Recipient.medoid = if_else(Recipient.silhouette < 0, "other", Recipient.medoid_label),
    Theme.medoid = if_else(Theme.silhouette < 0, "other", Theme.medoid_label)
    ) %>% 
  left_join(
    select(tree_df, Recipient.medoid_label = nodedata_Recipient.cluster,
           Theme.medoid_label = nodedata_Theme.cluster, tree_cluster = id),
    by = c("Theme.medoid_label", "Recipient.medoid_label")) %>% 
  mutate(tree_cluster = as.factor(tree_cluster))

### Coordinates ----
get_coords <- function(filename) {
  set.seed(8541)
  mat <- focdistsFromCsv(filename) %>% 
    transformMats(TRUE)
  coords <- Rtsne(as.dist(mat), dims=2, perplexity=30,
          theta=0.0, check.duplicates = FALSE,
          max_iter = 1000, is_distance = TRUE)$Y
  tibble(
    cw = rownames(mat),
    x = coords[,1],
    y = coords[,2]
  )
}
recipient_coords <- get_coords(fnames$wwmx_r) %>% 
  left_join(recipient_clusters, by = c("cw" = "type_word")) %>% 
  mutate(final_cluster = if_else(silhouette < 0, NA_character_, medoid_label))
theme_coords <- get_coords(fnames$wwmx_t) %>% 
  left_join(theme_clusters, by = c("cw" = "type_word")) %>% 
  mutate(final_cluster = if_else(silhouette < 0, NA_character_, medoid_label))

# Plots =========================================================================

## Functions -------------------------------------------------------------------
tsne_plot <- function(coords, seed = 7) {
  # Scatterplot with t-sne coordinates, colored by cluster (and negative silhouette
  # in gray). Instead of a legend the name of the medoid is indicated in its
  # approximate position (and the same color).
  ggplot(coords, aes(x = x, y = y, color = final_cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_label_repel(data = filter(coords, cw == medoid_label),
                     aes(label = cw),
                     seed = seed) +
    theme_void() +
    scale_color_viridis_d(na.value = "gray", guide = "none") +
    coord_fixed()
}
add_point <- function(p, lemma_data, add_x = 1, add_y = 1, curvature = -0.2,
                      hjust = 1, vjust = -0.1, font_size = 6) {
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

clusters_barplot <- function(column) {
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
  enframe(varimp_data, name = "Predictor", value = "Varimp") %>% 
    mutate(
      Predictor = str_remove(Predictor, ".bin"),
      Slot = str_extract(Predictor, "^[RT][a-z]+")) %>% 
    ggplot(aes(x = Varimp, y = reorder(Predictor, Varimp),
               color = Slot)) +
    geom_point(size = 2) +
    scale_color_viridis_d(begin = 0.2, end = 0.8, option = "A") +
    labs(x = "Variable importance", y = "Predictor") +
    geom_vline(xintercept = abs(min(varimp_bin)),
               color = "gray30", linetype = 2) +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.2))
}

## Actual plotting -------------------------------------------------------------
tsne_plot(recipient_coords) %>% 
  add_point(recipient_coords %>% filter(cw == "government/nn"), -6, 0) %>% 
  add_point(recipient_coords %>% filter(cw == "mother/nn"), -7)

tsne_plot(theme_coords) %>% 
  add_point(theme_coords %>% filter(cw == "it/pp"), -7, -10, vjust = 1)

clusters_barplot(Recipient.medoid)
clusters_barplot(Theme.medoid)
clusters_barplot(tree_cluster)
varimp_dotplot(varimp_bin)
# add a freq plot of the distribution of LEMMAS

# Experiment ----
# check the main clusters in tree results
if (FALSE) {
  cluster_distrib <- tokens %>% 
    select(Token.ID, ends_with("label"), tree_cluster) %>% 
    unite(cluster_combo, ends_with("medoid_label"), sep = "-", remove = FALSE) %>% 
    add_count(tree_cluster, name = "N")
  get_top <- function(column_name) {
    cluster_distrib %>% add_count(tree_cluster, {{ column_name }}) %>% 
      select(tree_cluster, {{ column_name }}, N, n) %>% 
      distinct() %>% 
      group_by(tree_cluster) %>% 
      filter(n == max(n))
  }
  combo_top <- get_top(cluster_combo) %>% 
    rename(name = cluster_combo) %>% 
    mutate(what = "combo")
  r_top <- get_top(Recipient.medoid_label) %>% 
    rename(name = Recipient.medoid_label) %>% 
    mutate(what = "Recipient")
  t_top <- get_top(Theme.medoid_label) %>% 
    rename(name = Theme.medoid_label) %>% 
    mutate(what = "Theme")
  bind_rows(combo_top, r_top, t_top) %>% 
    group_by(tree_cluster) %>% 
    filter(n == max(n))
  
}

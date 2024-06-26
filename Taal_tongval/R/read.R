library(tidyverse)
library(semcloud)
library(Rtsne)
library(ggparty)
library(here)

data_folder <- here::here("Taal_tongval", "data")
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
  unite(cluster_combo, nodedata_Recipient.cluster, nodedata_Theme.cluster, sep = "-") %>% 
  select(cluster_combo, id) %>% 
  deframe()

### Tokens with clustering ----
recipient_clusters <- read_tsv(fnames$clusters_r, show_col_types = FALSE) 
theme_clusters <- read_tsv(fnames$clusters_t, show_col_types = FALSE)

tokens <- read_tsv(fnames$tokens, show_col_types = FALSE) %>% 
  filter(Theme.inlist, Recipient.inlist, Speaker.sex %in% c("F", "M")) %>% 
  select(Token.ID, Context, Variant = Response.variable,
         ends_with("lemma"), ends_with("type")) %>% 
  mutate(Variant = fct_recode(Variant, prepositional = "P", ditransitive = "D")) %>%
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
  unite(cluster_combo, Recipient.medoid_label, Theme.medoid_label, sep = "-", remove = FALSE) 

#%>% mutate(tree_cluster = as.factor(tree_df[cluster_combo]))

### Coordinates ----
get_coords <- function(filename) {
  set.seed(7543)
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


# Write data ===============================================================================

output_dir <- here::here("Taal_tongval", "data_processed")
output_files <- file.path(output_dir, paste0(c("dat_clusters", "recipient_coords", "theme_coords", "varimp_bin"),".tsv"))
walk2(
  list(tokens, recipient_coords, theme_coords, enframe(varimp_bin, name = "Predictor", value = "Varimp")),
  output_files,
  write_tsv
)
# # Experiment ----
# # check the main clusters in tree results
# if (FALSE) {
#   cluster_distrib <- tokens %>% 
#     select(Token.ID, ends_with("label"), tree_cluster) %>% 
#     unite(cluster_combo, ends_with("medoid_label"), sep = "-", remove = FALSE) %>% 
#     add_count(tree_cluster, name = "N")
#   get_top <- function(column_name) {
#     cluster_distrib %>% add_count(tree_cluster, {{ column_name }}) %>% 
#       select(tree_cluster, {{ column_name }}, N, n) %>% 
#       distinct() %>% 
#       group_by(tree_cluster) %>% 
#       filter(n == max(n))
#   }
#   combo_top <- get_top(cluster_combo) %>% 
#     rename(name = cluster_combo) %>% 
#     mutate(what = "combo")
#   r_top <- get_top(Recipient.medoid_label) %>% 
#     rename(name = Recipient.medoid_label) %>% 
#     mutate(what = "Recipient")
#   t_top <- get_top(Theme.medoid_label) %>% 
#     rename(name = Theme.medoid_label) %>% 
#     mutate(what = "Theme")
#   bind_rows(combo_top, r_top, t_top) %>% 
#     group_by(tree_cluster) %>% 
#     filter(n == max(n))
#   
# }

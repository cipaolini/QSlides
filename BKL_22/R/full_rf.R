library(tidyverse)
library(partykit)

traditional_predictors <- read_tsv(here::here("BKL_22", "data", "dat_us.tsv"), show_col_types = FALSE) %>% 
  filter(Speaker.sex %in% c("F", "M")) %>% 
  mutate(
    Recipient.type.bin = if_else(Recipient.type == 'N', 'N', 'P'),
    Theme.type.bin = if_else(Theme.type == 'N', 'N', 'P'),
    Recipient.definiteness.bin = if_else(Recipient.definiteness == 'Indefinite', 'Indefinite', 'Definite'),
    Theme.definiteness.bin = if_else(Theme.definiteness == 'Indefinite', 'Indefinite', 'Definite'),
    Recipient.animacy.bin = if_else(Recipient.animacy == 'A', 'A', 'I'),
    Theme.animacy.bin = if_else(Theme.animacy == 'A', 'A', 'I'),
    Length.difference = log(Recipient.length) - log(Theme.length),
    Theme.lemma.pruned = fct_lump_min(Theme.lemma, 2, other_level = "other"),
    Recipient.lemma.pruned = fct_lump_min(Recipient.lemma, 2, other_level = "other")
  ) %>% 
  select(
    Token.ID, Semantics,
    Recipient.type.bin, Recipient.definiteness.bin, Recipient.animacy.bin,
    Theme.definiteness.bin, Length.difference, ends_with("pruned")
  )

# A few lines from the tokens are lost with the merge
binarize <- function(df) {
  df %>% mutate(cluster = str_replace(cluster, "([^/]+)/.*", "\\1") %>% str_replace("-", "_")) %>% 
    pivot_wider(names_from = cluster, values_from = cluster, values_fill = "no_predictor") %>% 
    rename_with(~paste0(.x, ".bin"), -Token.ID)
}

significant_dp <- read_tsv(here::here("BKL_22", "data_processed", "varimp_bin.tsv"),
         show_col_types = FALSE) %>% 
  filter(Varimp > 0) %>% 
  pull(Predictor)
clusters <- read_tsv(here::here("BKL_22", "data_processed", "dat_clusters.tsv"),
         show_col_types = FALSE) %>% 
  filter(!is.na(Recipient.medoid), !is.na(Theme.medoid)) %>% 
  select(Token.ID, Variant, ends_with("medoid"))
recipient_bin <- select(clusters, Token.ID, cluster = Recipient.medoid) %>% 
  binarize()
theme_bin <- select(clusters, Token.ID, cluster = Theme.medoid) %>% 
  binarize()
dataset <- clusters %>% 
  select(-ends_with("label")) %>% 
  left_join(
    rename_with(recipient_bin, ~ paste0("Recipient.", .x), -Token.ID),
    by = "Token.ID"
  ) %>% 
  left_join(
    rename_with(theme_bin, ~ paste0("Theme.", .x), -Token.ID),
    by = "Token.ID"
  ) %>% 
  select(Token.ID, Variant, all_of(significant_dp)) %>% 
  inner_join(traditional_predictors, by = "Token.ID") %>% 
  select(-Token.ID, -ends_with("pruned")) %>% 
  mutate(across(where(is.character), as.factor))

forest <- cforest(Variant ~ ., data = dataset)
varimp_full <- varimp(forest, conditional = FALSE)
varimp_full %>% enframe(name = "Predictor", value = "Importance") %>% 
  write_tsv(here::here("BKL_22", "data_processed", "varimp_full.tsv"))
#varimp_dotplot(varimp_full)

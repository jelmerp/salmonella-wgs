read_sample <- function(mut_file, meta_df) {
  read_tsv(mut_file, show_col_types = FALSE) %>%
    rename(sample_id = Name,
           gene = `Gene symbol`,
           description = `Sequence name`) %>%
    clean_names() %>%
    select(sample_id, gene, description, class, subclass) %>% 
    mutate(mutation = sub("(.*)_([^_]+$)", "\\2", gene),
           gene = sub("(.*)_([^_]+$)", "\\1", gene)) %>%
    left_join(meta_df, by = "sample_id") %>%
    relocate(serotype, .after = "sample_id")
}

## Function to get a per-gene, per-subclass count of mutations
get_counts <- function(mut_df) {
  mut_df %>%
    filter(!is.na(class)) %>%
    group_by(gene, subclass) %>%
    add_count() %>%
    distinct(gene, subclass, .keep_all = TRUE) %>%
    ungroup()
}
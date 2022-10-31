## Packages
library(tidyverse)
library(readxl)
library(janitor)

## Scripts
source("scripts/resistance_fun.R")

## Define input files
mut_dir <- "results/amrfinderplus"
mut_files <- list.files(mut_dir, pattern = "mutation-report.txt", full.names = TRUE)
meta_file <- "results/summary of the serotype_WGS_salmonella.xlsx"

## Read and prep metadata
meta <- read_xlsx(meta_file) %>%
  select(sample_id = `Sample number`, serotype = Serotype)

## Read Amrfinder plus results
mut <- map_dfr(mut_files, read_sample, meta_df = meta) 

## Get per-gene, per-class mutation count for each sample
mut %>%
  select(-sample_id, -serotype) %>% 
  filter(!is.na(class)) %>%
  group_by(gene, subclass) %>%
  add_count() %>%
  distinct(gene, subclass, .keep_all = TRUE) %>%
  ungroup()

## Per-mutation count
per_mut <- mut %>%
  select(-sample_id, -serotype) %>% 
  group_by(gene, mutation) %>%
  add_count() %>%
  distinct(gene, subclass, .keep_all = TRUE) %>%
  ungroup()

per_mut %>% filter(n != 51)


## Packages
if (!"pacman" %in% installed.packages())
  install.packages("pacman", repos='http://cran.us.r-project.org')
if (!"BiocManager" %in% installed.packages())
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
if (!"ggtree" %in% installed.packages()) BiocManager::install("ggtree")
packages <- c("tidyverse", "here", "ape", "ggtree", "argparse", "readxl")
pacman::p_load(char = packages, install = TRUE, repos = 'http://cran.us.r-project.org')

## Define input files
meta_file <- "results/summary of the serotype_WGS_salmonella.xlsx"
meta_ref_file <- "results/patric_refs/lookup.txt"
tree_file <- "results/gubbins/all/all.final_tree.tre"

tree <- read.tree(tree_file)
tree$tip.label <- sub("AG21-0", "", tree$tip.label)

## Read and prep metadata
meta_ref <- read_tsv(meta_ref_file, show_col_types = FALSE) %>% 
  select(sampleID = accession, serotype) %>%
  mutate(serotype = str_to_title(serotype),
         serotype = sub("Branderup", "Brenderup", serotype))

meta <- read_xlsx(meta_file) %>%
  mutate(sampleID = sub("AG21-0", "", `Sample number`)) %>%
  filter(sampleID != "Average") %>% 
  select(sampleID, serotype = Serotype) %>%
  bind_rows(meta_ref)

## Read the tree file and prep the tree
tree <- read.tree(tree_file)
tree$tip.label <- sub("AG21-0", "", tree$tip.label)
nseqs <- length(tree$tip.label)
#tree$node.label <- round(as.numeric(tree$node.label) * 100)

## Define colors
#randomcoloR::distinctColorPalette(20)
mycols <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
            "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
            "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736")
sero_cols <- mycols[1:length(unique(meta$serotype))]



p <- ggtree(tree, layout = "rectangular") %<+% meta +    # or layout = "circular"
  geom_tiplab(aes(color = serotype), size = 2.5) +
  geom_treescale(x = 0, y = nseqs - 1, color = "grey50") +
  geom_rootedge(rootedge = 0.005) +
  scale_color_manual(values = sero_cols) +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"),
        legend.position = "left") +
  coord_cartesian(clip = "off")
p


tree_ksnp_file <- "results/ksnp3/with_refs/tree.core.tre"
tree_ksnp <- read.tree(tree_file)

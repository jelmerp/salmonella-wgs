# SET UP ------------------------------------------------------------------------
## Docs
# - https://github.com/xavierdidelot/BactDating
# - https://xavierdidelot.github.io/BactDating/articles/exampleRec.html

## Packages
#devtools::install_github("xavierdidelot/BactDating")
library(BactDating)
library(ape)

## Input files
#gubbins_prefix <- "results/gubbins/all/all"
#dates_file <- "results/dates/all_dates.tsv"
gubbins_prefix <- "results/gubbins/cerro/cerro"
dates_file <- "results/dates/cerro_dates.tsv"


# PREP THE TREE ----------------------------------------------------------------
## Load the Gubbins tree
tree <- loadGubbins(prefix = gubbins_prefix)

## Load the dates
dates_df <- read.table(dates_file, header = TRUE)
dates_df <- dates_df[match(tree$tip.label, dates_df$name), ]
dates <- lubridate::decimal_date(as.Date(dates_df$collection_date))
names(dates) <- dates_df$name

## Root the tree
tree_rooted <- initRoot(tree, dates, useRec = TRUE)
plot(tree_rooted)


# TESTING TEMPORAL SIGNAL ------------------------------------------------------
## Root-to-tip regression analysis
res <- roottotip(tree_rooted, dates)
#> R2=0.09, p=1.94e-02 => Too weak [ALL PATHOVARS]

## Testing temporal signal - https://xavierdidelot.github.io/BactDating/articles/exampleTest.html
dates_same <- rep(mean(dates), length(dates))
tree_dated_mock <- bactdate(tree_rooted, dates_same, useRec = TRUE)
modelcompare(tree_dated, tree_dated_mock)
#> The first model has DIC=1152.45 and the second model has DIC=961.57. [ALL PATHOVARS]
#> Model 2 is definitely better. [ALL PATHOVARS]

## Clustered permutation test - https://xavierdidelot.github.io/BactDating/articles/clustered.html
#res <- clusteredTest(tree_rooted, dates) # Doesn't work


# DATING -----------------------------------------------------------------------
## Date the tree without a specified rate
tree_dated <- bactdate(tree_rooted, dates, useRec = TRUE,
                       nbIts = 10000, showProgress = TRUE)
plot(tree_dated, "treeCI")

## Use the 2.4e-7 rate
tree_dated_rg <- bactdate(tree_rooted, dates, useRec = TRUE, nbIts = 10000,
                          initMu = 2.4e-7 * 50e6, updateMu = TRUE,
                          model = "relaxedgamma")
plot(tree_dated_rg, "treeCI")

tree_dated_arc <- bactdate(tree_rooted, dates, useRec = TRUE, nbIts = 10000,
                          initMu = 2.4e-7 * 50e6, updateMu = TRUE,
                          model = "arc")
plot(tree_dated_arc, "treeCI")

tree_dated_carc <- bactdate(tree_rooted, dates, useRec = TRUE, nbIts = 10000,
                           initMu = 2.4e-7 * 50e6, updateMu = TRUE,
                           model = "carc")
plot(tree_dated_carc, "treeCI")


## Check MCMC traces
plot(tree_dated_rooted, "trace")

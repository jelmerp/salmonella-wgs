## Load packages
library(corrplot)

## Define input and output files
infile <- "results/resistance/mock.csv"
plotfile <- "results/resistance/corplot.png"

## Read data
mock_data <- read.csv(infile, row.names = 1)
## Get rid of periods (.) in column names
colnames(mock_data) <- gsub("\\.", "", colnames(mock_data))
## Remove columns with only 0s or only 1s
## These are non-informative for this analysis, and are the reasons for the error with cor.test
mock_data <- mock_data[, -which(colSums(mock_data) %in% c(0, nrow(mock_data)))]

## Get p-value and correlation matrices
p_mat <- cor.mtest(mock_data)$p
cor_matrix <- cor(mock_data)

## Correlation plot
png(plotfile)
corrplot(cor_matrix,
         p.mat = p_mat,        # Use the p-value matrix
         insig = "blank",      # Make cells with non-significant correlations blank
         diag = FALSE,         # Don't show the diagonal ("self-correlations")
         type = "upper",    
         tl.col = "black",
         tl.srt = 45)
dev.off()

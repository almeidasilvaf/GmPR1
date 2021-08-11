# Create a data frame of GmPR-1 gene expression on atlas and recently downloaded stress-related samples
# Jan 2021

library(tidyverse)

#----Create character vector of GmPR-1 genes----
pr1 <- read.csv(here("products", "tables", "table1.txt"))[,1]

#----Create character vector of abiotic and biotic stress-related samples----
load(here("data", "stress_samples_for_DGE_analysis.Rdata"))
abiotic_samples <- Reduce(rbind, abiotic_samplelist)$BioSample
biotic_samples <- Reduce(rbind, biotic_samplelist)$BioSample
stress_samples <- c(abiotic_samples, biotic_samples)

#----Create data frame of raw counts for atlas samples----
atlas <- load(here("data", "atlas.rda"))
atlas <- atlas[pr1, ]

#----Load count matrix with new samples (on server)----
counts <- load(here("data", "new_samples.rda"))
counts <- as.data.frame(counts)
counts <- counts[counts$gene_id %in% pr1, ]
counts <- counts[, colnames(counts) %in% c("gene_id", stress_samples)]

#----Create data frame of global expression: atlas + new samples----
globalexp <- merge(atlas, counts, by.x="row.names", by.y="gene_id")
rownames(globalexp) <- globalexp$Row.names
globalexp$Row.names <- NULL

#----Create data frame of stress-related expression----
stressexp <- globalexp[, colnames(globalexp) %in% stress_samples]

#----Save objects----
save(globalexp, stressexp, abiotic_samples, biotic_samples, stress_samples,
     file = here("data", "expression.RData"),
     compress = TRUE)






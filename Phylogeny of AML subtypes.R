library(Biobase)
library(maphylogeny)
#Riester et al., package
library(dplyr)
library(genefilter)

df <- read.csv("Normalized expression data.csv", fileEncoding="UTF-8-BOM")
df <- tibble::column_to_rownames(df, var = "X")
## import TCGA patient expression data.frame

var <- read.csv("Table of samples with AML subtype", fileEncoding="UTF-8-BOM")
var <- tibble::column_to_rownames(var, var = "X")
var <- AnnotatedDataFrame(var)
##REQUIRED for grouping TCGA patients by AML subtype 


mx <- as.matrix(df)
## convert data.frame to matrix (format required for ExpressionSet)


es <- ExpressionSet(mx,
                    phenoData = var,
                    annotation = "org.Hs.eg.db")
## generate expression set (format required for filtering and bootstrapping)


fes <- varFilter(es,
                 var.func = IQR, 
                 var.cutoff = 0.5, 
                 filterByQuantile = TRUE)
## Remove equally expressed genes from data with low variance filter


dists <- maphylo_bootstrap(fes,
                           group = var2,
                           r = c(1),
                           bootstrap = 1000,
                           dm = "pearson")
## bootstraping, repetitions beyond 1,000 did not change alter consensus tree

trees = maphylo_reconstruct(dists)
## use Neighbor-Joining to reconstruct trees

plot(trees[[1]])

maphylo_consensus_phylip(trees, outfile = "")
## export phylogeny as a .tree file 





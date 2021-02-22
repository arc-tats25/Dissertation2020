library(Biobase)
library(maphylogeny)
library(dplyr)
library(genefilter)

df <- read.csv("Normalized expression data.csv")
df <- tibble::column_to_rownames(df, var = "X")
## import assay data as a data.frame

var <- read.csv("Table of samples with AML subtype", fileEncoding="UTF-8-BOM")
var <- tibble::column_to_rownames(var, var = "X")
var <- AnnotatedDataFrame(var)
var2 <- as.factor(var$AML_subtype)
## correlate sample from assay data to AML subtype and transform as a factor for bootstrapping

mx <- as.matrix(df,
                nrow = 173,
                ncol = 19448)
## convert data.frame to matrix


es <- ExpressionSet(mx,
                    phenoData = var,
                    annotation = "org.Hs.eg.db")
## generate expression set


fes <- varFilter(es,
                 var.func = IQR, 
                 var.cutoff = 0.5, 
                 filterByQuantile = TRUE)
## filter uninformative genes with low variance filter


fes2 <- as(fes, "ExpressionSet")
fes2 <- as(fes2, "data.frame")
##for ensuring it worked

dists <- maphylo_bootstrap(fes,
                           group = var2,
                           r = c(1),
                           bootstrap = 1000,
                           dm = "pearson")
## bootstrap

trees = maphylo_reconstruct(dists)
## use Neighbor-Joining to reconstruct trees

plot(trees[[1]])

maphylo_consensus_phylip(trees, outfile = "test2")
## export tree



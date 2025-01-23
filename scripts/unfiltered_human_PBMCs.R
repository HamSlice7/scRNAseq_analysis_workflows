library(DropletTestFiles)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
##Loading in the data
raw.path <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.tar.gz")
out.path <- file.path(tempdir(), "pbmc4k")
untar(raw.path, exdir = out.path)


fname <- file.path(out.path, "raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

#Changing row names
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

#Getting gene locations
location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce.pbmc)$ID, column = "SEQNAME", keytype = "GENEID")


##Quality control - removing droplets with ambient RNA
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

unfiltered <- sce.pbmc

#mitochondria - removing cells with too much mtRNA
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
sce.pbmc <- sce.pbmc[,!high.mito]

summary(high.mito)

##Plotting the cells that are removed or not with mitochondrial QC
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- high.mito

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Totoal count"),
  plotColData(unfiltered, y="detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent", colour_by = "discard") +
    ggtitle("Mito percent"), ncol = 2
)

#Plotting proportion of mitochondrial reads in each cell of the PMBC dataset compared to its total count
plotColData(unfiltered, x="sum", y="subsets_Mito_percent", colour_by = "discard") + scale_x_log10()

##Normalization
set.seed(1000)
#Pooling cells before normalization
clusters <- quickCluster(sce.pbmc)
#Scaling normalization by deconvolving size factors from cell pools
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

summary(sizeFactors(sce.pbmc))
summary(librarySizeFactors(sce.pbmc))

#Relationship between the library size factors and deconvolution size factors in the PBMC
plot(librarySizeFactors(sce.pbmc), sizeFactors(sce.pbmc), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")

##Variance modelling
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)

#Per-gene variance as a function of the mean for the log-expression values in the PBMC dataset
plot(dec.pbmc$mean, dec.pbmc$total, pch=16, cex=0.5,
     xlab = "Mean of log-expression", ylab = "Variance of log-expression")
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

##Dimensionality reduction
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row = top.pbmc, technical = dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")

#Verifying that there is a reasonable number of PCs retained
ncol(reducedDim(sce.pbmc,"PCA"))


##Clustering
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)

table(colLabels(sce.pbmc))

plotTSNE(sce.pbmc, colour_by = "label")


##Interpretation
markers <- findMarkers(sce.pbmc, pval.type="some", direction = "up")
marker.set <- markers[["2"]]
as.data.frame(marker.set[1:30,1:3])

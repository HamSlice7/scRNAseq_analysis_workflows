library(DropletTestFiles)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)

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

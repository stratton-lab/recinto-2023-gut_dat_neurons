suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))

SEURAT_VERBOSE <- FALSE

seurat_pipeline <- function(so) {
  so <- FindVariableFeatures(
    so,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = SEURAT_VERBOSE
  )
  so <- ScaleData(so, features = rownames(so), verbose = SEURAT_VERBOSE)
  so <- RunPCA(so, verbose = SEURAT_VERBOSE)
  so <- FindNeighbors(so, verbose = SEURAT_VERBOSE)
  so <- FindClusters(so, resolution = 1, verbose = SEURAT_VERBOSE)
  so <- RunUMAP(so, dims = 1:10, verbose = SEURAT_VERBOSE)

  so
}

tdt_pos <- CreateSeuratObject(
  Read10X("data/GSE237170"),
  min.cells = 3,
  min.features = 200
)

tdt_pos[["orig.ident"]] <- "tdTom+"

tdt_pos <- PercentageFeatureSet(
  tdt_pos,
  pattern = "^mt-",
  col.name = "percent.mt"
)
tdt_pos <- subset(tdt_pos, subset = nFeature_RNA > 200 & percent.mt < 5)
tdt_pos <- NormalizeData(
  tdt_pos,
  scale.factor = 1000,
  verbose = SEURAT_VERBOSE
)
tdt_pos <- seurat_pipeline(tdt_pos)

tdt_pos <- doubletFinder_v3(
  tdt_pos,
  pN = 0.25,
  pK = 0.09,
  nExp = round(ncol(tdt_pos) * 0.05),
  PCs = 1:10
)
tdt_pos <- subset(tdt_pos, DF.classifications_0.25_0.09_329 == "Singlet")
tdt_pos <- PercentageFeatureSet(
  tdt_pos,
  pattern = "td",
  col.name = "percent.td"
)

# NOTE: the following markers were used to annotate clusters
neuron <- c("Elavl4", "Tubb3")
neuroblast <- "Ascl1"
progenitor <- "Sox10"
glia <- c("S100b", "Plp1")
immune <- "Ptprc"
endothelial <- c("Pdgfra", "Pdgfrb")
smc <- c("Cnn1", "Mylk", "Tpm2", "Tpm1", "Des", "Myh11")
icc <- "Cd117"
mm <- c("Cx3cr1", "Fcgr1a", "Adgre1")

# NOTE: when reproducing results, ensure that cluster numbers have not changed
neurons <- subset(tdt_pos, idents = c(11, 14))

neurons <- seurat_pipeline(neurons)

# NOTE: the following markers were used to determine neuronal subpopulations
cholinergic <- c("Chat", "SlC28a3")
nitrergic <- c("Nos1")
gabaergic <- c("Gad2")
glutamatergic <- c("Slc17a6")
catecholaminergic <- c("Th", "Dbh", "Ddc", "Slc18a2")
other <- c("Gal", "Vip", "Calcb", "Calb1", "Nefm")

markers <- FindAllMarkers(
  neurons,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = SEURAT_VERBOSE
)

saveRDS(neurons, "out/neurons.rds")
write.csv(markers, "out/markers.csv")

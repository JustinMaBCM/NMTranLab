#identifying 80% confidence from cell ranger based multiplets
CR0.8_Assign <- function(DF, CMOnum = 6){
  DFtemp <- DF[,c(2:(CMOnum+1),(CMOnum+3):(CMOnum+6), CMOnum+2)]
  DFtemp$CR0.8 <- 0
  for (i in 1:nrow(DFtemp)){
    if (DFtemp$Assignment_Probability[i] >= 0.8) {DFtemp$CR0.8[i] <- colnames(DFtemp)[which(match(DFtemp[i,1:(CMOnum+2)],DFtemp[i,(CMOnum+4)]) == 1)]
    } else {DFtemp$CR0.8[i] <- "Unassigned"}
  }
  return(DFtemp$CR0.8)
}

#default non-SCT Processing
QuickProcess <- function(SO, resolution = 0.5, scale = TRUE,PCA = TRUE, harmony = FALSE){
  SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^mt-")
  SO <- NormalizeData(object = SO, normalization.method = "LogNormalize", scale.factor = 10000)
  SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = 5000); gc()
  
  if (scale == TRUE) {SO <- ScaleData(SO, features = rownames(SO)); gc()}
  
  if (PCA == TRUE){SO <- RunPCA(SO, verbose = FALSE, npcs = 100); gc()}
  
  if (harmony == FALSE){
    
    SO <- RunUMAP(SO, reduction = "pca", dims = 1:100, verbose = FALSE)
    SO <- FindNeighbors(SO, reduction = "pca", dims = 1:100, verbose = TRUE)
    SO <- FindClusters(SO, verbose = TRUE, resolution = resolution)
  }
  
  if (harmony == TRUE){
    options(repr.plot.height = 2.5, repr.plot.width = 6)
    SO <- SO %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, assay.use = "RNA")
    SO <- RunUMAP(SO, reduction = "harmony", dims = 1:100, verbose = FALSE)
    SO <- FindNeighbors(SO, reduction = "harmony", dims = 1:100, verbose = TRUE)
    SO <- FindClusters(SO, verbose = TRUE, resolution = resolution)
  }
  return(SO)
}

#identification of adult retinal cells mouse
RetinaID <- function(SO, group.by = "seurat_clusters", width = 4, height = 8){
  require(ggplot2)
  require(ragg)
  Rods <- c("Rho","Pdc","Nrl","Sag","Gnat1","Gngt1","Nr2e3","Gnb1")
  Cones <- c("Arr3","Rcvrn","Gnat2","Gngt2","Opn1sw","Opn1mw","Guca1a")
  BP <- c("Vsx2","Otx2","Grm6","Prkca","Trpm1","Grik1","Vsx1","Cabp5")
  RGCs <- c("Rbpms","Pou4f1","Pou4f2","Pou4f3","Thy1", "Slc17a6","Nefl","Nefm","Sncg","Opn4")
  AC <- c("Tfap2a","Tfap2b","Tfap2c","Gad1","Gad2","Slc6a9","C1ql1","C1ql2","Lgr5")
  HZ <- c("Lhx1","Onecut1","Onecut2","Calb1")
  MuGl <- c("Rlbp1","GluI","Apoe","Crabp1","Clu","Slc1a3","Dkk3","Crym")
  MiGl <- c("C1qa","C1qb","C1qc","Hexb","Ctss","P2ry12","Tmem119","B2m")
  Endo <- c("Cldn5","Igfbp7","Col4a1","Vegfa","Vegfb","Vegfc","Pecam1","Pdgfrb","Cspg4","Anpep","Acta2")
  Endo2 <- c("Kcnj8","Atp1a2","Itih5","Cspg4","Ramp2","Pecam1","Cldn5","Acta2","Myh11","Crip1","Parm1","Timp3")
  GeneList <- list(Rods,Cones,BP,RGCs,AC,HZ,MuGl,MiGl,Endo,Endo2)
  names(GeneList) <- c("Rods","Cones","BP","RGCs","AC","HZ","MuGl","MiGl","Endo","Endo2")
  for(i in 1:length(GeneList)){
    plot <- (DotPlot(SO, features = c(GeneList[i]), assay = "RNA",group.by = group.by, cols = c("grey", "black"), dot.scale = 6, dot.min = 0, cluster.idents = FALSE, scale.by = 'size') + RotatedAxis() + theme())
    ggsave(plot, filename = sprintf("Dotplot_%s.png", names(GeneList)[i]),width = width, height = height, scaling = 0.9, bg = "white")
  }
}

#Processing HTODemux for CMO data
HTODemux_f <- function(SO, quantile = 0.99){
  SO <- NormalizeData(SO, normalization.method = "CLR")
  SO <- HTODemux(SO, assay = "RNA", positive.quantile = quantile)
  SO$hash.ID <- factor(SO$hash.ID, levels = sort(levels(SO$hash.ID)))
  return(SO)
}

CMO_Process <- function(SO, CMO){
  SO <- subset(SO, features = CMO)
  SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000)
  SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = length(CMO))
  SO <- ScaleData(SO, features = rownames(SO))
  SO <- RunPCA(SO, npcs = length(CMO), approx = FALSE)
  SO <- RunUMAP(SO, dims = 1:length(CMO), verbose = FALSE)
  SO <- RunTSNE(SO, dims = 1:length(CMO), verbose = FALSE, check_duplicates = FALSE)
  SO <- FindNeighbors(SO, dims = 1:length(CMO), verbose = TRUE)
  SO <- FindClusters(SO, verbose = TRUE, dims = 1:length(CMO), resolution = 1)
  #HTODemux
  SO <- HTODemux_f(SO, quantile = 0.99)
  return(SO)
}

#DietExport
DietExport <- function(SO){
  export = CreateSeuratObject(counts = LayerData(SO, assay = "RNA", layer = "counts"), project = "Export", meta.data = SO@meta.data)
  export@reductions$pca <- SO@reductions$pca; export@reductions$umap <- SO@reductions$umap
  SO <- export; rm(export)
  return(SO)
}

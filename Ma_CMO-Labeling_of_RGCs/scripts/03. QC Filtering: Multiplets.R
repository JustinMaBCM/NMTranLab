#Here we perform filters for multiplets and MLCs

#Filter 3: we remove c0 which is primarily RGC/non-RGC multiplets
ExpAll <- subset(ExpAll, idents = 0, invert = TRUE)
ExpAll <- QuickProcess(ExpAll)

#Identify MLCs from RGC subclusters to remove them. RGC subclusters were defined using known markers and cluster proximity
#to summary for each subgroup isolated, we validate the clusters with the appropriate markers and based on the dotplots, evaluate which clusters are MLCs (labeled XX_M)
#these are then mapped back to the master cluster as $newidents

#i3
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_i3 <- subset(ExpAll, idents = c(36,20,29)); ExpAll_i3 <- QuickProcess(ExpAll_i3)
DimPlot(ExpAll_i3, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_i3, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_i3, features = unique(c("Irx3","Prkcg","Ceacam10","Pcdh20","Pcdh11x","Pcdh20","Prkg2")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
DotPlot(ExpAll_i3, features = unique(c(i3genes)),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_i3$newident <- mapvalues(ExpAll_i3$seurat_clusters, from = c(0,6,1,3,2,4,5), to = c("27_Novel_1","27_Novel_2","i3_M","37_Novel","18_Novel_b","18_Novel_a1","18_Novel_a2"))
DimPlot(ExpAll_i3, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_i3@active.ident <- factor(ExpAll_i3$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_i3, Drop = TRUE)
ExpAll$newidents <- ExpAll$SubclusterID

#ooDS
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_ooDS <- subset(ExpAll, idents = c(4,12,14)); ExpAll_ooDS <- QuickProcess(ExpAll_ooDS)
ExpAll_ooDS <- FindClusters(ExpAll_ooDS, resolution = 10)
DimPlot(ExpAll_ooDS, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_ooDS, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_ooDS, features = unique(c("Bnc2","Gpr88","Pde11a","Gpr101","Tafa4","Ptprt","Calb1","Col25a1","Serpinb8","Efnb2","Fam129a")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_ooDS$newident <- mapvalues(ExpAll_ooDS$seurat_clusters, from = c(0,1,8,2,3,4,6,7,5), to = c("ooDS_M","24_Novel","24_Novel","16_ooDS_D","16_ooDS_V",rep("10_Novel",3),"16_ooDS_DVmix"))
DimPlot(ExpAll_ooDS, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_ooDS@active.ident <- factor(ExpAll_ooDS$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_ooDS, Drop = TRUE)
ExpAll$newidents <- ExpAll$SubclusterID

#MK
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_MK <- subset(ExpAll, idents = c(24,28,21)); ExpAll_MK <- QuickProcess(ExpAll_MK)
DimPlot(ExpAll_MK, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_MK, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_MK, features = unique(c("Mafb","Kcnd2","Prokr1","Vgf","Fzd6","Tpbg","Spp1","Chrm1","Kit","Fes")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_MK$newident <- mapvalues(ExpAll_MK$seurat_clusters, from = c(0,6,1,4,2,5,3,7), to = c(rep("34_Novel",2),rep("30_Novel",2),rep("23_W3D2",2),"MK_M","45_AlphaOFFT"))
DimPlot(ExpAll_MK, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_MK@active.ident <- factor(ExpAll_MK$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_MK, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#ipAE
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_ipAE <- subset(ExpAll, idents = c(9,26,2,6,39,34)); ExpAll_ipAE <- QuickProcess(ExpAll_ipAE)
ExpAll_ipAE <- FindClusters(ExpAll_ipAE, resolution = 10)
DimPlot(ExpAll_ipAE, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_ipAE, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_ipAE, features = unique(c("Eomes","Cdh6","Gm11100","Irx1","Opn4","Adra2a","Pnoc","Nmb","Gpc5","Cdhr1","Spp1","Il1rapl2","Kit","Fes","Tpbg","Chrm1")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_ipAE$newident <- mapvalues(ExpAll_ipAE$seurat_clusters, from = c(0,1,7,2,8,3,4,5,6,9,10,11), to = c("ipAE_M",rep("8_Novel",2),rep("7_Novel",2),"22_M5","42_AlphaONS","31_M2","33_M1","40_M1dup","43_AlphaONS","41_AlphaONT"))
DimPlot(ExpAll_ipAE, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_ipAE@active.ident <- factor(ExpAll_ipAE$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_ipAE, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#TFp
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_TFp <- subset(ExpAll, idents = c(10,13,38,17)); ExpAll_TFp <- QuickProcess(ExpAll_TFp)
ExpAll_TFp <- FindClusters(ExpAll_TFp, resolution = 2)
DimPlot(ExpAll_TFp, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_TFp, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_TFp, features = unique(c("Foxp2","Sema5a","Tbr1","Slc7a11","Satb2","Pcdh20","Calca","Zeb2","Slc35f3","Irx4")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_TFp$newident <- mapvalues(ExpAll_TFp$seurat_clusters, from = c(0,7,3,13,1,2,4,5,6,8,11,10,9,14,12), to = c(rep("9_Tbr1_Novel",4),rep("5_J-RGC",4),rep("17_Tbr1_S1",4),rep("28_FmidiOFF",2),"TFp_M"))
DimPlot(ExpAll_TFp, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_TFp@active.ident <- factor(ExpAll_TFp$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_TFp, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#T5_F
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_T5 <- subset(ExpAll, idents = c(1,31,0,8,23)); ExpAll_T5 <- QuickProcess(ExpAll_T5)
ExpAll_T5 <- FindClusters(ExpAll_T5, resolution = 2)
DimPlot(ExpAll_T5, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_T5, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_T5, features = unique(c("Trarg1","Amigo2","Cmtm8","Lypd1","Ldb2","Foxp2","Anxa3","Gria1")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
DotPlot(ExpAll_T5, features = unique(c(rgcgenes)),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_T5$newident <- mapvalues(ExpAll_T5$seurat_clusters, from = c(0,6,10,9,16,1,14,7,11,5,2,3,12,13,8,4,15), to = c(rep("3_FminiOFF",4),"3_FminiOFF_e",rep("2_W3D1.2",5),rep("1_W3D1.1",6),"T5F_M"))
DimPlot(ExpAll_T5, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_T5@active.ident <- factor(ExpAll_T5$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_T5, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#F
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_F <- subset(ExpAll, idents = c(3,32)); ExpAll_F <- QuickProcess(ExpAll_F)
ExpAll_F <- FindClusters(ExpAll_F, resolution = 5)
DimPlot(ExpAll_F, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_F, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_F, features = unique(c("Foxp2","Coch","Irx4","Gria1","Anxa3","Sema5a")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_F$newident <- mapvalues(ExpAll_F$seurat_clusters, from = c(0,18,11,15,6,16,7,13,9,14,12,17,3,4,1,2,5,8,19,10), to = c(rep("4_FminiOFF",15),rep("32_F_Novel",3),"F_M","F_M"))
DimPlot(ExpAll_F, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_F@active.ident <- factor(ExpAll_F$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_F, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#NC1
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_NC1 <- subset(ExpAll, idents = c(11,5,27,35,15)); ExpAll_NC1 <- QuickProcess(ExpAll_NC1)
ExpAll_NC1 <- FindClusters(ExpAll_NC1, resolution = 3)
DimPlot(ExpAll_NC1, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_NC1, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_NC1, reduction = "umap", label = TRUE, group.by = "big.ident", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_NC1, features = unique(c("Neurod2","Mmp17","Ccer2","Trhr","Chrm2","Cidea","Cpne5","Slc17a7","S100b","Penk","Aga","Prdm8","Gal","Ipcef1","Cxcl13","Kctd4","Luzp2","Itga4","Tafa1","Cartpt","Dab1","Cck","Tent5a","Coch","Stxbp6","Cd44")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
DotPlot(ExpAll_NC1, features = unique(c("Neurod2")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
DotPlot(ExpAll_NC1, features = unique(c("Neurod2","Mmp17","Chrm2","Cidea","Cpne5","Slc17a7","S100b","Penk","Aga","Prdm8","Gal","Ipcef1","Cxcl13","Kctd4","Luzp2","Itga4","Tafa1","Cartpt","Dab1","Cck","Tent5a","Coch","Stxbp6","Cd44")),
        group.by = "newident", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_NC1$newident <- mapvalues(ExpAll_NC1$seurat_clusters, from = c(0,14,13,1,7,23,2,3,16,10,5,17,4,29,6,8,12,15,20,9,28,11,19,27,18,21,24,22,25,26,30,31), to = c(rep("25_Novel",3),rep("26_Novel",3),"14_M",rep("14_ooDS_Cck",5),rep("12_M",2),rep("12_ooDS_NT",5),rep("36_Novel",2),rep("29_Novel",3),rep("35_Novel",2),"20_Novel","26_Novel_M","25_Novel_M","NC1_M","25_Novel_Unk","NC1_M"))
DimPlot(ExpAll_NC1, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_NC1@active.ident <- factor(ExpAll_NC1$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_NC1, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#N
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_N <- subset(ExpAll, idents = c(37,30,25)); ExpAll_N <- QuickProcess(ExpAll_N)
ExpAll_N <- FindClusters(ExpAll_N, resolution = 3)
DimPlot(ExpAll_N, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_N, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_N, reduction = "umap", label = TRUE, group.by = "big.ident", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_N, features = unique(c("Neurod2","Chrm2","Cidea","Slc17a7","S100b","Penk","Aga","Prdm8","Gal","Prdm8","Ipcef1","Itga4","Tafa1")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_N$newident <- mapvalues(ExpAll_N$seurat_clusters, from = c(0,1,4,6,2,3,5), to = c("19_Novel",rep("20_Novel_2",3),"39_Novel","19_Novel_2","N_M"))
DimPlot(ExpAll_N, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_N@active.ident <- factor(ExpAll_N$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_N, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#NC2
ExpAll@active.ident <- factor(ExpAll$seurat_clusters)
ExpAll_NC2 <- subset(ExpAll, idents = c(7,33,22,19,16,18)); ExpAll_NC2 <- QuickProcess(ExpAll_NC2)
ExpAll_NC2 <- FindClusters(ExpAll_NC2, resolution = 2)
DimPlot(ExpAll_NC2, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_NC2, reduction = "umap", label = FALSE, group.by = "hash.ID_broad", pt.size = 1) + theme(aspect.ratio = 1)
DimPlot(ExpAll_NC2, reduction = "umap", label = TRUE, group.by = "big.ident", pt.size = 1) + theme(aspect.ratio = 1)
DotPlot(ExpAll_NC2, features = unique(c("Foxp2","Anxa3","Gria1","Sema5a","Cartpt","Cck","Dab","Coch","Stxbp6","Runx1","Ntrk1","Tent5a","Irx4","Tbr1","Slc7a11","Pcdh20","Calca","Zeb2","Sema3a","Penk","Zic1","Meis2","Col11a1","Apela","Cacna1e","Gm17750")),
        group.by = "seurat_clusters", cols = c("grey", "black"), dot.scale = 6, dot.min = 0,
        cluster.idents = TRUE, scale.by = 'size') + RotatedAxis()
ExpAll_NC2$newident <- mapvalues(ExpAll_NC2$seurat_clusters, from = c(0,18,19,1,4,20,10,11,2,3,14,22,21,5,6,13,7,8,15,16,  9,17,12), to = c(rep("11_Novel_re",3),rep("6_W3B",5),"38_FmidiON_re",rep("15_Novel",4),rep("21_Tbr1_S2",3),rep("13_Novel",4),   "13_M","21_M","11_Novel_M"))
DimPlot(ExpAll_NC2, reduction = "umap", label = TRUE, group.by = "newident", pt.size = 1) + theme(aspect.ratio = 1)
ExpAll@active.ident <- factor(ExpAll$newidents); ExpAll_NC2@active.ident <- factor(ExpAll_NC2$newident)
ExpAll <- Subcluster_ID(ExpAll, ExpAll_NC2, Drop = TRUE);ExpAll$newidents <- ExpAll$SubclusterID

#revised and consolidated some labels
ExpAll$newidents <- mapvalues(ExpAll$newidents, from = c("38_FmidiON_re","3_FminiOFF","Unk1","11_Novel_re","3_FminiOFF_e","27_Novel_1","27_Novel_2","18_Novel_a1","18_Novel_a2","19_Novel_2","25_Novel_Unk"), to = c("38_FmidiON","3_FminiON","44_Novel","11_Novel","3_FminiON_e","27_Novel","27_Novel","18_Novel_a","18_Novel_a","19_Novel","Unk_LQ"))
ExpAll$newidents <- mapvalues(ExpAll$newidents, from = c("20_Novel","20_Novel_2"), to = c("UNK","20_Novel"))

#Filter 4: remove multiplets as defined by subclustering
ExpAll@active.ident <- ExpAll$newidents
# ExpAll <- subset(ExpAll, idents = c("13_M","21_M","N_M","14_M","12_M","NC1_M","Unk_LQ","F_M","T5F_M","TFp_M","ipAE_M","MK_M","ooDS_M","i3_M","26_Novel_M","11_Novel_M","3_FminiON_e","25_Novel_M"), invert = TRUE)
# ExpAll <- DietExport(ExpAll); gc()
# ExpAll <- QuickProcess(ExpAll)
ExpAll$newidents <- droplevels(ExpAll$newidents)
DimPlot(ExpAll, reduction = "umap", label = TRUE, group.by = "newidents") + theme(aspect.ratio = 1)

#Filter 5: remove CMO-defined multiplets
ExpAll@active.ident <- factor(ExpAll$hash.ID)
ExpAll <- subset(ExpAll, idents = "Doublet", invert = TRUE)
ExpAll <- DietExport(ExpAll)
ExpAll <- QuickProcess(ExpAll)

#Fixing some labels
ExpAll$newidents <- mapvalues(ExpAll$newidents, from = "42_AlphaONS", to = "42_AlphaOFFS")
ExpAll$newidents_rev <- mapvalues(ExpAll$newidents, from = c("UNK"), to = c("35_Novel"))
#Note 16_ooDS_D, 16_ooDS_V, 18_Novel_a, and 18_Novel_b, had lower quality cells that weren't properly assigned by louvain clustering and needs to be manually assigned based on marker expression
#example for how to do show is shown below
plot <- DimPlot(ExpAll, reduction = "umap", label = TRUE, group.by = "newidents_rev", pt.size = 0.15, label.size = 5) + theme(aspect.ratio = 1)
c18_Novel_b <- CellSelector(plot, ExpAll$newidents_rev, ident = "18_Novel_b")
ExpAll$newidents_rev[c18_Novel_b] <- "18_Novel_b"

options(future.globals.maxSize=128 * 1000 * 1024^2)

library(scater)
library(scran)
library(DropletUtils)
library(scDblFinder)
library(dplyr)
library(httr)
library(glue)
library(Seurat)
library(stringr)
library(tidyr)

source('utils.R')

set.seed(42)

dir.create(processed.data.folder, showWarnings = FALSE)

load.sample <- function(timepoint, condition, experiment) {
    if (experiment == 1) {
        file.name <- 'GSM9924041_liver-exp1_'
    } else {
        file.name <- 'GSM9924042_liver-exp2_'
    }
    file.name <- paste0(file.name, timepoint)
    if (condition == 'saline') {
        file.name <- paste0(file.name, '_sal_')
    } else {
        file.name <- paste0(file.name, '_pcRBC_')
    }
    file.name <- paste0(file.name, 'matrix_cellbender_filtered.h5')
    path <- file.path(raw.data.folder, file.name)
    print(path)
    sample.sob <- sample.qc(path)
    sample.sob$timepoint <- timepoint
    sample.sob$condition <- condition
    sample.sob$experiment <- experiment
    sample.sob$sample <- paste0(sample.sob$condition, '_', sample.sob$timepoint, '_s', sample.sob$experiment)

    return(sample.sob)
}

if (!file.exists(file.path(processed.data.folder, 'processed_data_all_cells.rds'))) {
    sob.combined <- merge(
        load.sample('2h', 'RBC', 1), 
        y = c(
            load.sample('17h', 'RBC', 1), 

            load.sample('2h', 'saline', 1), 
            load.sample('17h', 'saline', 1),

            load.sample('2h', 'RBC', 2),
            load.sample('17h', 'RBC', 2), 

            load.sample('2h', 'saline', 2), 
            load.sample('17h', 'saline', 2)
        ), 
        add.cell.ids = c(
            'RBC_2h_s1', 
            'RBC_17h_s1', 

            'saline_2h_s1', 
            'saline_17h_s1', 

            'RBC_2h_s2', 
            'RBC_17h_s2', 

            'saline_2h_s2', 
            'saline_17h_s2'
        ), 
        project = 'liver'
    )

    sob.combined <- split(sob.combined, f = sob.combined$sample)
    sob.combined <- SCTransform(sob.combined, vars.to.regress = 'mito.fraction', verbose = TRUE)
    sob.combined <- PrepSCTFindMarkers(object = sob.combined)
    sob.combined <- RunPCA(sob.combined, verbose = FALSE)

    sob.combined.int <- IntegrateLayers(
        object = sob.combined, 
        method = HarmonyIntegration, 
        normalization.method = 'SCT', 
        new.reduction = 'harmony.integrated', 
        verbose = FALSE
    )

    sob.combined.int <- FindNeighbors(sob.combined.int, reduction = 'harmony.integrated', k.param = 20)
    set.seed(42)
    sob.combined.int <- RunUMAP(
        sob.combined.int, 
        dims = 1:23, 
        reduction = 'harmony.integrated', 
        verbose = F, 
        reduction.name = 
        'UMAP.integrated'
    )

    sob.combined.int$sample <- factor(sob.combined.int$sample, levels = c(
        'RBC_2h_s1', 
        'RBC_17h_s1', 
        'RBC_2h_s2', 
        'RBC_17h_s2', 
        'saline_2h_s1', 
        'saline_17h_s1',
        'saline_2h_s2', 
        'saline_17h_s2'
    ))
    sob.combined.int$timepoint <- factor(sob.combined.int$timepoint, levels = c('2h', '17h'))
    sob.combined.int$condition <- factor(sob.combined.int$condition, levels = c('saline', 'RBC'))
    sob.combined.int <- FindClusters(sob.combined.int, resolution = 0.25, algorithm = 4, method = 'igraph')
    Idents(sob.combined.int) <- sob.combined.int$SCT_snn_res.0.25
    
    sob.combined.int <- sob.combined.int[, sob.combined.int$seurat_clusters != '13']

    sob.combined.int$seurat_clusters <- factor(sob.combined.int$seurat_clusters, levels = sort(unique(sob.combined.int$seurat_clusters)))
    Idents(sob.combined.int) <- sob.combined.int$seurat_clusters

    sob.combined.int <- RenameIdents(
        sob.combined.int, 
        '1' = 'T cells 1',
        '2' = 'KCs', 
        '3' = 'T cells 2', 
        '4' = 'Neutrophils.1',
        '5' = 'T cells 3',
        '6' = 'B cells',
        '7' = 'Mono/moMF',
        '8' = 'DCs',
        '9' = 'Neutrophils.2',
        '10' = 'LCMs',
        '11' = 'Cycling',
        '12' = 'pDCs'
    )

    sob.combined.int <- RenameIdents(
        sob.combined.int, 
        'T cells 2' = 'Naive T', 
        'T cells 1' = 'Eff T CD4', 
        'T cells 3' = 'Eff T CD8'
    )

    sob.combined.int$seurat_clusters_annotated <- Idents(sob.combined.int)

    sob.combined.int$seurat_clusters_annotated <- factor(
        sob.combined.int$seurat_clusters_annotated, 
        levels = c('B cells', 'Naive T', 'Eff T CD4', 'Eff T CD8', 'Neutrophils.1', 'Neutrophils.2', 'DCs', 'pDCs', 'Mono/moMF', 'KCs', 'LCMs', 'Cycling')
    )

    Idents(sob.combined.int) <- sob.combined.int$seurat_clusters_annotated
    sob.combined.int$condition.timepoint <- factor(
        paste0(sob.combined.int$condition, '_', sob.combined.int$timepoint), 
        levels = c('saline_2h', 'saline_17h', 'RBC_2h', 'RBC_17h')
    )
    
    saveRDS(sob.combined.int, file.path(processed.data.folder, 'processed_data_all_cells.rds'))
} else {
    if (!file.exists(file.path(processed.data.folder, 'processed_data_myeloid.rds'))) {
        sob.combined.int <- readRDS('processed_data_all_cells.rds')
    }
}

set.seed(42)

if (!file.exists(file.path(processed.data.folder, 'processed_data_myeloid.rds'))) {
    if (!exists("sob.combined.int")) {
        sob.combined.int <- readRDS(file.path(processed.data.folder, 'processed_data_all_cells.rds'))
    }
    
    myeloid <- sob.combined.int[, sob.combined.int$seurat_clusters_annotated %in% c(
        'Neutrophils.1',
        'Neutrophils.2',
        'DCs',
        'pDCs',
        'Mono/moMF',
        'KCs',
        'LCMs'
    )]
    
    set.seed(42)
    myeloid <- RunUMAP(
        myeloid, 
        dims = 1:23, 
        reduction = 'harmony.integrated', 
        verbose = F, 
        reduction.name = 'UMAP.integrated'
    )
    myeloid <- FindClusters(myeloid, resolution = 0.15, algorithm = 4, method = 'igraph')
    myeloid.sub <- FindSubCluster(
        myeloid,
        3,
        graph.name = 'SCT_nn',
        subcluster.name = 'sub.clusters',
        resolution = 0.15,
        algorithm = 4  
    )
    
    ident.levels <- sort(as.numeric(gsub('_', '.', unique(myeloid.sub$sub.clusters))))
    Idents(myeloid.sub) <- factor(as.numeric(gsub('_', '.', myeloid.sub$sub.clusters)), levels = ident.levels)
    myeloid.sub@meta.data$sub.clusters <- Idents(myeloid.sub)
    myeloid.sub@meta.data[['sub.clusters_numeric']] <- Idents(myeloid.sub)

    myeloid.sub@meta.data$sub.clusters.clec4f <- as.numeric(as.character(myeloid.sub@meta.data$sub.clusters_numeric))

    clec4f.indices <- LayerData(myeloid.sub, assay = "SCT", layer = "data")['Clec4f', ] > 0
    myeloid.sub@meta.data[clec4f.indices, 'sub.clusters.clec4f'] <- 'KCs'
   
    myeloid.sub@meta.data[myeloid.sub@meta.data$sub.clusters.clec4f == '1', 'sub.clusters.clec4f'] <- 'moMF.2'
    myeloid.sub@meta.data[myeloid.sub@meta.data$sub.clusters.clec4f == '8', 'sub.clusters.clec4f'] <- 'moMF.2'
    
    Idents(myeloid.sub) <- factor(myeloid.sub@meta.data$sub.clusters.clec4f, levels=sort(unique(myeloid.sub@meta.data$sub.clusters.clec4f)))
    
    myeloid.sub <- RenameIdents(myeloid.sub, '3.1' = 'MonoLy6c-', '3.3' = 'MonoLy6c+', '3.2' = 'moMF.1')
    myeloid.sub <- RenameIdents(myeloid.sub, '4' = 'DCs')
    
    myeloid.sub <- RenameIdents(myeloid.sub, '9' = 'pDCs')
    myeloid.sub <- RenameIdents(myeloid.sub, '2' = 'Neutrophils.1', '7' = 'Neutrophils.2', '5' = 'Neutrophils.3')
    myeloid.sub <- RenameIdents(myeloid.sub, '6' = 'LCMs')

    myeloid.sub$sub.clusters.clec4f.annotated <- factor (
        Idents(myeloid.sub), 
        levels = c('Neutrophils.1', 'Neutrophils.2', 'Neutrophils.3', 'DCs', 'pDCs', 'MonoLy6c-', 'MonoLy6c+', 'moMF.1', 'moMF.2', 'KCs', 'LCMs', 'Cycling')
    )
    Idents(myeloid.sub) <- myeloid.sub$sub.clusters.clec4f.annotated

    myeloid.sub$condition.timepoint <- factor(
        paste0(myeloid.sub$condition, '_', myeloid.sub$timepoint), 
        levels = c('saline_2h', 'saline_17h', 'RBC_2h', 'RBC_17h')
    )
    
    myeloid.sub$condition.timepoint.sal.unif = myeloid.sub$condition.timepoint

    levels(myeloid.sub$condition.timepoint.sal.unif) <- list(
        saline = 'saline_2h', 
        saline = 'saline_17h', 
        RBC_2h = 'RBC_2h', 
        RBC_17h = 'RBC_17h'
    )

    combinations.df <- expand.grid(levels(myeloid.sub$condition.timepoint.sal.unif), levels(Idents(myeloid.sub)))

    myeloid.sub$condition.timepoint.sal.unif.cluster <- factor(paste0(myeloid.sub$condition.timepoint.sal.unif, '_', myeloid.sub$sub.clusters.clec4f.annotated), levels=apply(combinations.df, 1, {function (row) paste0(row[[1]], '_', row[[2]])}))
    
    saveRDS(myeloid.sub, file.path(processed.data.folder, 'processed_data_myeloid.rds'))
}
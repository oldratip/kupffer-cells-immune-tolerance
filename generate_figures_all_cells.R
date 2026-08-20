options(future.globals.maxSize=128 * 1000 * 1024^2)

library(dplyr)
library(httr)
library(ggalluvial)
library(glue)
library(EnhancedVolcano)
library(Seurat)
library(stringr)
library(tidyr)

source('utils.R')

res.folder <- './paper_figures'
dir.create(res.folder, showWarnings = FALSE)

set.seed(42)

sob.combined.int <- readRDS(file.path(processed.data.folder, 'processed_data_all_cells.rds'))

all.clusters.feats <- c(
    'Mafb',
    'Fcgr1',
    'Adgre1',
    'Clec4f',
    'Timd4',
    'Csf1r',
    'Cd5l',
    'Ly6c1',
    'Ly6c2',
    'Itgam',
    'Spn',
    'Irf8',
    'Ly6g',
    'S100a8',
    'S100a9',
    'Itgax',
    'Sirpa',
    'Clec10a',
    'Zeb2',
    'Xcr1',
    'Clec9a',
    'H2-Aa',
    'Ccr7',
    'Tcf4',
    'Siglech',
    'Cd19',
    'Ms4a1',
    'Cd79a',
    'Cd3e',
    'Cd5',
    'Cd28',
    'Cd4',
    'Cd8a',
    'Sell',
    'Cd44',
    'Cd69',
    'Il2ra',
    'Il7r',
    'Klrb1c',
    'Cxcr6',
    'Gzmb',
    'Mki67',
    'Stmn1',
    'Pclaf'
)

plt <- dot.plot.ordered(sob.combined.int, all.clusters.feats)

gene.order <- ggplot_build(plt)$layout$panel_params[[1]]$x$breaks

dot.plot <- DotPlot(sob.combined.int, gene.order, cluster.idents = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    scale_colour_gradient2(low = "#0000FF", mid = "#AAAAAA", high = "#FF0000") + scale_y_discrete(limits=rev)

pdf(
    file.path(res.folder, 'fig_s6c.pdf'),  
    width = 13, 
    height = 5
)
print(dot.plot)
dev.off()

plot <- DimPlot(
    subset(sob.combined.int, condition.timepoint != 'saline_17h'),
    label = FALSE, 
    reduction = 'UMAP.integrated', 
    split.by = 'condition.timepoint',
    cols = palette.all,
    ncol = 3
) + xlab('UMAP_1') + ylab('UMAP_2')
pdf(
    file.path(res.folder, 'fig_2a.pdf'), 
    width = 8, 
    height = 3.5
)
print(plot)
dev.off()

plot <- DimPlot(
    sob.combined.int, 
    label = FALSE, 
    reduction = 'UMAP.integrated', 
    split.by = 'condition.timepoint',
    cols = palette.all,
    ncol = 2
) + xlab('UMAP_1') + ylab('UMAP_2')
pdf(
    file.path(res.folder, 'fig_s6a.pdf'), 
    width = 7, 
    height = 5.5
)
print(plot)
dev.off()

all.markers.filtered <- FindAllMarkers(
    sob.combined.int, 
    min.pct = 0.35, 
    only.pos = TRUE,
    test.use = 'MAST',
    verbose = FALSE,
    recorrect_umi = FALSE, 
    latent.vars = 'experiment'
  )

all.markers.filtered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

all.cluster.markers.heatmap <- DoHeatmap(
    subset(sob.combined.int, downsample = 1000), 
    features = top10$gene,
    angle = 90,
    group.colors = palette.all
) + NoLegend()

pdf(
    file.path(res.folder, 'fig_s6b.pdf'), 
    width = 8, 
    height = 13
)
print(all.cluster.markers.heatmap)
dev.off()

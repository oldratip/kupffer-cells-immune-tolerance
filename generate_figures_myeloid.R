options(future.globals.maxSize=128 * 1000 * 1024^2)

library(dplyr)
library(httr)
library(cowplot)
library(ggalluvial)
library(glue)
library(EnhancedVolcano)
library(Seurat)
library(stringr)
library(tidyr)
library(msigdbr)
library(fgsea)
library(KEGGREST)

source('utils.R')

res.folder <- './paper_figures'
dir.create(res.folder, showWarnings = FALSE)

set.seed(42)

myeloid.sub <- readRDS(file.path(processed.data.folder, 'processed_data_myeloid.rds'))

plot <- DimPlot(
    subset(myeloid.sub, condition.timepoint != 'saline_17h'),
    label = FALSE, 
    reduction = 'UMAP.integrated', 
    split.by = 'condition.timepoint', 
    ncol = 3,
    cols = palette.all
) + xlab('UMAP_1') + ylab('UMAP_2')
pdf(
    file.path(res.folder, 'fig_2b.pdf'),
    width = 8,
    height = 3.5
)
print(plot)
dev.off()

frequencies <- subset(myeloid.sub, condition.timepoint != 'saline_17h')@meta.data %>% 
group_by(condition.timepoint, sub.clusters.clec4f.annotated) %>% summarise(n = n()) %>% mutate(percent.cells = 100 * n / sum(n))
frequencies$sub.clusters.clec4f.annotated <- factor(frequencies$sub.clusters.clec4f.annotated, levels = levels(Idents(myeloid.sub)))

cluster.condition.timepoint.barplot <- ggplot(frequencies, aes(stratum=sub.clusters.clec4f.annotated, alluvium=sub.clusters.clec4f.annotated, fill=sub.clusters.clec4f.annotated, y=percent.cells, x=!!sym('condition.timepoint'))) + 
    xlab('Condition and timepoint') + 
    ylab('Frequency') + 
    geom_flow() +
    geom_stratum(alpha = .8) +
    scale_fill_manual(name = 'Cluster', values = palette.all) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    theme_bw() + theme(
        panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
    )

pdf(
    file.path(res.folder, 'fig_2c.pdf'), 
    width = 4, 
    height = 4.5
)
print(cluster.condition.timepoint.barplot)
dev.off()

myeloid.dotplot.markers <- c(
    'Ly6g',
    'S100a8',
    'S100a9',
    'Cxcr2', 
    'Ly6c1',
    'Ly6c2',
    'Ccr2',
    'Itgam',
    'Spn',
    'F13a1',
    'Cd14',
    'Fcgr1',
    'Fcgr3',
    'Spic',
    'Ccl2',
    'Slc40a1',
    'Adgre1',
    'Clec4f',
    'Timd4',
    'Marco',
    'Mertk',
    'Ms4a7',
    'Apoe',
    'Xcr1',
    'Cxcr3',
    'Itgax',
    'H2-Aa',
    'H2-Eb1',
    'H2-Ab1'
    
)

dot.plot <- DotPlot(myeloid.sub, myeloid.dotplot.markers, cluster.idents = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    scale_colour_gradient2(low = "#0000FF", mid = "#AAAAAA", high = "#FF0000") + 
    scale_y_discrete(limits=rev)

pdf(
    file.path(res.folder, 'fig_2d.pdf'),  
    width = 11, 
    height = 5
)
print(dot.plot)
dev.off()

receptor.clusters <- rev(c(
    'Neutrophils.1',
    'Neutrophils.2',
    'Neutrophils.3',
    'KCs',
    'moMF.1',
    'moMF.2'
))

violin.phagocytes <- function(sob, split.by) {
    p1 <- VlnPlot(sob, 'Hmox1', pt.size = 0, split.by = split.by) + xlab('Cluster')
    p2 <- VlnPlot(sob, 'Blvrb', pt.size = 0, split.by = split.by) + xlab('Cluster')
    p3 <- VlnPlot(sob, 'Slc40a1', pt.size = 0, split.by = split.by) + xlab('Cluster')
    plot <- plot_grid(p1, p2, p3, ncol = 3)
    return(plot)
}

subset.myloid.sub <- subset(myeloid.sub, sub.clusters.clec4f.annotated %in% receptor.clusters)
plot <- violin.phagocytes(subset.myloid.sub, 'condition.timepoint.sal.unif')

pdf(
    file.path(res.folder, 'fig_2e.pdf'), 
    width = 16, 
    height = 3.5
)
print(plot)
dev.off()

receptors <- c(
    'Stab1',
    'Stab2',
    'Mertk',
    'Cd300c2',
    'Sirpa',
    'Sra1',
    'Fcgr1',
    'C1qa',
    'C3',
    'Cd36',
    'Mrc1', 
    'Mfge8',
    'Msr1'
)

plot <- DotPlot(myeloid.sub[, myeloid.sub$sub.clusters.clec4f.annotated %in% receptor.clusters], features = receptors, group.by = c('condition.timepoint.sal.unif.cluster')) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    scale_colour_gradient2(low = "#0000FF", mid = "#AAAAAA", high = "#FF0000") + scale_y_discrete(limits=rev)

for (i in 1:length(receptor.clusters)) {
    plot <- plot + geom_hline(yintercept=3.5 + (i-1) * 3)
}

plot <- plot + scale_y_discrete(limits=rev, labels=gsub('_[^2|17].+$', '', levels(myeloid.sub$condition.timepoint.sal.unif.cluster))) +
    theme(axis.title.y=element_blank(), plot.margin=margin(0, 0, 0, 2, "cm")) +
    coord_cartesian(xlim = c(0.5, length(receptors) + 0.5),
                      clip = 'off')
for (i in 1:length(receptor.clusters)) {
    plot <- plot + annotate("text", x=-3, y=2 + (i-1) * 3, label = gsub('Neutrophils', 'Neutrop', receptor.clusters[[i]]), size = unit(5.2, "pt"), angle=90)
}
pdf(
    file.path(res.folder, 'fig_2f.pdf'), 
    width = 9, 
    height = 7
)
print(plot)
dev.off()


# (ident.group__ident.1__ident.2__subset.ident__filter.condition__filter.timepoint, figure.name)
settings <- list(
    c('condition.timepoint__RBC_2h__saline_2h__KCs______', 'fig_s14a'),
    c('condition.timepoint__RBC_17h__saline_2h__KCs______', 'fig_s14b'),
    c('timepoint__17h__2h__KCs__saline__', 'fig_s14d'),
    c('__KCs__moMF.2____RBC__17h', 'fig_5h')
)

for (s in settings) {
    filename <- s[2]
    setting <- str_split(s[1], '__')[[1]]
    dge.results <- calc.save.DGE(
        myeloid.sub, 
        setting[[2]], 
        setting[[3]], 
        ident.group = if (setting[[1]] != '') setting[[1]] else NULL,
        subset.ident = if (setting[[4]] != '') setting[[4]] else NULL,
        filter.condition = if (setting[[5]] != '') setting[[5]] else NULL,
        filter.timepoint = if (setting[[6]] != '') setting[[6]] else NULL,
    )
    
    title <- paste(paste0(setting[[5]], setting[[6]]), setting[[4]], 'DEGs -', setting[[2]], 'vs', setting[[3]])
    plot <- plot.volcano.dge(dge.results, title)

    pdf(
        file.path(res.folder, paste0(filename, '.pdf')),
        width = 10, 
        height = 10
    )
    print(plot)
    dev.off()
}

dge.res <- per.cluster.DGE(myeloid.sub, 'RBC_17h', 'saline_2h', 'condition.timepoint', c('KCs'))
filtered.res <- dge.res[['KCs']] %>% arrange(avg_log2FC) %>% filter(p_val_adj <= 0.01 & abs(avg_log2FC) >= 1.5)
up <- filtered.res %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj)
down <- filtered.res %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj)
top.up <- rownames(head(up, 35))
top.up <- top.up[top.up %in% rownames(myeloid.sub[['SCT']]@scale.data)]
top.down <- rownames(head(down, 35))
top.down <- top.down[top.down %in% rownames(myeloid.sub[['SCT']]@scale.data)]
all.feats <- c(top.up, top.down)

sob.subset <- subset(myeloid.sub, (condition.timepoint == 'RBC_2h' | condition.timepoint == 'RBC_17h' | condition.timepoint == 'saline_2h') & sub.clusters.clec4f.annotated == 'KCs')
sob.subset$condition.timepoint <- factor(sob.subset$condition.timepoint, levels = c('saline_2h', 'RBC_2h', 'RBC_17h'))

Idents(sob.subset) <- sob.subset@meta.data[['condition.timepoint']]
sob.subset <- subset(sob.subset, downsample = 500)

plot <- DoHeatmap(
    sob.subset,
    c(top.up, top.down), 
    slot = 'scale.data',
    group.by = 'condition.timepoint'
) + geom_hline(yintercept=length(top.down)+0.5, color = "white", linewidth = 1.5)

comparison.level <- levels(sob.subset$condition.timepoint)[[length(levels(sob.subset$condition.timepoint))]]
additional.dge.res <- calc.save.DGE(myeloid.sub, comparison.level, 'RBC_2h', subset.ident = 'KCs', ident.group = 'condition.timepoint')

additional.dge.res <- additional.dge.res %>% arrange(avg_log2FC) %>% filter(p_val_adj <= 0.01 & abs(avg_log2FC) >= 1.15)
additional.dge.res.up <- additional.dge.res %>% filter(avg_log2FC > 0)
additional.dge.res.down <- additional.dge.res %>% filter(avg_log2FC < 0)
significant.up <- ifelse(top.up %in% rownames(additional.dge.res.up), "red", "black")
significant.down <- ifelse(top.down %in% rownames(additional.dge.res.down), "blue", "black")
plot <- plot + theme(axis.text.y = element_text(colour = rev(c(significant.up, significant.down))))

pdf(
    file.path(res.folder, 'fig_4c.pdf'),
    width = 9, 
    height = length(c(top.up, top.down)) * 13 / 70
)
print(plot)
dev.off()
                                                     
msigdb.mm <- msigdbr(db_species = "MM", species = "Mus musculus")
                                
interesting.pathways <- c(
    'COATES_MACROPHAGE_M1_VS_M2_UP',
    'COATES_MACROPHAGE_M1_VS_M2_DN',
    'GOBP_INTRACELLULAR_IRON_ION_HOMEOSTASIS',
    'GOBP_CELLULAR_RESPONSE_TO_IRON_ION',
    'GOBP_TOLERANCE_INDUCTION',
    'BIOCARTA_IL10_PATHWAY',
    'GOBP_HEME_METABOLIC_PROCESS',
    'REACTOME_HEME_SIGNALING',
    'GOBP_HEME_TRANSPORT',
    'WP_APOPTOSIS',
    'GOBP_FERROPTOSIS',
    'GOBP_PHAGOCYTOSIS_ENGULFMENT',
    'GOBP_PHAGOCYTOSIS',
    'GOMF_FERROUS_IRON_TRANSMEMBRANE_TRANSPORTER_ACTIVITY',
    'GOBP_PROSTAGLANDIN_SECRETION',
    'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',
    'GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION',
    'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
    'WP_INFLAMMATORY_RESPONSE_PATHWAY',
    'GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
    'GOBP_POSITIVE_REGULATION_OF_MYELOID_LEUKOCYTE_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE'
)
     
filtered.sets <- msigdb.mm %>% filter(gs_name %in% interesting.pathways | grepl('HALLMARK', gs_name))

kegg.ferroptosis.id <- "mmu04216"
p <- keggGet(kegg.ferroptosis.id)
kegg.ferroptosis.genes <- p[[1]]$GENE
kegg.ferroptosis.genes <- gsub(';.+', '', kegg.ferroptosis.genes[seq(2, length(kegg.ferroptosis.genes), 2)])

additional.gene.sets <- list(
    list(
        rep('KEGG_FERROPTOSIS', length(kegg.ferroptosis.genes))
    ), 
    list(
        kegg.ferroptosis.genes
    )
)                                                    
                                                     
# (initial.subset__initial.subset.val__ident.1__ident.2__group.by__subset.ident, figure.name)
settings <- list(
    c('*__*__RBC_2h__saline_2h__condition.timepoint__KCs', 'fig_4d1'),
    c('*__*__RBC_17h__saline_2h__condition.timepoint__KCs', 'fig_4d2'),
    c('condition.timepoint__RBC_17h__KCs__moMF.1__*__*', 'fig_5j'),
    c('condition.timepoint__RBC_17h__KCs__moMF.2__*__*', 'fig_5i')
)            
                                                     
for (s in settings) {
    filename <- s[2]
    setting <- str_split(s[1], '__')[[1]]
    pvalue.limit <- c(0, 0.05)
    
    subset.obj <- myeloid.sub
    if (setting[[1]] != '*') {
        subset.obj <- subset(myeloid.sub, !!sym(setting[[1]]) == setting[[2]])
    }
    
    gsea.result <- run.gsea(
        subset.obj,
        setting[[3]], 
        setting[[4]],
        group.by = if (setting[[5]] == '*') NULL else setting[[5]],
        subset.ident = if (setting[[6]] == '*') NULL else setting[[6]],
    )
    title <- paste(gsub('RBC', 'pcRBC', gsub('_', ' ', setting[[2]])), gsub('\\*', '-', setting[[6]]), setting[[3]], 'vs.', setting[[4]])
    plot <- plot.gsea(gsea.result, title, paste0(filename, '.pdf'), pvalue.limit)
    plot
}
        
alb.plot.all <- FeaturePlot(myeloid.sub, c('Alb'), reduction = 'UMAP.integrated', label = TRUE, order = TRUE)
alb.plot.all <- alb.plot.all + ggtitle('Alb - all samples')

pdf(
    file.path(res.folder, 'fig_s14c.pdf'), 
    width = 7, 
    height = 6
)
print(alb.plot.all)
dev.off()

markers.filtered <- FindAllMarkers(
    myeloid.sub, 
    min.pct = 0.35, 
    only.pos = TRUE,
    test.use = 'MAST',
    verbose = FALSE,
    recorrect_umi = FALSE, 
    latent.vars = 'experiment'
)

markers.filtered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

cluster.markers.heatmap <- DoHeatmap(
    subset(myeloid.sub, downsample = 1000), 
    features = top10$gene,
    angle = 90,
    group.colors = palette.all
) + NoLegend()

pdf(
    file.path(res.folder, 'fig_s6d.pdf'), 
    width = 7.2, 
    height = 11
)
print(cluster.markers.heatmap + theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm")))
dev.off()

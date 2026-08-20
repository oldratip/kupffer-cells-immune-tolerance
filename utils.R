library(RColorBrewer)

raw.data.folder <- './raw_data'
processed.data.folder <- './processed_data'

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], 'Set3')
general.palette <- f("Set3")
general.palette[[2]] <- '#d7d7a5'
general.palette[[7]] <- '#9bc15a'
palette.all <- rev(general.palette)
palette.myeloid <- general.palette

sample.qc <- function(file) {
    sce <- read10xCounts(file, col.names=TRUE)

    sce <- scDblFinder(sce)
    singlets <- sce[, sce$scDblFinder.class == 'singlet']
    
    m <- counts(singlets)
    umi.sum <- Matrix::colSums(m)
    umi.criteria <- isOutlier(umi.sum, type='lower', log=TRUE) 
    print(paste('Removing', sum(umi.criteria), 'cells due to number of UMIs.'))

    detected.criteria <- isOutlier(Matrix::colSums(m > 1), type='lower', log=TRUE) 
    print(paste('Removing', sum(detected.criteria), 'cells due to number of detected genes.'))

    mito.fraction <- Matrix::colSums(m[grep('^mt-', rowData(singlets)$Symbol), ]) / umi.sum
    mito.criteria <- isOutlier(mito.fraction, type='higher')
    print(paste('Removing', sum(mito.criteria), 'cells due to fraction of mito genes.'))

    colData(singlets)$mito.fraction <- mito.fraction

    overall.criteria <- !umi.criteria & !umi.criteria & !mito.criteria
    singlets.qcd <- singlets[, overall.criteria]
    print(paste('Removing a total of', ncol(m) - sum(overall.criteria), 'cells.'))
    print(paste('Remaining cells:', ncol(singlets.qcd)))
    
    rownames(singlets.qcd) <- rowData(singlets.qcd)$Symbol
        
    sob <- RenameAssays(object = as.Seurat(singlets.qcd, data = NULL), originalexp = 'RNA')
    
    return(sob)
}

dot.plot.ordered <- function(sob, features) {
    options(repr.plot.width = 13.5, repr.plot.height = 5)
    
    dot.plot <- DotPlot(sob, features) + 
        theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
        scale_colour_gradient2(low = "#0000FF", mid = "#AAAAAA", high = "#FF0000") +
         scale_y_discrete(limits=rev)
    
    dp.data <- dot.plot$data

    wide.data <- dp.data %>%
      select(features.plot, id, avg.exp.scaled) %>%
      pivot_wider(names_from = id, values_from = avg.exp.scaled)

    gene.matrix <- as.matrix(wide.data[,-1])
    rownames(gene.matrix) <- wide.data$features.plot
    dend <- hclust(dist(gene.matrix))

    gene.order <- dend$labels[dend$order]

    dot.plot <- DotPlot(sob, features = gene.order) + 
        theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
        scale_colour_gradient2(low = "#0000FF", mid = "#AAAAAA", high = "#FF0000") +
         scale_y_discrete(limits=rev)

    return(dot.plot)

}

calc.save.DGE <- function(
    sob,
    ident.1, 
    ident.2 = NULL, 
    ident.group = NULL, 
    subset.ident = NULL, 
    filter.condition = NULL, 
    filter.timepoint = NULL,
    min.pct = 0.15,
    verbose = FALSE
) {    
    if (!is.null(filter.condition)) {
        sob <- subset(sob, condition == filter.condition)
    }
    
    if (!is.null(filter.timepoint)) {
        sob <- subset(sob, timepoint == filter.timepoint)
    }
    
    markers <- FindMarkers(
        sob,
        ident.1 = ident.1, 
        ident.2 = ident.2, 
        subset.ident = subset.ident, 
        group.by = ident.group, 
        test.use = 'MAST', 
        verbose = verbose,
        recorrect_umi = FALSE, 
        latent.vars = 'experiment',
        min.pct = min.pct
    )
    
    return(markers)
}

filter.markers <- function(markers, p.val.thresh = 0.01, fc.thresh = 1.5) {
    markers <- markers %>% 
        filter(p_val_adj < p.val.thresh) %>% 
        arrange(desc(avg_log2FC)) %>% 
        filter(abs(avg_log2FC) >= fc.thresh)
    
    return(markers)
}

plot.volcano.dge <- function(dge.res, title, selectLab = NULL) {
    point.size = c(ifelse(rownames(dge.res) == 'Marco', 4, 2))

    plot <- EnhancedVolcano(
        dge.res,
        lab = rownames(dge.res),
        x = 'avg_log2FC',
        FCcutoff = 1.5,
        pCutoff = 0.01,
        pointSize = 3,
        drawConnectors = TRUE,
        widthConnectors = 0.75,
        subtitle = '',
        caption = NULL,
        title = title,
        selectLab = selectLab,
        y = 'p_val_adj'
    )
    
    return(plot)
}

per.cluster.DGE <- function(sob, ident.1, ident.2, ident.group, interesting.clusters = NULL) {
    results = list()
    if (is.null(interesting.clusters)) {
        interesting.clusters <- unique(Idents(sob))
    }
    for (cluster in interesting.clusters) {
        results[[cluster]] <- calc.save.DGE(sob, ident.1, ident.2, subset.ident = cluster, ident.group = ident.group)
    }
    return(results)
}


plot.gsea <- function(result, name, filename, custom.limits = NULL) {
    df <- result$res %>% arrange(NES) %>% filter(padj < 0.05)
    if (nrow(df) == 0) {
        print(paste('No significantly enriched gene sets for', name))
        return()
    }
    df$pathway <- factor(df$pathway, levels = df$pathway)


    if (is.null(custom.limits)) {
        custom.limits <- c(min(df$padj), max(df$padj))
    }
    
    df <- result$res %>% arrange(-log10(padj)) %>% filter(padj < 0.05)
    df$pathway <- factor(df$pathway, levels = df$pathway)
    
    df.pos <- df[df$NES > 0, ] %>% arrange(desc(NES))
    df.neg <- df[df$NES < 0, ] %>% arrange(desc(abs(NES)))

    df <- rbind(df.pos, df.neg)
    df$pathway <- factor(df$pathway, levels = rev(c(as.character(df.pos$pathway), as.character(df.neg$pathway))))

    reordered <- ggplot(df,
       aes(
           x = NES,
           y = pathway
        )
    ) + geom_col(aes(fill=padj), width = 0.8) + ggtitle(paste(name, '- enriched')) +
        theme_bw() + theme(
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
    ) + ylab('Pathway') 
    
    reordered$data$label <- sub('WP_APOPTOSIS', 'APOPTOSIS (WP)', reordered$data$pathway)
    reordered$data$label <- sub('FERROUS_IRON_TRANSMEMBRANE_TRANSPORTER_ACTIVITY', 'FERROUS_IRON_TRANSMEMBRANE_TRANSP.', sub('^[^_]+_', '', reordered$data$label))
    reordered$data$label <- sub('REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION', 'REGULATION_OF_ANTIGEN_PROC_AND_PRES.', reordered$data$label)
    reordered$data$label <- sub('LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE', 'LEUK_ACTIV_INVOLVED_IN_INFLAMM_RESP', reordered$data$label)

    if (is.null(custom.limits)) {
        custom.limits <- c(min(df$padj), max(df$padj))
    }
    reordered <- reordered + scale_y_discrete(labels=sub('^[^_]+_', '', reordered$data$pathway)) + 
        scale_fill_gradient(name = 'Adjusted p', low = 'red', high = 'blue', limits = custom.limits) + 
        theme(legend.position = "right") +
        xlim(-2.5, 2.5) +
        geom_text(aes(label = label,
                    x = ifelse(NES > 0, -0.05, 0.05),
                    hjust = ifelse(NES > 0, 1, 0),
                    vjust = 0.5),
                size = 3) + theme_bw() + theme(, 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                           )

    pdf(
        file.path(res.folder, paste0(filename, '.pdf')), 
        width = 7.5, 
        height = 1 + nrow(df) * 0.2 
    )
    print(reordered)
    dev.off()
    
    return(df)
}

run.gsea <- function(data, ident.1, ident.2, group.by = NULL, subset.ident = NULL, top = 100, fdr.thresh = 0.05) {
    message(paste('Running GSEA between', ident.1, ident.2, 'grouped by', group.by, 'subsetted on', subset.ident))
    degs <- FindMarkers(
        data,
        ident.1 = ident.1,
        ident.2 = ident.2,
        min.pct = 0,
        group.by = group.by,
        subset.ident = subset.ident,
        logfc.threshold = 0.0,
        test.use = 'MAST', 
        verbose = FALSE,
        recorrect_umi = FALSE,
        latent.vars = 'experiment'
    )
    
    degs$gene <- rownames(degs) 
    degs$p_val_nonzero <- degs$p_val
    if (nrow(degs[degs$p_val_nonzero == 0, ]) > 0) {
        degs[degs$p_val_nonzero == 0, ]$p_val_nonzero <- .Machine$double.xmin
    }
        
    degs$fcsign <- sign(degs$avg_log2FC)
    degs$logP <- -log10(degs$p_val_nonzero)
    degs$metric <- degs$logP * degs$fcsign    
    
    sorted.degs <- degs %>%
        dplyr::arrange(desc(metric)) %>%
        dplyr::select(gene, metric)

    ranks <- tibble::deframe(sorted.degs)
    
    msigdbr.list <- split(x = c(filtered.sets$gene_symbol, unlist(additional.gene.sets[[2]])), f = c(filtered.sets$gs_name, unlist(additional.gene.sets[[1]])))
    fgsea.res <- fgsea(pathways = msigdbr.list, stats = ranks)
    
    top.pathways.up <- fgsea.res[ES > 0 & padj < fdr.thresh][head(order(pval), n=top), pathway]
    top.pathways.down <- fgsea.res[ES < 0 & padj < fdr.thresh][head(order(pval), n=top), pathway]
    top.pathways <- c(top.pathways.up, rev(top.pathways.down))
    
    return(list(top.pathways=top.pathways, res=fgsea.res))
}

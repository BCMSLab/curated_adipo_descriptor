library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(EnrichedHeatmap)

WGCNA::allowWGCNAThreads(4)

peak_counts <- read_rds('data/peak_counts.rds')

hm <- c(H3K27ac_0 = 'GSM1412512',
        H3K27ac_1 = 'GSM1370470',
        H3K27ac_3 = 'GSM535765',
        H3K4me1_0 = 'GSM2515940',
        H3K4me1_1 = 'GSM1370471',
        H3K4me1_3 = 'GSM535764',
        H3K4me3_0 = 'GSM1017630',
        H3K4me3_1 = 'GSM535755',
        H3K4me3_3 = 'GSM535762')
fls <- paste0('data/bws/', hm, '_EF.bw')
names(fls) <- names(hm)
all(file.exists(fls))

goi <- AnnotationDbi::select(org.Mm.eg.db,
                             c('GO:0045599', 'GO:0045600'),
                             c('SYMBOL', 'ENTREZID'),
                             'GO')

proms <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,  
                   filter = list(gene_id = unique(goi$ENTREZID)),
                   upstream = 3000,
                   downstream = 3000)

covr <- map(fls, function(x) import.bw(x, selection = BigWigSelection(proms)))

tss <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,  
                 filter = list(gene_id = unique(goi$ENTREZID)),
                 upstream = 0,
                 downstream = 1,
                 columns = c('tx_id', 'gene_id'))

mats <- map(covr, function(x) {
  normalizeToMatrix(x,
                    value_column = 'score',
                    target = tss,
                    extend = 3000,
                    mean_mode = 'w0',
                    w = 10)
})
# smooth = true

group <- goi$GO[match(as.character(tss$gene_id), as.character(goi$ENTREZID))]
ta <- HeatmapAnnotation(
  enriched = anno_enriched(
    gp = gpar(col = 2:4, lty = 1:3),
    axis = FALSE
  )
)

hm_list <- NULL
for (s in 1:length(mats)) {
  x <- mats[[s]]
  y <- names(mats)[s]
  
  h <- EnrichedHeatmap(x,
                       column_title = y,
                       col = c("white", "red"),
                       top_annotation = ta,
                       row_split = group,
                       cluster_rows = TRUE,
                       axis_name = c("-3kb", "TSS", "3kb"))
  hm_list <- hm_list + h
}

png(filename = 'manuscript/figures/hm_coverage.png',
    height = 15, width = 30, units = 'cm', res = 300)
draw(hm_list,
     show_heatmap_legend = FALSE,
     ht_gap = unit('1', 'mm'))
dev.off()

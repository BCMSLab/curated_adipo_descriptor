# loading requried librars
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(ComplexHeatmap)

# loading data
adipo_count <- read_rds('data/adipo_counts.rds')
se <- adipo_count[, order(adipo_count$time)]

dds <- DESeqDataSet(se,design = ~stage)
dds_transform <- vst(dds)

d <- as.matrix(dist(t(assay(se))))

ra <- rowAnnotation(Factor1 = anno_mark(at = which(!duplicated(se$time)),
                                        labels = unique(se$time)))
ca <- columnAnnotation(Factor2 = anno_mark(at = which(!duplicated(se$time)),
                                           labels = unique(se$time)))

hm1 <- Heatmap(d,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE,
        top_annotation = ca,
        right_annotation = ra,
        column_split = se$stage)

# chip
peak_counts <- read_rds('data/peak_counts.rds')

# select samples for the factor
sample_ind <- (peak_counts$factor %in% c('CEBPB', 'PPARG')) & (!is.na(peak_counts$factor))
sample_ids <- colnames(peak_counts)[sample_ind]

# select peaks from the selected samples
peak_ind <- lapply(mcols(peak_counts)$peak, function(x) sum(sample_ids %in% x))
peak_ind <- unlist(peak_ind) > 2

# subset the object
se_chip <- peak_counts[peak_ind, sample_ind]
se_chip <- se_chip[, !is.na(se_chip$time)]
se_chip <- se_chip[, order(se_chip$factor, se_chip$time)]

dds_chip <- DESeqDataSet(se_chip,design = ~stage)
dds_chip_transform <- vst(dds_chip)

d_chip <- as.matrix(dist(t(assay(dds_chip_transform))))

ra <- rowAnnotation(Factor1 = anno_mark(at = which(!duplicated(se_chip$time)),
                                        labels = unique(se_chip$time)))
ca <- columnAnnotation(Factor2 = anno_mark(at = which(!duplicated(se_chip$time)),
                                           labels = unique(se_chip$time)))

hm2 <- Heatmap(d_chip,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE,
        top_annotation = ca,
        right_annotation = ra,
        column_split = se_chip$factor)

png(filename = 'manuscript/figures/replicates_similarity.png',
    width = 18, height = 9, units = 'cm', res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1,
                      layout.pos.row = 1))
grid.text("A", x = .02, y = .97)

draw(hm1, newpage = FALSE, padding = unit(rep(.5, 4), 'cm'))

upViewport()

pushViewport(viewport(layout.pos.col = 2,
                      layout.pos.row = 1))
grid.text("B", x = 0, y = .97)
draw(hm2, newpage = FALSE, padding = unit(rep(.5, 4), 'cm'))
upViewport()

dev.off()

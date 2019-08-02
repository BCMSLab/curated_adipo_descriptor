# loading requried librars
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# loading data
adipo_count <- read_rds('data/adipo_counts.rds')
pd <- tibble(time = adipo_count$time,
             stage = adipo_count$stage,
             study = adipo_count$study)

dds <- DESeqDataSet(adipo_count,design = ~stage)
dds_transform <- vst(dds)
mat <- assay(dds_transform)

markers <- list('Adipogenic' = c('Pparg', 'Cebpb', 'Cebpa'),
                'Lipogenic' = c('Lpl', 'Acly', 'Fasn'))

col_fun <- colorRamp2(c(0, 20), c('white', 'darkblue'))

ind <- rownames(mat) %in% unlist(markers)
hm1 <- Heatmap(mat[ind,],
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_split = dds_transform$stage,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)

# chip
peak_counts <- read_rds('data/peak_counts.rds')

# select samples for the factor
sample_ind <- (peak_counts$factor %in% c('CEBPB', 'PPARG')) & (!is.na(peak_counts$factor))
sample_ids <- colnames(peak_counts)[sample_ind]

# select peaks from the selected samples
peak_ind <- lapply(mcols(peak_counts)$peak, function(x) sum(sample_ids %in% x))
peak_ind <- unlist(peak_ind) > 2

# subset the object
se <- peak_counts[peak_ind, sample_ind]
#se$factor <- ifelse(is.na(se$factor), 'Input', se$factor)

pd_chip <- tibble(time = se$time,
                  stage = se$stage,
                  study = se$study,
                  factor = se$factor)

dds_chip <- DESeqDataSet(se,design = ~factor)
dds_chip_transform <- vst(dds_chip)
mat_chip <- assay(dds_chip_transform)

targets <- c('Acly', 'Fasn', 'Lpl')
ind <- mcols(dds_chip)$geneId %in% targets

hm2 <- Heatmap(mat_chip[ind,],
        row_labels = ifelse(duplicated(mcols(dds_chip)$geneId[ind]),
                            '',
                            mcols(dds_chip)$geneId[ind]),
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        cluster_rows = FALSE,
        column_split = dds_chip$factor,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)

png(filename = 'manuscript/figures/markers.png',
    width = 18, height = 9, units = 'cm', res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1,
                      layout.pos.row = 1))
grid.text("A", x = .02, y = .97)

draw(hm1, newpage = FALSE)

upViewport()

pushViewport(viewport(layout.pos.col = 2,
                      layout.pos.row = 1))
grid.text("B", x = 0, y = .97)
draw(hm2, newpage = FALSE)
upViewport()

dev.off()

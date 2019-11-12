# load required libraries
library(GEOquery)
library(tidyverse)
library(reshape2)
library(WGCNA)
library(sva)
library(rtracklayer)
library(ChIPseeker)
library(ComplexHeatmap)
library(SummarizedExperiment)

allowWGCNAThreads(4)

# if (!file.exists('data/GSE98680_series_matrix.txt.gz')) {
#   getGEO('GSE98680',
#          destdir = 'data/')
# }
# 
# eset <- getGEO('GSE98680', destdir = 'data/')[[1]][, 1:24]
# 
# mat <- collapseRows(exprs(eset),
#                     rowID = featureNames(eset),
#                     rowGroup = fData(eset)$GENE_SYMBOL)[[1]]
# 
# newset <- ExpressionSet(mat)
# 
# newset$group <- str_split(eset$title, '_', simplify = TRUE)[, 3]
# newset$patient <- str_split(eset$title, '_', simplify = TRUE)[, 2]
# 
# exprs(newset) <- ComBat(exprs(newset),
#                         batch = newset$patient,
#                         mod = model.matrix(~group, data = pData(newset)))
# 
# write_rds(newset, 'data/primary_adipocyte.rds')
#
# if(!file.exists('data/GSE68864_RAW.tar')) {
#   download.file(
#     'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68864&format=file',
#     destfile = 'GSE68864_RAW.tar'
#   )
#   untar('data/GSE68864_RAW.tar',
#         exdir = 'GSE68864_RAW')
#   file.copy('data/GSE68864_RAW/GSM1684635_CEBPb.HDDMI.bed.gz',
#             'data/CEBPB_hMSC_6h.bed.gz')
# }
# 
# if(!file.exists('data/PPARG_hMADS_19d.bed.gz')) {
#   download.file(
#     'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE59703&format=file&file=GSE59703%5FhMADSPPARgPeakfile%2Ebed%2Egz',
#     destfile = 'data/PPARG_hMADS_19d.bed.gz'
#   )
# }
#
# peaks <- c('data/CEBPB_hMSC_6h.bed.gz', 'data/PPARG_hMADS_19d.bed.gz') %>%
#   map(function(x) {
#     import.bed(x) %>% 
#       annotatePeak(level = 'gene',
#                    verbose = FALSE,
#                    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#                    annoDb = 'org.Hs.eg.db') %>%
#       as.GRanges()
#   })
# 
# write_rds(peaks, 'data/primary_adipocyte_chip.rds')

targets <- list(Adipogenic = c('PPARG', 'CEBPB', 'CEBPA'),
                Lipogenic = c('LPL', 'ACLY', 'FASN'))

primary_adipocyte <- read_rds('data/primary_adipocyte.rds')
se1 <- primary_adipocyte[rownames(primary_adipocyte) %in% unlist(targets)]
mat1 <- exprs(se1)[unlist(targets), ]

marker1 <- columnAnnotation(cell = anno_mark(
  at = c(2, 15),
  labels = c('Primary Cell\n None-differentiated\n Pre-adipocyte',
             'Primary Cell\n Differentiating\n Adipocyte')
))

hm1 <- Heatmap(mat1,
               top_annotation = marker1,
               column_title = rep('', 2),
               show_column_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_split = se1$group,
               row_split = rep(names(targets), lengths(targets)),
               show_heatmap_legend = FALSE)
hm1

primary_adipocyte_chip <- read_rds('data/primary_adipocyte_chip.rds') %>%
  set_names(c('CEBPB', 'PPARG')) %>%
  map(function(x) {
    mat <- as_tibble(x) %>%
      filter(SYMBOL %in% unlist(targets)) %>%
      mutate(annotation = str_split(annotation, ' ', simplify = TRUE)[, 1]) %>%
      filter(annotation %in% c('Distal', 'Exon', 'Intron', 'Promoter')) %>%
      acast(SYMBOL ~ annotation)
    mat[mat > 0] <- 1
    mat[unlist(targets),]
  })

mat2 <- cbind(primary_adipocyte_chip$CEBPB, primary_adipocyte_chip$PPARG)
marker2 <- columnAnnotation(cell = anno_mark(at = c(1, 5),
                                             labels = c('CEBPB\n hMSC\n + MDI (6 hours)',
                                                        'PPARG\n hMADS\n + MDI (19 days)')))
hm2 <- Heatmap(mat2,
               top_annotation = marker2,
               col = c('white', 'black'),
               cluster_rows = FALSE,
               cluster_row_slices = FALSE,
               row_split = rep(names(targets), lengths(targets)),
               cluster_columns = FALSE,
               cluster_column_slices = FALSE,
               column_split = rep(c('CEBPB', 'PPARG'), each = 4),
               column_title = rep('', 2),
               show_heatmap_legend = FALSE,
               column_names_rot = 45)


png(filename = 'manuscript/figures/primary_adipocytes.png',
    width = 20, height = 13, units = 'cm', res = 300)
hm2 + hm1
dev.off()

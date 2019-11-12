# loading requried librars
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(goseq)
library(org.Mm.eg.db)
library(cowplot)

# loading data
adipo_count <- read_rds('data/adipo_counts.rds')

low_counts <- apply(assay(adipo_count), 1, function(x) length(x[x>10])>=2)
se <- adipo_count[low_counts,]
se$stage <- factor(se$stage)

dds <- DESeqDataSet(se, design = ~stage)
dds <- DESeq(dds)
res <- map(resultsNames(dds)[-1],
           function(x) {
             results(dds,
                     name = x,
                     tidy = TRUE)
           }) %>%
  bind_rows(.id = 'contrast') %>%
  as_tibble() %>%
  na.omit()

deg <- split(res, res$contrast) %>%
  map(function(x) {
    deg <- as.integer((x$padj < .2) & (abs(x$log2FoldChange) > 1))
    names(deg) <- x$row
    deg
  })

# prepare gene to gene ontology terms
gene_go <- AnnotationDbi::select(org.Mm.eg.db,
                                 unique(res$row),
                                 'GO',
                                 'SYMBOL') %>%
  dplyr::select(gene = SYMBOL, cat = GO) %>%
  unique()

# perform gene erichment analysis
deg_go_res <- map(deg, function(y) {
  nullp(y, 'mm10', 'geneSymbol', plot.fit = FALSE) %>%
    goseq(gene2cat = gene_go)
}) %>%
  bind_rows(.id = 'contrast') %>%
  as_tibble() 

go <- c('GO:0060612', 'GO:0019915', 'GO:0032869', 'GO:0016042', 'GO:0006006')

p1 <- deg_go_res %>%
  filter(category %in% go) %>%
  mutate(fraction = numDEInCat/numInCat,
         label = ifelse(over_represented_pvalue < .05, '*', '')) %>%
  ggplot(aes(x = category, y = fraction)) +
  geom_col() +
  geom_text(aes(x = category, y = fraction + .05, label = label)) +
  lims(y = c(0, 1)) +
  facet_wrap(~contrast, nrow = 1) +
  theme_bw() +
  labs(x = '', y = 'Fraction of DE genes') +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# chip
peak_counts <- read_rds('data/peak_counts.rds')

# select samples for the factor
dds <- map(list(PPARG = 'PPARG', CEBPB = 'CEBPB'),
           function(x) {
             sample_ind <- (peak_counts$factor == x) & (!is.na(peak_counts$factor))
             sample_ids <- colnames(peak_counts)[sample_ind]
             
             # select peaks from the selected samples
             peak_ind <- lapply(mcols(peak_counts)$peak, function(x) sum(sample_ids %in% x))
             peak_ind <- unlist(peak_ind) > 2
             
             # subset the object
             se <- peak_counts[peak_ind, sample_ind]
             se$stage <- factor(se$stage)
             d <- DESeqDataSet(se,design = ~stage)
             DESeq(d)
           })

res <- list()
res$PPARG <- results(dds$PPARG, name = 'stage_3_vs_0', tidy = TRUE)
res$CEBPB <- results(dds$CEBPB, name = 'stage_1_vs_0', tidy = TRUE)

res <- bind_rows(res, .id = 'factor') %>% as_tibble()

gene_go_res <- tibble(gene = as.character(mcols(peak_counts)$geneId),
                      row = rownames(peak_counts)) %>%
  right_join(gene_go) %>%
  right_join(res) %>%
  na.omit()

factor_labels <- c('CEBPB (Stage 1 vs 0)', 'PPARG (Stage 3 vs 0)')
names(factor_labels) <- c('CEBPB', 'PPARG')

set.seed(1223)
random_res <- gene_go_res[sample(1:nrow(gene_go_res), 50),]
random_res$cat <- '(Random Set)'

p2 <- gene_go_res %>%
  filter(cat %in% go,
         padj < .2) %>%
  rbind(random_res) %>%
  ggplot(aes(x = cat, y = abs(log2FoldChange))) +
  geom_boxplot() +
  facet_wrap(~factor, labeller = labeller(factor = factor_labels)) +
  theme_bw() +
  labs(x = '', y = 'Absolute fold-change (Log2)') +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"),
        axis.text.x = element_text(angle = 45, hjust = 1))

plot_grid(p1, p2, 
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/go_enrichment.png',
         height = 9, width =20, units = 'cm')

# loading requried librars
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(cowplot)

# loading data
adipo_count <- read_rds('data/adipo_counts.rds')
pd <- tibble(time = adipo_count$time,
             stage = adipo_count$stage,
             study = adipo_count$study)

dds <- DESeqDataSet(adipo_count,design = ~stage)
dds_transform <- vst(dds)

mat <- assay(dds_transform)

dims <- cmdscale(dist(t(mat))) %>%
  as_tibble() %>%
  setNames(c('dim1', 'dim2')) %>%
  bind_cols(pd)

p1 <- dims %>%
  ggplot(aes(x = dim1, y = dim2, color = as.factor(study)), alpha = .5) +
  geom_text(aes(label = as.factor(stage)), alpha = .5) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'Dim1', y = 'Dim 2')

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

pd_chip <- tibble(time = se$time,
             stage = se$stage,
             study = se$study,
             factor = se$factor)

dds_chip <- DESeqDataSet(se,design = ~stage)
dds_chip_transform <- vst(dds_chip)

mat_chip <- assay(dds_chip_transform)

dims_chip <- cmdscale(dist(t(mat_chip))) %>%
  as_tibble() %>%
  setNames(c('dim1', 'dim2')) %>%
  bind_cols(pd_chip)

p2 <- dims_chip %>%
  ggplot(aes(x = dim1, y = dim2, color = as.factor(study)), alpha = .5) +
  geom_text(aes(label = as.numeric(as.factor(factor))), alpha = .5) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'Dim1', y = 'Dim 2')

plot_grid(p1, p2,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/mds.png',
         height = 8, width = 17, units = 'cm')

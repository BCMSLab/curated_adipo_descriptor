# loading requried librars
library(tidyverse)
library(xtable)
library(SummarizedExperiment)
library(curatedAdipoRNA)

# loading data
#data("adipo_counts")
#write_rds(adipo_counts, 'data/adipo_counts.rds')

adipo_counts <- read_rds('data/adipo_counts.rds')

# generating table
phenotype_data <- colData(adipo_counts)

phenotype_data %>%
  as.data.frame() %>%
  group_by(study_name) %>%
  summarise(pmid = unique(pmid),
            nsamples = n(),
            time = paste(unique(time), collapse = '/'),
            stages = paste(unique(stage), collapse = '/'),
            reference = paste0('\\cite{', unique(bibtexkey), '}')) %>%
  setNames(c('GEO ID', 'PMID', '(N)',
             'Time (hr)', 'Stage', 'Ref.')) %>%
  xtable(align = 'cllcccc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/datasets_rna.tex')
  
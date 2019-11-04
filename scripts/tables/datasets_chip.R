# loading requried librars
library(tidyverse)
library(xtable)
library(SummarizedExperiment)
library(ExperimentHub)
library(curatedAdipoChIP)

# loading data
# query package resources on ExperimentHub
#eh <- ExperimentHub()
#query(eh, "curatedAdipoChIP")

# load data from ExperimentHub
#peak_counts <- query(eh, "curatedAdipoChIP")[[1]]
#write_rds(peak_counts, 'data/peak_counts.rds')

peak_counts <- read_rds('data/peak_counts.rds')

# generating table
phenotype_data <- colData(peak_counts)[, c(-2, -12)]

phenotype_data %>%
  as_tibble() %>%
  filter(!is.na(control_id)) %>%
  group_by(study) %>%
  summarise(pmid = as.character(unique(pmid)),
            nsamples = n(),
            time = paste(unique(time), collapse = '/'),
            stages = paste(unique(stage), collapse = '/'),
            factor = paste(unique(factor), collapse = '/ '),
            reference = paste0('\\cite{', unique(bibtexkey), '}')) %>%
  setNames(c(c('SRA ID', 'PMID', '(N)',
               'Time (hr)', 'Stage', 'Factor', 'Ref.'))) %>%
  xtable(align = 'cllcccp{.2\\textwidth}c') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/datasets_chip.tex')

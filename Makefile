#!/bin/bash

# Define directory structure; for scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

# Define directory structure; for intermediates
DATA=data
LOG=log
LOG_FIG=log/figures
LOG_TAB=log/tables

# Define commands
RFIG=R CMD BATCH --vanilla $< $(LOG_FIG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG_TAB)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: figures tables clean

# Directories
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || exit 0
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || exit 0
dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG) && mkdir $(LOG) || exit 0
	test ! -d $(LOG_FIG) && mkdir $(LOG_FIG) || exit 0
	test ! -d $(LOG_TAB) && mkdir $(LOG_TAB) || exit 0
	
figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(FIG_DIR)/mds.png \
	$(FIG_DIR)/markers.png \
	$(FIG_DIR)/replicates_similarity.png
	
tables: ## Generate the tables
tables: dir_manuscript \
	dir_logs \
	$(TAB_DIR)/datasets_rna.tex \
	$(TAB_DIR)/datasets_chip.tex \
	$(FIG_DIR)/go_enrichment.png

# Figures
$(FIG_DIR)/mds.png: $(FIG_SRC)/mds.R \
	$(DATA)/adipo_counts.rds \
	$(DATA)/peak_counts.rds
	$(RFIG)
$(FIG_DIR)/markers.png: $(FIG_SRC)/markers.R \
	$(DATA)/adipo_counts.rds \
	$(DATA)/peak_counts.rds
	$(RFIG)
$(FIG_DIR)/go_enrichment.png: $(FIG_SRC)/go_enrichment.R \
	$(DATA)/adipo_counts.rds \
	$(DATA)/peak_counts.rds
	$(RFIG)
$(FIG_DIR)/replicates_similarity.png: $(FIG_SRC)/replicates_similarity.R \
	$(DATA)/adipo_counts.rds \
	$(DATA)/peak_counts.rds
	$(RFIG)	
	
# Tables
$(TAB_DIR)/datasets_rna.tex: $(TAB_SRC)/datasets_rna.R \
	$(DATA)/adipo_counts.rds
	$(RTAB)
$(TAB_DIR)/datasets_chip.tex: $(TAB_SRC)/datasets_chip.R \
	$(DATA)/peak_counts.rds
	$(RTAB)
	
# Clean Up
.PHONY: clean
clean: ## Clean up
clean:
	rm -f *.pdf
	rm -f *.RData

# Source: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## Print the current page
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help

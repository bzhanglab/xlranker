alias build := build-all
default: build-all

build-coloc:
  cd scripts/colocalization && \
  uv run create_colocalization_db.py && \
  cp -R output/ ../../src/xlranker/data/species

build-homologs:
  cd scripts/homologs && \
  uv run create_homolog_db.py && \
  cp -R output/ ../../src/xlranker/data/species

build-all: build-coloc build-homologs
  @echo "Built all files"

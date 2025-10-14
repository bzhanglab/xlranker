alias build := build-all
default: build-all

build-gmt:
  cd scripts/gmt && \
  uv run create_gmt.py && \
  cp -R output/ ../../src/xlranker/data/species

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

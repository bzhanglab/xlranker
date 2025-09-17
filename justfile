build-coloc:
  cd scripts/colocalization && \
  uv run create_colocalization_db.py && \
  mv coloc.pkl.gz ../../src/xlranker/data
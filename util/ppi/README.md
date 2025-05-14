# PPI Database

`create_ppi_db.py` loads `ppi.tsv` which is a two column TSV file where each row is a known PPI. The script then loads the data into Polars DataFrame such that P1 is always alphabetically before P2. Outputs a parquet file `ppi.parquet`. This file should then be moved to `src/xlranker/data`.

## Example Usage

```bash
cd util/ppi # if not already in ppi folder
python create_ppi_db.py
mv ppi.parquet ../../src/xlranker/data # moves parquet file into package
```

### Using `uv`

It is recommended to use [uv](https://github.com/astral-sh/uv) for managing the virtual environment. You can modify the script above by changed `python` to `uv run` which will make sure you have the correct packages installed.

```bash
cd util/ppi # if not already in ppi folder
uv run create_ppi_db.py # uv will handle package installation
mv ppi.parquet ../../src/xlranker/data # moves parquet file into package
```
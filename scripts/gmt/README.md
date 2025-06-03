# gmt

Python scripts to build gmt objects that are used by the machine learning pipeline. Will build database from all `gmt` files in the gmts folder.

Run `create_gmt_db.py`

```bash
cd scripts/gmt # if not already in this folder
python create_gmt_db.py
```

This creates `gmt.pkl.gz`. This needs to be placed in the `src/xlranker/data` folder.

# Colocalization

This code creates the pkl objects required for determining if two genes/proteins are located in the same cellular component. It uses data from [COMPARTMENTS](https://compartments.jensenlab.org/Search) and [Human Protein Atlas](https://www.proteinatlas.org/humanproteome/subcellular/data#locations)

## Example Usage

```bash
cd scripts/colocalization # if not already in colocalization folder
python create_colocalization_db.py
mv coloc.pkl.gz ../../src/xlranker/data # moves pkl.gz file into package
```

## Collection Dates

| XLRanker Version | COMPARTMENTS | Human Protein Atlas |
|:----------------:|:------------:|:-------------------:|
|      > 0.2.1     |  2025-09-17  |      2025-09-17     |

## Input Data

- data/human_compartment_integrated_full.tsv

This is the dataset from COMPARTMENTS that is located on the downloads page.

data/subcellular_location.tsv

This is the subcellular location information from the human protein atlas.


## References

```
Janos X. Binder, Sune Pletscher-Frankild, Kalliopi Tsafou, Christian Stolte, Seán I. O’Donoghue, Reinhard Schneider, Lars Juhl Jensen, COMPARTMENTS: unification and visualization of protein subcellular localization evidence, Database, Volume 2014, 2014, bau012, https://doi.org/10.1093/database/bau012

Mathias Uhlén et al. The human secretome.Sci. Signal.12,eaaz0274(2019).DOI:10.1126/scisignal.aaz0274
```

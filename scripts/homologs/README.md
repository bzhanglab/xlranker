# Colocalization

This code creates the pkl objects for mapping different species to its human homolog. It uses data from [Mouse Genome Informatics](https://www.informatics.jax.org/homology.shtml).

## Example Usage

```bash
cd scripts/homologs # if not already in homologs folder
python create_homoglog_db.py
```

## Collection Dates

| XLRanker Version | MGI |
|:----------------:|:------------:|
|      > 0.2.1     |  2025-09-23  |

## Input Data

- data/HOM_ALLOrganism.rpt

This is the table containing all of the homologs. It should be available at https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt.

## References

```
Baldarelli RM, Smith CL, Ringwald M, Richardson JE, Bult CJ, Mouse Genome Informatics Group. 2024. Mouse Genome Informatics: an integrated knowledgebase system for the laboratory mouse. Genetics. 2024 May 7;227(1):iyae031.
```

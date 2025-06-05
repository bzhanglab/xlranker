# Omic Data

The machine learning portion of the pipeline requires measurements that is used to prioritize ambiguous pairs. The input to the pipeline is a path to a folder containing data matrices. Traditionally, RNAseq and proteome data is provided.

## Data Format

The file should be a tab-separated file, with the first column being the protein/gene in Gene Symbol. The following columns are measurements in each sample. The value used for each gene is the mean measurement across all the included samples.

Accepted values for missing values are: `NA` and blank.

### Example

```text
idx	11BR047	11BR043
TSPAN6	23.5362298956281	21.6672110897463
TNMD	NA	NA
DPM1	26.0377091834862	25.6517806500539
SCYL3	23.5747379833942	23.9216878914892
FIRRM	18.22840670697	18.9470251111963
```

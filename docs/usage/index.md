# Usage

Install `xlranker` using `pip` or other Python package managers.

```bash
pip install xlranker
```

This will install a `xlranker` command which can be used to run the pipeline. You can also use the library if you are using a Jupyter Notebook. For notebook users, please see the [notebook example]().

## Input Data

The required input data for `xlranker` are:

[Mapping Table](input_data/mapping_table.md)
:   TSV file that maps peptide sequences to the protein matches

[Peptide Pairs](input_data/peptide_pairs.md)
:   TSV file showing all of the identified Peptide Pairs in the dataset. All sequences must appear in the mapping table.

## Running the Pipeline

??? example "Example Data"

    To test the pipeline or view the input data formatting, download the example data below
    
    [:octicons-download-16: Download example.tar.gz](#){ .md-button .md-button--secondary }


For most users, you would want to run the full pipeline. This can be achieved by running the following command:

```bash
xlranker start mapping_table.tsv peptide_pairs.tsv
```

This example assumes `mapping_table.tsv` and `peptide_pairs.tsv` are already prepared according to the instructions above and are in the current working directory.

The CLI contains multiple feature flags, such as only using the parsimony selection, saving more data, and custom filtering options. To view all of the options, please see [CLI option documentation](./CLI_options/index.md)

## Output

The output of the pipeline contains two files and a folder

[network.tsv](#)
:   TSV file of the final xl network, with all of the accepted pairs

[info.tsv](#)
:   TSV file showing all of the protein pairs with information about their prioritization status.

[plots/](#)
:   Folder containing plots showing model performance and feature importance.
# Usage

Install `xlranker` using `pip` or other Python package managers.

```bash
pip install xlranker
```

This will install a `xlranker` command which can be used to run the pipeline. You can also use the library if you are using a Jupyter Notebook. For notebook users, please see the [notebook example]().

## Input Data

The input data for `xlranker` are:

[Peptide Pairs](input_data/peptide_pairs.md)
:   TSV file showing all of the identified Peptide Pairs in the dataset.

(*Optional*) [Custom Mapping Table](input_data/custom_mapping_table.md)
:   `xlranker` uses gencode v48 to map peptide sequences to proteins. You can provide a custom peptide sequence mapping table in FASTA format or a TSV table. Please read documentation for requirements.

## Running the Pipeline

??? example "Example Data"

    To test the pipeline or view the input data formatting, download the example data below
    
    [:octicons-download-16: Download example.tar.gz](#){ .md-button .md-button--secondary }


For most users, you would want to run the full pipeline. This can be achieved by running the following command:

```bash
xlranker start peptide_pairs.tsv
```

This example assumes `peptide_pairs.tsv` is already prepared according to the instructions above and is in the current working directory.

The CLI contains multiple feature flags, such as only using the parsimony selection, saving more data, and custom filtering options. To view all of the options, please see [CLI option documentation](./CLI_options/index.md)

## Output

The output of the pipeline contains two files and a folder

[network.tsv](#)
:   TSV file of the final xl network, with all of the accepted pairs

[info.tsv](#)
:   TSV file showing all of the protein pairs with information about their prioritization status.

[plots/](#)
:   Folder containing plots showing model performance and feature importance.
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

[Omic Data](input_data/omic_data.md)
:   Omic data used by the machine learning model for prioritizing ambiguous pairs

Custom Sequence Mapping (**Strongly Recommended**, *Optional*)
:   By default, `xlranker` uses the human UNIPROT (accessed 5-30-2025) one sequence per gene to map peptide sequences to proteins. It is strongly recommended that you provide the same database used for mapping the proteomics data. You can provide either a [FASTA file](input_data/fasta.md) or a [TSV table](input_data/custom_mapping_table.md) with mapping pre-computed Please read documentation for requirements.

The typical file structure for the input looks like

```text
omic_data/
├── protein.tsv
└── rna.tsv
peptide_network.tsv
```

## Running the Pipeline

??? example "Example Data"

    To test the pipeline or view the input data formatting, download the example data below

    [:octicons-download-16: Download example.tar.gz](../downloads/example_data.tar.gz){ .md-button .md-button--secondary }


For most users, you would want to run the full pipeline. This can be achieved by running the following command:

```bash
xlranker start peptide_network.tsv omic_data/
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

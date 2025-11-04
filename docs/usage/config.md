# Config

Instead of setting multiple CLI options, you can create a config file.

You can create a config using the `xlranker init` command. This will run an interactive prompt that will create a custom config. If you just want the default configs, run `xlranker init --default`.

## CLI Configuration

If you are using the CLI you can provide a config file (JSON or YAML format) by using the `-c` option.

```
xlranker start -c xlranker_config.yml
```

## Jupyter Notebook/Library

If using a Jupyter notebook or in a Python script, you can import the config with

```py
from xlranker.config import config as xlr_config
```

You can then set various options with the `xlr_config` object.
 **Configure the config before running any other XLRanker code**.

```py
xlr_config.network = "peptide_network.tsv" # Set network path

xlr_config.advanced.save_model_files = True # Save more detailed ML related files
```

## Config Options

There are many config options, with the full description detailed at [/API/config](../API/config.md).

The key parameters that need to be set by the user include:

- **network**: the path to the [peptide network](input_data/peptide_pairs.md).
- **omic_data**: folder that contains the omic data needed to build the machine learning model. See [Omic Data](input_data/omic_data.md) for more information.
- **mapping**: Multiple peptide sequence mapping options. See [usage/mapping](mapping.md) for more details.


??? abstract "Default Config"

    Below is the default config created by running `xlranker init --default`

    ```yaml
    additional_null_values: []
    advanced:
        binary_compartments: false
        intra_in_training: false
        pair_separator: +
        save_model_files: false
    detailed: false
    fragile: false
    mapping:
        custom_table: null
        fasta_type: UNIPROT
        is_fasta: true
        reduce_fasta: false
        split_by: null
        split_index: null
    network_path: network.tsv
    omic_data_folder: omic_data/
    output: xlranker_output/
    primary_column: null
    reduce_fasta: false
    seed: null
    species: hsapiens
    threshold: 0.5
    use_homologs: false
    ```

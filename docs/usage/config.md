# Config

Instead of setting multiple CLI options, you can create a config file.

You can create a config using the `xlranker init` command. This will run an interactive prompt that will create a custom config. If you just want the default configs, run `xlranker init --default`.

## Config Options

\* = required if not set in CLI

network*
:  path to the network containing peptide sequences. Described in [input_data/peptide_pairs](input_data/peptide_pairs.md).

data_folder*
:   path to a directory containing multi-omic data used by the machine learning model. Described in [input_data/omic_data](input_data/omic_data.md).

seed
:   integer to seed random number generators. If not set, random seed selected.

custom_mapping_table
:    path to a custom mapping table (**recommended**). Can either be a FASTA file or a TSV file. Described in [input_data/fasta](input_data/fasta.md) and [input_data/custom_mapping_table](input_data/custom_mapping_table.md). If not set, uses UNIPROT human one sequence per gene acquired on May 29, 2025.

is_fasta
:   `true` if the `custom_mapping_table` is a FASTA file.

fasta_type
:   Valid options `GENCODE` or `UNIPROT`. Type of FASTA file used. _Must be set if custom_mapping_table is set_.

only_human
:   defaults to `true`. `true` if the data in the pipeline only contains human data. If `true`, allows for better negative pair generation and PPI information.

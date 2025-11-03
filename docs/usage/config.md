# Config

Instead of setting multiple CLI options, you can create a config file.

You can create a config using the `xlranker init` command. This will run an interactive prompt that will create a custom config. If you just want the default configs, run `xlranker init --default`.

## Config Options

There are many config options, with the full description detailed at [/API/config](../API/config.md).

The key parameters that need to be set by the user include:

- **network**: the path to the [peptide network](input_data/peptide_pairs.md).
- **omic_data**: folder that contains the omic data needed to build the machine learning model. See [Omic Data](input_data/omic_data.md) for more information.
- **mapping**: See [usage/mapping](mapping.md).

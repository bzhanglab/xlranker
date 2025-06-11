# CLI Options

There are many customizable parameters for the `xlranker` CLI. To view all of the options, you can run

```bash
xlranker start --help
```

```txt
Usage: xlranker start [ARGS] [OPTIONS]

Run the full prioritization pipeline

Requires input file to be in the format specified in the project documentation.

╭─ Parameters ────────────────────────────────────────────────────────────────────╮
│ NETWORK --network               -n  path to TSV file containing peptide         │
│                                     network.                                    │
│ DATA-FOLDER --data-folder       -d  folder containing the omics data for the    │
│                                     model prediction.                           │
│ CONFIG --config                 -c  if set, read and load options from config   │
│                                     file. Can be in JSON or YAML format.        │
│ SEED --seed                     -s  seed for machine learning pipeline. If not  │
│                                     set, seed is randomly selected.             │
│ VERBOSE --verbose --no-verbose  -v  enable verbose logging. [default: False]    │
│ LOG-FILE --log-file             -l  if set, saves logging to path               │
│ MAPPING-TABLE --mapping-table   -m  path to custom mapping table for peptide    │
│                                     sequences                                   │
│ SPLIT --split                       character used for splitting the FASTA file │
│                                     header                                      │
│ GS-INDEX --gs-index                 index in the FASTA file that contains the   │
│                                     gene symbol. Index starts at 0.             │
│ IS-FASTA --is-fasta                 Enable if mapping table is a FASTA file.    │
│   --no-is-fasta                     [default: False]                            │
╰─────────────────────────────────────────────────────────────────────────────────╯
```

Instead of setting these parameters each run, you can create a config file. To create the config file, run `xlranker init` or [view the documentation](../config.md).

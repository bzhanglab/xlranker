<div align="center">
<img src="images/logo.svg" alt="Logo" width="80" height="80">
<h1 style="margin-top: -1pt; margin-bottom: 0pt">XLRanker</h1>
<!-- Badges -->
<a target="_blank" style="margin-top: 0pt" href="https://colab.research.google.com/github/bzhanglab/xlranker/blob/master/notebooks/xlranker_example.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>
</div><br>

XLRanker is a tool to rank and prioritize protein to protein cross linking data. The tool uses parsimonious selection as well as an XGBoost model to select the best protein pairs to represent ambiguous peptides.

## Installation

You can install `xlranker` using pip:

```bash
pip install xlranker
```

## Usage

To analyze a TSV file containing a network of peptide sequences run:

```
xlranker start peptide_network.tsv
```

More details about the method and usage instructions can be found in [Usage](./usage/index.md)

## Example Notebook

For a quick start, you can check out the example notebook in the `notebooks` directory or [launch a google colab notebook](https://colab.research.google.com/github/bzhanglab/xlranker/blob/master/notebooks/xlranker_example.ipynb) to see how to use the package interactively.
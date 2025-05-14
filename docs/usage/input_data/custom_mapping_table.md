!!! note "For Advanced Users"
    xlranker uses GenCode v48 to map peptide sequences, which should be sufficient for most users. However, some analysis may require 

# Custom Mapping Table

## FASTA File

`xlranker` uses gencode for mapping peptide sequences to gene symbols, but you can provide your own FASTA file, the FASTA file must provide a mechanism to getting the gene symbol.

FASTA files are not necessarily consistent and `xlranker` provides simple tools to try to allow for multiple formats. You can always perform the mapping yourself, and use the [TSV table format](#tsv-table-format) as input instead of the FASTA file.

### Example

Below is an example FASTA entry from gencode v48

```fasta
>ENSP00000485175.1|ENST00000623578.3|ENSG00000169224.13|OTTHUMG00000040648.6|-|GCSAML-211|GCSAML|103
MTTFERKLQDQDKKSQEVSSTSNQENENGSGSEEVCYTVINHIPHQRSSLSSNDDGYENI
DSLTRKVRQFRERSETEYALLRTSVSRPCSCTHEHDYEVVFPH
```

`xlranker` needs to know how to extract the gene symbol from the fasta description. You can provide the split character, and the index for the gene symbol. In the above example, the split character is `|` and the index of the gene symbol is 6 (**0-based indexing**).

To test your custom mapping script, you can use the `xlranker test-fasta --split "|" --gs-index 6 mapping.fa` command. You should replace the arguments with your desired inputs. It will output the mapping for three peptide sequences using the provided parameters. Make sure the output is in gene symbol.

## TSV Table Format

The mapping table should be a tab-separated file where the first column in a line is the peptide sequence with the following columns being proteins (**GENE SYMBOL**) that map to that sequence. There is no restrictions on length, but the sequences in the mapping table must match the sequences given in the list of identified peptide pairs. 

### Example

```tsv
PEPTIDE_SEQUENCE1	PROTEIN1	PROTEIN2
PEPTIDE_SEQUENCE2	PROTEIN3
```

In the above example, `PEPTIDE_SEQUENCE1` maps to two proteins, while the second sequence is unambiguous.

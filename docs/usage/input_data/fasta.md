# FASTA File

By default, `xlranker` uses UNIPROT One Sequence Per Gene for mapping peptide sequences to gene symbols. However, it is **strongly recommended** to provide your own mapping file that matches what was used for generating the XL data. The FASTA file must provide a mechanism to getting the gene symbol. `xlranker` supports FASTA files in two formats: GENCODE and UNIPROT. If you provide your own FASTA file, provide the FASTA format with the `--fasta-type` flag in the CLI or in the config.

FASTA files are not necessarily consistent and `xlranker` provides simple tools to try to allow for multiple formats. You can always perform the mapping yourself and use the [TSV table format](custom_mapping_table.md) as input instead of the FASTA file.

### FASTA Reduction

FASTA files may contain multiple sequences for the same gene symbol. `xlranker` will reduce the FASTA file to one sequence per gene symbol by default. This is done by accepting the largest sequence only. If you want to disable this behavior, you can use the `--no-reduce-fasta` flag in the CLI or set `reduce_fasta: false` in the config.

## Example - UNIPROT Format

```fasta
>sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens OX=9606 GN=YWHAB PE=1 SV=3
MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS
WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY
LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY
YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD
AGEGEN
```

The UNIPROT FASTA format contains a `GN=` string in the description, which is set to the Gene Symbol representing a match that that FASTA entry.

## Example - GENCODE Format

Below is an example FASTA entry from GENCODE v48

```FASTA
>ENSP00000485175.1|ENST00000623578.3|ENSG00000169224.13|OTTHUMG00000040648.6|-|GCSAML-211|GCSAML|103
MTTFERKLQDQDKKSQEVSSTSNQENENGSGSEEVCYTVINHIPHQRSSLSSNDDGYENI
DSLTRKVRQFRERSETEYALLRTSVSRPCSCTHEHDYEVVFPH
```

`xlranker` needs to know how to extract the gene symbol from the fasta description. You can provide the split character, and the index for the gene symbol. In the above example, the split character is `|` and the index of the gene symbol is 6 (**0-based indexing**).

To test your custom mapping settings, you can use the `xlranker test-fasta --split "|" --gs-index 6 mapping.fa --fasta-type gencode` command. You should replace the arguments with your desired inputs. It will output the mapping for three peptide sequences using the provided parameters. Make sure the output is in gene symbol.

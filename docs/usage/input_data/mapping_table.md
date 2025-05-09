# Mapping Table

## Format

The mapping table should be a tab-separated file where the first column in a line is the peptide sequence with the following columns being proteins that map to that sequence. There is no restrictions on length, but the sequences in the mapping table must match the sequences given in the list of identified peptide pairs. 

### Example

```tsv
PEPTIDE_SEQUENCE1	PROTEIN1	PROTEIN2
PEPTIDE_SEQUENCE2	PROTEIN3
```

In the above example, `PEPTIDE_SEQUENCE1` maps to two proteins, while the second sequence is unambiguous.

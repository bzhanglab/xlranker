# TSV Table Format

!!! note "For Advanced Users"
    xlranker provides peptide sequence mapping using FASTA files. This is typically good enough for most users, but for more complex use cases, you can perform the mapping yourself and follow the below instructions.

The mapping table should be a tab-separated file where the first column in a line is the peptide sequence with the following columns being proteins (in Gene Symbol if possible) that map to that sequence. There are no restrictions on length, but the sequences in the mapping table must match the sequences given in the list of identified peptide pairs (see [Peptide Pairs](peptide_pairs.md) for more information).

### Example

```tsv
AVAWTLGVSHS 	OR4F16	OR4F29
PLLALPPQGPPG	SAMD11
```

In the above example, `AVAWTLGVSHS` maps to two proteins, while the second sequence is unambiguous.

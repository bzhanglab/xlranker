import glob
import gzip
import pickle


def process_compartments(comp_path: str = "data/human_compartment_integrated_full.tsv") -> dict[str, set[str]]:
    gene_to_comp: dict[str, set[str]] = {}
    has_header = False
    with open(comp_path) as r:
        rows = r.readlines()
    for row in rows:
        if "\t" in row:
            if has_header:
                has_header = False
                continue
            parts = row.split("\t")
            gs = parts[1]
            compartment = parts[2]
            if gs not in gene_to_comp:
                gene_to_comp.update({gs: set()})
            gene_to_comp[gs].add(compartment)
    return gene_to_comp


def process_hpa(hpa_path: str = "data/subcellular_location.tsv") -> dict[str, set[str]]:
    gene_to_comp: dict[str, set[str]] = {}
    has_header = True
    with open(hpa_path) as r:
        rows = r.readlines()
    for row in rows:
        if "\t" in row:
            if has_header:
                has_header = False
                continue
            parts = row.split("\t")
            gs = parts[1]
            main_compartments = parts[3]
            secondary_compartments = parts[4]
            if gs not in gene_to_comp:
                gene_to_comp.update({gs: set()})
            for val in main_compartments.split(";"):
                gene_to_comp[gs].add(val.upper())
            for val in secondary_compartments.split(";"):
                gene_to_comp[gs].add(val.upper())
    return gene_to_comp


with gzip.open("coloc.pkl.gz", "wb") as w:
    pickle.dump({"compartments": process_compartments(), "HumanProteinAtlas": process_hpa()}, w)

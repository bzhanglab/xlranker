"""Script that builds a pickled gzipped file that contains cell location information from COMPARTMENTS and Human Protein Atlas."""  # noqa: E501

import gzip
import pickle


def process_compartments(
    go_terms: set[str],
    comp_path: str = "data/human_compartment_integrated_full.tsv",
) -> dict[str, set[str]]:
    """Process the COMPARTMENTS database.

    Args:
        go_terms (set[str]): list of allowable go terms to filter for.
            Typically from GO Terms seen in HPA.
        comp_path (str): path to the compartments database file in TSV format

    Returns:
       dict[str, set[str]]: dictionary where keys are gene symbols and the
           values are a set of GO terms that gene is annotated to.

    """
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
            if compartment in go_terms:
                gene_to_comp[gs].add(compartment)
    return gene_to_comp


def process_hpa(
    hpa_path: str = "data/subcellular_location.tsv",
) -> tuple[set[str], dict[str, set[str]]]:
    """Process the Human Protein Atlas (HPA) database.

    Args:
        hpa_path (str): path to the compartments database file in TSV format

    Return:
       dict[str, set[str]]: dictionary where keys are gene symbols and the
           values is a set of GO terms that gene is annotated to.

    """
    gene_to_comp: dict[str, set[str]] = {}
    has_header = True
    go_terms: set[str] = set()
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
            gos = parts[13]
            if gs not in gene_to_comp:
                gene_to_comp.update({gs: set()})
            for val in main_compartments.split(";"):
                gene_to_comp[gs].add(val.upper())
            for val in secondary_compartments.split(";"):
                gene_to_comp[gs].add(val.upper())
            for val in gos.split(";"):
                # format of val is Cell Junctions (GO:0030054);Cytosol (GO:0005829);Nucleoli fibrillar center (GO:0001650) # noqa: E501
                if "GO:" in val:
                    go_id = val.split("GO:")[-1].strip(" )")
                    go_id = f"GO:{go_id}"
                    go_terms.add(go_id)
    return (go_terms, gene_to_comp)


dirs_to_create = ["output", "output/hsapiens", "output/mmusculus"]
for d in dirs_to_create:
    import os

    if not os.path.exists(d):
        os.mkdir(d)

valid_go_terms: set[str] = set()
with gzip.open("output/hsapiens/coloc.pkl.gz", "wb") as w:
    (go_terms, hpa) = process_hpa()
    valid_go_terms = go_terms
    pickle.dump(
        {
            "compartments": process_compartments(go_terms),
            "HumanProteinAtlas": hpa,
        },
        w,
    )

with gzip.open("output/mmusculus/coloc.pkl.gz", "wb") as w:
    pickle.dump(
        {
            "compartments": process_compartments(
                valid_go_terms, "data/mouse_compartment_integrated_full.tsv"
            )
        },
        w,
    )  # No HPA for mouse

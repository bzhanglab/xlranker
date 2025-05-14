import xlranker as xlr

SMALL_PROTEIN = xlr.bio.Protein(name="Small", abundance=1.0)
BIG_PROTEIN = xlr.bio.Protein(name="Big", abundance=2.0)
MISSING_A = xlr.bio.Protein(name="Missing A", abundance=None)
MISSING_B = xlr.bio.Protein(name="Missing B", abundance=None)


def test_protein_order_with_one_null():
    # small protein should always come first
    assert xlr.bio.protein.sort_proteins(SMALL_PROTEIN, MISSING_A) == (
        SMALL_PROTEIN,
        MISSING_A,
    )
    assert xlr.bio.protein.sort_proteins(MISSING_A, SMALL_PROTEIN) == (
        SMALL_PROTEIN,
        MISSING_A,
    )


def test_protein_order_with_both_null():
    # Output should be same order as input
    assert xlr.bio.protein.sort_proteins(MISSING_B, MISSING_A) == (
        MISSING_B,
        MISSING_A,
    )
    assert xlr.bio.protein.sort_proteins(MISSING_A, MISSING_B) == (
        MISSING_A,
        MISSING_B,
    )


def test_protein_order_no_nulls():
    # Big protein should always come first
    assert xlr.bio.protein.sort_proteins(SMALL_PROTEIN, BIG_PROTEIN) == (
        BIG_PROTEIN,
        SMALL_PROTEIN,
    )
    assert xlr.bio.protein.sort_proteins(BIG_PROTEIN, SMALL_PROTEIN) == (
        BIG_PROTEIN,
        SMALL_PROTEIN,
    )
    same_val_as_small = xlr.bio.Protein("Same as Small", SMALL_PROTEIN.abundance)
    assert xlr.bio.protein.sort_proteins(same_val_as_small, SMALL_PROTEIN) == (
        same_val_as_small,
        SMALL_PROTEIN,
    )
    assert xlr.bio.protein.sort_proteins(SMALL_PROTEIN, same_val_as_small) == (
        SMALL_PROTEIN,
        same_val_as_small,
    )
